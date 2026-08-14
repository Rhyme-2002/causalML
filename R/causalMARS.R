myridge <- function(x, y, int = TRUE, lambda = 1e-6){
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  if(length(y) == 0){
    stop("No observations available.")
  }
  
  if(any(!is.finite(x))){
    stop("x contains NA/NaN/Inf.")
  }
  
  if(any(!is.finite(y))){
    stop("y contains NA/NaN/Inf.")
  }
  
  if(int){
    X <- cbind(1, x)}
  else{
    X <- x
  }
  
  p <- ncol(X)
  penalty <- diag(lambda, p)
  
  
  if(int){
    penalty[1, 1] <- 0
  }
  
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  
  coef <- tryCatch(solve(XtX + penalty, Xty), error = function(e) {
    qr.solve(XtX + penalty, Xty)
  }
  )
  
  fitted <- as.vector(X %*% coef)
  residuals <- y - fitted
  
  list(coef = as.vector(coef), fitted = fitted, res = residuals)
}

truncpow <- function(x, knot, direction = 1){
  
  if(direction == 1){
    return(pmax(x - knot, 0))
  }
  
  if(direction == 2){
    return(pmax(knot - x, 0))
  }
  stop("direction must be 1 or 2.")
}

#' Causal Multivariate Adaptive Regression Splines
#'
#' Estimates heterogeneous treatment effects using a causal
#' multivariate adaptive regression splines (MARS) algorithm.
#'
#' @param x A model matrix containing the confounder
#' variables.
#' @param treatment Binary treatment variable coded as numeric 0 and 1.
#' Must not be a factor.
#' @param y Numeric outcome variable. Must be numeric and not a factor.
#' @param maxterms Maximum number of basis functions allowed in the model.
#' @param nquant Number of quantile-based candidate knot points considered
#' for each covariate.
#' @param degree Maximum interaction degree of the basis functions.
#' @param eps Shrinkage parameter controlling the update of the residuals.
#' @param lambda Ridge penalty parameter used in the fitting procedure.
#' @param verbose Logical; if `TRUE`, displays progress messages during
#' model fitting.
#'
#' @return An object of class `causalMARS` containing:
#' \describe{
#'   \item{basis}{Final MARS basis-function matrix.}
#'   \item{parent}{Parent basis-function indices.}
#'   \item{variable}{Covariate indices used to construct new basis functions.}
#'   \item{knot}{Knot values used in the basis functions.}
#'   \item{direction}{Direction of the truncated power basis functions.}
#'   \item{basis_degree}{Degree of each basis function.}
#'   \item{basis_variables}{Variables used in each basis function.}
#'   \item{x}{Processed covariate matrix used for fitting.}
#'   \item{x_scaled}{Centered covariate matrix.}
#'   \item{y}{Processed outcome variable.}
#'   \item{treatment}{Processed binary treatment variable.}
#'   \item{x_center}{Column means used to center the covariates.}
#'   \item{variable_names}{Names of the covariates.}
#'   \item{quantiles}{Candidate knot values for each covariate.}
#'   \item{fit0}{Fitted model for the control potential outcome Y(0).}
#'   \item{fit1}{Fitted model for the treated potential outcome Y(1).}
#'   \item{y0_hat}{Estimated control potential outcomes.}
#'   \item{y1_hat}{Estimated treated potential outcomes.}
#'   \item{cate}{Estimated conditional average treatment effects (CATE).}
#'   \item{ate}{Estimated average treatment effect (ATE).}
#'   \item{cate_history}{CATE estimates recorded during model building.}
#'   \item{mse}{Factual prediction mean squared error at each step.}
#' }
#'
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' X <- model.matrix(~ x1 + x2)
#' Y <- rnorm(100)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- causalMARS(x = X, treatment = treatment, y = Y)
#'
#' result$cate
#' result$ate
#' }
#'
#' @export
causalMARS <- function(x, treatment, y, maxterms = 11, nquant = 5,
                       degree = 2, eps = 1, lambda = 1e-6, verbose = TRUE){
  tx <- treatment
  x <- as.matrix(x)
  tx <- as.numeric(tx)
  y <- as.numeric(y)
  
  if(nrow(x) != length(y)){
    stop("x and y have different numbers of observations.")
  }
  
  if(nrow(x) != length(tx)){
    stop("x and tx have different numbers of observations.")
  }
  
  if (!all(tx %in% c(0, 1))) {
    stop("tx must contain only 0 and 1.")
  }
  
  complete <- complete.cases(x, y, tx)
  if(!all(complete)){
    if (verbose) {
      cat("Removing", sum(!complete), "incomplete observations.\n")
    }
    
    x <- x[complete, , drop = FALSE]
    y <- y[complete]
    tx <- tx[complete]
  }
  
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  
  if(sum(tx == 0) < 2){
    stop("Too few control observations.")
  }
  
  if(sum(tx == 1) < 2){
    stop("Too few treated observations.")
  }
  
  
if(!is.null(colnames(x)) && "(Intercept)" %in% colnames(x)){
    intercept_col <- which(colnames(x) == "(Intercept)")
    keep <- apply(x, 2, function(z){
      length(unique(z)) > 1
    })
    keep[intercept_col] <- TRUE
    
  }else{
    keep <- apply(x, 2, function(z){
      length(unique(z)) > 1
    })
  }
  removed <- which(!keep)
  if(length(removed) > 0){
    if(verbose){
      cat("Removing", length(removed), "zero-variance variables.\n")
    }
    x <- x[, keep, drop = FALSE]
  }
  
  if(ncol(x) == 0){
    stop("No usable covariates remain.")
  }
  
  original_colnames <- colnames(x)
  if(is.null(original_colnames)){
    original_colnames <- paste0("X", seq_len(ncol(x)))
    colnames(x) <- original_colnames
  }
  
  x_center <- colMeans(x)
  x_scaled <- sweep(x, 2, x_center, "-")
  
  if(any(!is.finite(x_scaled))){
    stop("Non-finite values detected after centering.")
  }
  
  n <- nrow(x)
  p <- ncol(x)
  
  quantiles <- vector("list", p)
  
  for(j in seq_len(p)){
    vals <- sort(unique(x_scaled[, j]))
    if(length(vals) < 2){
      quantiles[[j]] <- numeric(0)
    } 
    else{probs <- seq(0, 1, length.out = nquant + 1)[-c(1, nquant + 1)]
    q <- unique(as.numeric(quantile(vals, probs = probs, na.rm = TRUE, type = 7)))
    q <- q[q > min(vals) & q < max(vals)]
    quantiles[[j]] <- q
    }
  }
  
  basis <- matrix(1, nrow = n, ncol = 1)
  colnames(basis) <- "Intercept"
  
  parent <- integer(0)
  variable <- integer(0)
  knot <- numeric(0)
  direction <- integer(0)
  basis_degree <- 0
  basis_variables <- list(integer(0))
  
  fit0 <- myridge(basis[tx == 0, , drop = FALSE], y[tx == 0], int = FALSE, lambda = lambda)
  fit1 <- myridge(basis[tx == 1, , drop = FALSE], y[tx == 1], int = FALSE, lambda = lambda)
  
  y0_hat <- as.vector(basis %*% fit0$coef)
  y1_hat <- as.vector(basis %*% fit1$coef)
  cate <- y1_hat - y0_hat
  ate <- mean(cate)
  
  residual <- numeric(n)
  residual[tx == 0] <- y[tx == 0] - mean(y[tx == 0])
  residual[tx == 1] <- y[tx == 1] - mean(y[tx == 1])
  
  cate_history <- list()
  mse_history <- numeric(0)
  cate_history[[1]] <- cate
  
  yhat_factual <- ifelse(tx == 1, y1_hat, y0_hat)
  mse_history[1] <- mean((y - yhat_factual)^2)
  
  while(ncol(basis) + 2 <= maxterms){
    if(verbose){
      cat("Current number of basis functions:", ncol(basis), "\n")
    }
    best_score <- -Inf
    best_parent <- NA_integer_
    best_variable <- NA_integer_
    best_knot <- NA_real_
    
    for(parent_col in seq_len(ncol(basis))){
      
      parent_degree <- basis_degree[parent_col]
      if(parent_degree >= degree){
        next
      }
      parent_basis <- basis[, parent_col]
      used_variables <- basis_variables[[parent_col]]
      
      for(j in seq_len(p)){
        if(j %in% used_variables){
          next
        }
        candidate_knots <- quantiles[[j]]
        if(length(candidate_knots) == 0){
          next
        }
        
        for(q in candidate_knots){
          h1 <- parent_basis * truncpow(x_scaled[, j], q, direction = 1)
          h2 <-parent_basis * truncpow(x_scaled[, j], q, direction = 2)
          if(sum(h1^2) < 1e-12 || sum(h2^2) < 1e-12){
            next
          }
          
          candidate <- cbind(h1, h2)
          fit_all <- tryCatch(myridge(candidate, residual, int = TRUE, lambda = lambda),
                              error = function(e) NULL)
          
          fit_control <- tryCatch(myridge(candidate[tx == 0, , drop = FALSE],
                                          residual[tx == 0], int = TRUE, lambda = lambda),
                                  error = function(e) NULL)
          
          fit_treated <- tryCatch(myridge(candidate[tx == 1, , drop = FALSE],
                                          residual[tx == 1], int = TRUE, lambda = lambda),
                                  error = function(e) NULL)
          
          if(is.null(fit_all) || is.null(fit_control) || is.null(fit_treated)){
            next
          }
          
          score <- sum(fit_all$res^2) - sum(fit_control$res^2) - sum(fit_treated$res^2)
          if(is.finite(score) && score > best_score){
            best_score <- score
            best_parent <- parent_col
            best_variable <- j
            best_knot <- q
          }
        }
      }
    }
    if(is.na(best_parent) || !is.finite(best_score)){
      if(verbose){cat("No further valid MARS split found.\n")}
      break
    }
    parent_basis <- basis[, best_parent]
    new1 <- parent_basis * truncpow(x_scaled[, best_variable], best_knot, direction = 1)
    new2 <- parent_basis * truncpow(x_scaled[, best_variable], best_knot, direction = 2)
    
    basis <- cbind(basis, new1, new2)
    new_col1 <- ncol(basis) - 1
    new_col2 <- ncol(basis)
    
    parent <- c(parent, best_parent, best_parent)
    variable <- c(variable, best_variable, best_variable)
    knot <- c(knot, best_knot, best_knot)
    direction <- c(direction, 1, 2)
    
    new_degree <- basis_degree[best_parent] + 1
    basis_degree <- c(basis_degree, new_degree, new_degree)
    parent_variables <- basis_variables[[best_parent]]
    new_variables <- c(parent_variables, best_variable)
    
    basis_variables[[new_col1]] <- new_variables
    basis_variables[[new_col2]] <- new_variables
    candidate <- cbind(new1, new2)
    
    fit_control_res <- myridge(candidate[tx == 0, , drop = FALSE],
                               residual[tx == 0], int = TRUE, lambda = lambda)
    fit_treated_res <- myridge(candidate[tx == 1, , drop = FALSE],
                               residual[tx == 1], int = TRUE, lambda = lambda)
    
    residual[tx == 0] <- (1 - eps) * residual[tx == 0] + eps * fit_control_res$res
    residual[tx == 1] <- (1 - eps) * residual[tx == 1] + eps * fit_treated_res$res
    
    fit0 <- myridge(basis[tx == 0, , drop = FALSE], y[tx == 0], int = FALSE, lambda = lambda)
    fit1 <- myridge(basis[tx == 1, , drop = FALSE], y[tx == 1], int = FALSE, lambda = lambda)
    
    y0_hat <- as.vector(basis %*% fit0$coef)
    y1_hat <- as.vector(basis %*% fit1$coef)
    
    cate <- y1_hat - y0_hat
    ate <- mean(cate)
    
    cate_history[[length(cate_history) + 1]] <- cate
    yhat_factual <- ifelse(tx == 1, y1_hat, y0_hat)
    mse_history <- c(mse_history, mean((y - yhat_factual)^2))
  }
  
  result <- list(basis = basis, parent = parent, variable = variable, knot = knot, 
                 direction = direction, basis_degree = basis_degree, basis_variables = basis_variables,
                 x = x, x_scaled = x_scaled, y = y, treatment = tx, x_center = x_center, 
                 variable_names = colnames(x), quantiles = quantiles, fit0 = fit0, fit1 = fit1,
                 y0_hat = y0_hat, y1_hat = y1_hat, cate = cate, ate = ate, cate_history = cate_history,
                 mse = mse_history, maxterms = maxterms, nterms = ncol(basis), nquant = nquant,
                 degree = degree, eps = eps, lambda = lambda)
  class(result) <- "causalMARS"
  return(result)
}
