split_criterion <- function(y, treatment, weights, left, right){
  tr_left <- which(treatment == 1 & left)
  co_left <- which(treatment == 0 & left)
  tr_right <- which(treatment == 1 & right)
  co_right <- which(treatment == 0 & right)
  
  if(length(tr_left) < 2 || length(co_left) < 2 || length(tr_right) < 2 || length(co_right) < 2){
    return(0)
  }
  
  weighted_mean <- function(x, w){
    sum(x * w) / sum(w)
  }
  
  weighted_variance <- function(x, w){
    if(length(x) <= 1){
      return(0)
    }
    
    m <- weighted_mean(x, w)
    sum(w * (x - m)^2) / (sum(w) - 1)
  }
  
  mean_treated_left <- weighted_mean(y[tr_left], weights[tr_left])
  
  mean_control_left <- weighted_mean(y[co_left], weights[co_left])
  
  var_treated_left <- weighted_variance(y[tr_left], weights[tr_left])
  
  var_control_left <- weighted_variance(y[co_left], weights[co_left])
  
  
  mean_treated_right <-weighted_mean(y[tr_right], weights[tr_right])
  
  mean_control_right <- weighted_mean(y[co_right], weights[co_right]) 
  
  var_treated_right <- weighted_variance(y[tr_right], weights[tr_right])
  
  var_control_right <- weighted_variance(y[co_right], weights[co_right])
  
  tau_left <- mean_treated_left - mean_control_left
  
  tau_right <- mean_treated_right - mean_control_right
  
  variance_left <- var_treated_left / length(tr_left) + var_control_left / length(co_left)
  
  variance_right <- var_treated_right / length(tr_right) + var_control_right / length(co_right)
  
  denominator <- sqrt(variance_left + variance_right)
  
  if(!is.finite(denominator) || denominator <= 1e-12){
    return(0)
  }
  
  score <- abs(tau_left - tau_right) / denominator
  
  if(!is.finite(score)){
    score <- 0
  }
  score
}



find_best_split <- function(X, y, treatment, weights, node_index, split_spread = 0.1){
  p <- ncol(X)
  best_score <- 0
  best_variable <- NA_integer_
  best_value <- NA_real_
  
  if(length(node_index) < 8){
    return(list(variable = NA_integer_, value = NA_real_, score = 0))
  }
  
  for(j in seq_len(p)){
    x <- X[node_index, j]
    if(anyNA(x)){
      next
    }
    
    unique_values <- sort(unique(x))
    
    if(length(unique_values) < 2){
      next
    }
    
    candidate_values <- (unique_values[-length(unique_values)] + unique_values[-1]) / 2
    
    if(split_spread > 0 && split_spread < 1){
      step <- max(1, floor(length(candidate_values) * split_spread))
      
      candidate_values <-candidate_values[seq(1, length(candidate_values), by = step)]
    }
    
    for(threshold in candidate_values){
      left <- rep(FALSE, nrow(X))
      right <- rep(FALSE, nrow(X))
      
      left[node_index] <- X[node_index, j] < threshold
      right[node_index] <- X[node_index, j] >= threshold
      
      n_treated_left <- sum(treatment[left] == 1)
      n_control_left <- sum(treatment[left] == 0)
      n_treated_right <- sum(treatment[right] == 1)
      n_control_right <- sum(treatment[right] == 0)
      
      if(n_treated_left < 2 || n_control_left < 2 || n_treated_right < 2 || n_control_right < 2){
        next
      }
      
      score <- split_criterion(y = y, treatment = treatment, weights = weights, left = left, right = right)
      
      if(is.finite(score) && score > best_score){
        best_score <- score
        best_variable <- j
        best_value <- threshold
      }
    }
  }
  
  list(variable = best_variable, value = best_value, score = best_score)
}


build_causal_tree <- function(X, y, treatment, weights, max_leaves = 4, split_spread = 0.1){
  n <- nrow(X)
  nodes <- list()
  nodes[[1]] <- list(left = NA_integer_, right = NA_integer_, variable = NA_integer_,
                     value = NA_real_, observations = seq_len(n), prediction0 = NA_real_,
                     prediction1 = NA_real_, is_leaf = TRUE)
  
  leaves <- 1
  next_node <- 2
  
  while(length(leaves) < max_leaves){
    best_score <- 0
    best_leaf <- NA_integer_
    best_split <- NULL
    
    for (leaf in leaves) {
      node_index <- nodes[[leaf]]$observations
      split <- find_best_split(X = X, y = y, treatment = treatment, weights = weights,
                               node_index = node_index, split_spread = split_spread)
      
      if(split$score > best_score){
        best_score <- split$score
        best_leaf <- leaf
        best_split <- split
      }
    }
    
    if(is.na(best_leaf) || is.null(best_split) || is.na(best_split$variable) || best_score <= 0){
      break
    }
    
    parent_index <- nodes[[best_leaf]]$observations
    variable <- best_split$variable
    threshold <- best_split$value
    left_index <- parent_index[X[parent_index, variable] < threshold]
    right_index <- parent_index[X[parent_index, variable] >= threshold]
    
    if(length(left_index) == 0 || length(right_index) == 0){
      break
    }
    
    left_id <- next_node
    right_id <- next_node + 1
    next_node <- next_node + 2
    nodes[[best_leaf]]$left <- left_id
    nodes[[best_leaf]]$right <- right_id
    nodes[[best_leaf]]$variable <- variable
    nodes[[best_leaf]]$value <- threshold
    nodes[[best_leaf]]$is_leaf <- FALSE
    
    nodes[[left_id]] <- list(left = NA_integer_, right = NA_integer_, variable = NA_integer_,
                             value = NA_real_, observations = left_index, prediction0 = NA_real_,
                             prediction1 = NA_real_, is_leaf = TRUE)
    
    nodes[[right_id]] <- list(left = NA_integer_, right = NA_integer_, variable = NA_integer_,
                              value = NA_real_, observations = right_index, prediction0 = NA_real_, prediction1 = NA_real_,
                              is_leaf = TRUE)
    
    leaves <- leaves[leaves != best_leaf]
    leaves <- c(leaves, left_id, right_id)
  }
  
  for (leaf in leaves) {
    index <- nodes[[leaf]]$observations
    treated <- index[treatment[index] == 1]
    control <- index[treatment[index] == 0]
    
    if(length(control) > 0){
      nodes[[leaf]]$prediction0 <- sum(weights[control] * y[control]) /sum(weights[control])
    } else {
      nodes[[leaf]]$prediction0 <- 0
    }
    if(length(treated) > 0){
      nodes[[leaf]]$prediction1 <- sum(weights[treated] * y[treated]) / sum(weights[treated])
    } else {
      nodes[[leaf]]$prediction1 <- 0
    }
  }
  nodes
}

predict_causal_tree <- function(tree, newdata){
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  pred0 <- numeric(n)
  pred1 <- numeric(n)
  for (i in seq_len(n)){
    node <- 1
    while(!tree[[node]]$is_leaf){
      variable <- tree[[node]]$variable
      threshold <- tree[[node]]$value
      if(newdata[i, variable] < threshold){
        node <- tree[[node]]$left
      } else {
        node <- tree[[node]]$right
      }
    }
    pred0[i] <- tree[[node]]$prediction0
    pred1[i] <- tree[[node]]$prediction1
  }
  list(pred0 = pred0, pred1 = pred1)
}

#' Causal Boosting
#'
#' Estimates heterogeneous treatment effects using a causal boosting algorithm.
#'
#' @param X A model matrix containing the confounder variables.
#' @param Y Numeric outcome variable. Must be numeric and not a factor.
#' @param treatment Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param n_trees Number of trees.
#' @param max_leaves Maximum number of leaves in each tree.
#' @param learning_rate Learning rate for boosting.
#' @param split_spread Proportion of candidate split points considered.
#' @param verbose Logical; whether to display progress.
#'
#' @return An object of class `causal_boosting` containing the fitted
#' model, CATE estimates, ATE, predicted potential outcomes, residuals,
#' and training MSE.
#'
#' @examples
#' \dontrun{
#' x1 = rnorm(100)
#' x2 = rnorm(100) 
#' X = model.matrix( ~ x1 + x2)
#' Y <- rnorm(100)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- causal_Additive_Forest(X = X, Y = Y,treatment = treatment)
#' result$fit
#' result$CATE
#' result$ATE
#' }
#'
#' @export
causal_boosting <- function(X, y, treatment, n_trees = 100, max_leaves = 4,
                            learning_rate = 0.01, split_spread = 0.1, verbose = TRUE){
  X <- as.matrix(X)
  y <- as.numeric(y)
  treatment <- as.numeric(treatment)
  
  if(nrow(X) != length(y)){
    stop("X and y have different number of observations.")
  }
  if(nrow(X) != length(treatment)){
    stop("X and treatment have different number of observations.")
    
  }
  if(!all(treatment %in% c(0, 1))){
    stop("Treatment must contain only 0 and 1.")
  }
  complete <- complete.cases(X, y, treatment)
  X <- X[complete, ,drop = FALSE]
  y <- y[complete]
  treatment <- treatment[complete]
  remove_cols <- integer(0)
  for(j in seq_len(ncol(X))){
    if(all(X[, j] == treatment)){
      remove_cols <- c(remove_cols, j)
    }
  }
  if(length(remove_cols) > 0){
    X <- X[ , -remove_cols, drop = FALSE]
    if(verbose){
      cat("Treatment column removed from X:", length(remove_cols), "\n")
    }
  }
  n <- nrow(X)
  weights <- rep(1, n)
  residual <- y
  G0 <- numeric(n)
  G1 <- numeric(n)
  trees <- vector("list", n_trees)
  tauhat <- matrix(NA_real_, nrow = n, ncol = n_trees)
  mse <- numeric(n_trees)
  for(k in seq_len(n_trees)){
    if(verbose && (k == 1 || k %% 10 == 0 || k == n_trees)){
      cat("Building tree:", k, "of", n_trees, "\n")
    }
    tree <- build_causal_tree(X = X, y = residual, treatment = treatment,
                              weights = weights, max_leaves = max_leaves, split_spread = split_spread)
    trees[[k]] <- tree
    fitted <- predict_causal_tree(tree = tree, newdata = X)
    fitted0 <- fitted$pred0
    fitted1 <- fitted$pred1
    G0 <- G0 + learning_rate * fitted0
    G1 <- G1 + learning_rate * fitted1
    
    control <- treatment == 0
    treated <- treatment == 1
    
    residual[control] <- residual[control] - learning_rate * fitted0[control]
    residual[treated] <- residual[treated] - learning_rate * fitted1[treated]
    
    tauhat[, k] <- G1 - G0
    yhat <- ifelse(treatment == 1, G1, G0)
    
    mse[k] <- mean((y - yhat)^2)
  }
  
  cate <- G1 - G0
  ate <- mean(cate)
  result <- list(trees = trees, cate = cate, ate = ate, y0_hat = G0, y1_hat = G1,
                 tauhat = tauhat, residuals = residual, mse = mse, n_trees = n_trees,
                 max_leaves = max_leaves, learning_rate = learning_rate, split_spread = split_spread)
  
  class(result) <- "causal_boosting"
  result
}
#' Predict from a Causal Boosting Model
#'
#' Generates potential outcome predictions, CATE estimates,
#' and the ATE from a fitted causal boosting model.
#'
#' @param object A fitted object of class `causal_boosting`.
#' @param newdata A matrix or data frame containing covariates
#' for prediction.
#'
#' @return A list containing:
#' \itemize{
#'   \item `CATE_results`: Predicted Y(0), Y(1), and CATE.
#'   \item `ATE`: Average treatment effect.
#' }
#'
#' @export
predict.causalBoosting <- function(object, newdata){
  if(!inherits(object, "causal_boosting")){
    stop("object must be a causal_boosting model.")
  }
  newdata <- as.matrix(newdata)
  n_new <- nrow(newdata)
  G0 <- numeric(n_new)
  G1 <- numeric(n_new)
  for (tree in object$trees){
    prediction <- predict_causal_tree(tree = tree, newdata = newdata)
    G0 <- G0 + object$learning_rate * prediction$pred0
    G1 <- G1 + object$learning_rate * prediction$pred1
  }
  CATE <- G1 - G0
  CATE_results <- data.frame(Y0_hat = G0, Y1_hat = G1, CATE = CATE)
  ATE <- mean(CATE)
  return(list(CATE_results = CATE_results, ATE = ATE))
}
