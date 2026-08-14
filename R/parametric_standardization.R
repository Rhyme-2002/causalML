#' Parametric Standardization for Causal Effect Estimation
#'
#' Estimates individual potential outcome probabilities, conditional
#' average treatment effects (CATE), and the average treatment effect
#' (ATE) using parametric standardization with a logistic regression model.
#'
#' @param X A model matrix containing the covariates used in the
#'   outcome model.
#' @param Y A binary numeric outcome variable coded as 0 and 1.
#' @param treatment A binary numeric treatment variable coded as 0 and 1.
#'
#' @return An object of class \code{parametric_standardization} containing:
#' \describe{
#'   \item{\code{model}}{
#'     The fitted logistic regression model used to estimate the potential
#'     outcome probabilities.
#'   }
#'   \item{\code{causal_result}}{
#'     A data frame containing the estimated treated potential outcome
#'     (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
#'     conditional average treatment effect (\code{CATE}) for each
#'     observation.
#'   }
#'   \item{\code{CATE}}{
#'     A numeric vector containing the estimated conditional average
#'     treatment effects.
#'   }
#'   \item{\code{ATE}}{
#'     The estimated average treatment effect, calculated as the mean
#'     of the estimated CATE values.
#'   }
#' }
#'
#' @details
#' The function fits a logistic regression model for the binary outcome
#' conditional on treatment and the observed covariates. The fitted model
#' is then used to estimate the potential outcome probabilities under
#' treatment (\code{Y(1)}) and control (\code{Y(0)}) for every observation.
#'
#' The conditional average treatment effect is calculated as
#' \deqn{
#' CATE(X) = P\{Y(1)=1 \mid X\} - P\{Y(0)=1 \mid X\}.
#' }
#'
#' The average treatment effect is calculated as
#' \deqn{
#' ATE = \frac{1}{n}\sum_{i=1}^{n} CATE(X_i).
#' }
#'
#' @examples
#' \dontrun{
#' x1 = rnorm(100)
#' x2 = rnorm(100) 
#' X = model.matrix( ~ x1 + x2)
#' Y <- rnorm(100)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- parametric_standardization(X = X, Y = Y, treatment = treatment)
#'
#' result$causal_result
#' result$CATE
#' result$ATE
#' }
#'
#' @export
parametric_standardization <- function(X, Y, treatment) {
  
  X <- as.data.frame(X)
  Y <- as.numeric(Y)
  treatment <- as.numeric(treatment)
  
  if(nrow(X) != length(Y)){
    stop("X and Y have different numbers of observations.")
  }
  
  if(nrow(X) != length(treatment)){
    stop("X and treatment have different numbers of observations.")
  }
  
  if(!all(treatment %in% c(0, 1))){
    stop("Treatment must contain only 0 and 1.")
  }
  
  if(!all(Y %in% c(0, 1))){
    stop("Y must contain only 0 and 1.")
  }
  
  dat <- cbind(treatment = treatment, X)
  dat$Y <- Y
  dat <- dat[, c("Y", "treatment", names(X))]
  
  outcome_model <- glm(Y ~ .- 1, data = dat, family = binomial(link = "logit"))
  
  X0 <- X
  X0$treatment <- 0
  X1 <- X
  X1$treatment <- 1
  
  Y0_hat <- predict(outcome_model, newdata = X0, type = "response")
  Y1_hat <- predict(outcome_model, newdata = X1, type = "response")
  
  CATE <- Y1_hat - Y0_hat
  ATE <- mean(CATE)
  causal_result <- data.frame(Yhat1 = Y1_hat, Yhat0 = Y0_hat, CATE = CATE)
  
  result <- list(model = outcome_model, causal_result = causal_result, 
                 CATE = CATE, ATE = ATE)
  class(result) <- "parametric_standardization"
  return(result)
}
