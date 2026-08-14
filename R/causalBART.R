#' Causal BART Estimation
#'
#' Estimates heterogeneous treatment effects using Bayesian Additive
#' Regression Trees (BART).
#'
#' The function estimates the potential outcomes under treatment and
#' control and calculates the conditional average treatment effect (CATE)
#' and average treatment effect (ATE). Optionally, the propensity score
#' can be estimated using BART and included in the outcome model.
#'
#' @param X A model matrix containing the covariates.
#' @param Y Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param treatment Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param PS_adjusted Logical indicating whether propensity-score adjustment
#'   should be used. The default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{causal_result}{
#'     A data frame containing the estimated treated potential outcome
#'     (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
#'     CATE. When \code{PS_adjusted = TRUE}, the estimated propensity
#'     score (\code{PS}) is also included.
#'   }
#'   \item{CATE}{
#'     A numeric vector containing the estimated conditional average
#'     treatment effects.
#'   }
#'   \item{ATE}{
#'     The estimated average treatment effect, calculated as the mean
#'     of the estimated CATE values.
#'   }
#'   \item{PS_adjusted}{
#'     Logical value indicating whether propensity-score adjustment
#'     was used.
#'   }
#'   \item{model}{
#'     The fitted BART outcome model.
#'   }
#'   \item{PS_model}{
#'     The fitted BART propensity-score model. Returned only when
#'     \code{PS_adjusted = TRUE}.
#'   }
#' }
#'
#' @details
#' When \code{PS_adjusted = FALSE}, BART is fitted using treatment and
#' covariates as predictors of the outcome.
#'
#' When \code{PS_adjusted = TRUE}, a separate BART model is first fitted
#' to estimate the propensity score from the covariates. The estimated
#' propensity score is then included together with treatment and the
#' covariates in the outcome BART model.
#'
#' The potential outcomes are estimated by setting treatment to 0 and 1
#' for every observation. The CATE is calculated as
#' \deqn{
#' CATE(X) = E[Y(1) \mid X] - E[Y(0) \mid X].
#' }
#'
#' The ATE is calculated as
#' \deqn{
#' ATE = \frac{1}{n}\sum_{i=1}^{n} CATE(X_i).
#' }
#'
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' X <- model.matrix(~ x1 + x2)
#' Y <- rbinom(100, 1, 0.5)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' # BART without propensity-score adjustment
#' result <- causalBART(X = X, Y = Y, treatment = treatment)
#'
#' result$causal_result
#' result$CATE
#' result$ATE
#'
#' # BART with propensity-score adjustment
#' result_ps <- causalBART(X = X, Y = Y, treatment = treatment, PS_adjusted = TRUE)
#'
#' result_ps$causal_result
#' }
#'
#' @importFrom BART pbart
#' @export
causalBART <- function(X, Y, treatment, PS_adjusted = FALSE){
  if(PS_adjusted == FALSE){
    X2 <- cbind(treatment, X)
    BART0_model <- BART::pbart(x.train = X2, y.train = Y)
    
    X0 <- cbind(rep(0, nrow(X2)), X)
    p0 <- predict(object = BART0_model, newdata = X0)
    Yhat0 <- p0$prob.test.mean
    
    X1 <- cbind(rep(1, nrow(X2)), X)
    p1 <- predict(object = BART0_model, newdata = X1)
    Yhat1 <- p1$prob.test.mean
    
    CATE = Yhat1 - Yhat0
    ATE <- mean(CATE)
    causal_result <- data.frame(Yhat1 = Yhat1, Yhat0 = Yhat0, CATE = CATE)
    return(list(causal_result = causal_result, CATE = CATE, ATE = ATE,
                PS_adjusted = FALSE, model = BART0_model, class = "BART"))
  }
  if(PS_adjusted == TRUE){
    PS_model <- BART::pbart(x.train = X, y.train = treatment)
    PS <- PS_model$prob.train.mean
    
    X2 <- cbind(treatment = treatment, X, PS = PS)
    BART1_model <- BART::pbart(x.train = X2, y.train = Y)
    
    X0 <- cbind(treatment = rep(0, nrow(X)), X, PS = PS)
    p0 <- predict(object = BART1_model, newdata = X0)
    Yhat0 <- p0$prob.test.mean
    
    X1 <- cbind(treatment = rep(1, nrow(X)), X, PS = PS)
    p1 <- predict(object = BART1_model, newdata = X1)
    Yhat1 <- p1$prob.test.mean
    
    CATE <- Yhat1 - Yhat0
    ATE <- mean(CATE)
    causal_result <- data.frame(Yhat1 = Yhat1, Yhat0 = Yhat0, PS = PS, CATE = CATE)
    
    return(list(PS_model = PS_model, model = BART1_model, causal_result = causal_result,
                PS = PS, CATE = CATE, ATE = ATE, PS_adjusted = TRUE, class = "BART"))
  }
}
