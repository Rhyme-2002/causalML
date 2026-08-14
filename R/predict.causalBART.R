#' Predict from a Causal BART Model
#'
#' Generates estimated potential outcomes, conditional average treatment
#' effects (CATE), and the average treatment effect (ATE) from a fitted
#' causal BART model.
#'
#' @param object A fitted object returned by \code{\link{causalBART}}.
#' @param newdata A numeric matrix or data frame containing the same
#'   covariates used to fit the causal BART model.
#'
#' @return A list containing:
#' \describe{
#'   \item{causal_result}{
#'     A data frame containing the estimated treated potential outcome
#'     (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
#'     CATE. When the fitted model uses propensity-score adjustment,
#'     the estimated propensity score (\code{PS}) is also included.
#'   }
#'   \item{CATE}{
#'     A numeric vector containing the estimated conditional average
#'     treatment effects.
#'   }
#'   \item{ATE}{
#'     The estimated average treatment effect, calculated as the mean
#'     of the estimated CATE values.
#'   }
#' }
#'
#' @details
#' For a model fitted with \code{PS_adjusted = FALSE}, the function
#' predicts the potential outcomes by setting treatment to 0 and 1
#' while keeping the covariates fixed.
#'
#' For a model fitted with \code{PS_adjusted = TRUE}, the propensity
#' score is first predicted for \code{newdata} using the fitted
#' propensity-score BART model. The predicted propensity score is then
#' included in the outcome BART model when predicting the potential
#' outcomes.
#'
#' The CATE is calculated as
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
#' result <- causalBART(X = X, Y = Y, treatment = treatment)
#'
#' prediction <- predict.causalBART(object = result, newdata = X)
#'
#' prediction$causal_result
#' prediction$CATE
#' prediction$ATE
#' }
#'
#' @importFrom BART pbart
#' @export
predict.causalBART <- function(object, newdata){
  X <- newdata
  if(object$PS_adjusted == FALSE){
    X0 <- cbind(rep(0, nrow(X)), X)
    p0 <- predict(object = object$model, newdata = X0)
    Yhat0 <- p0$prob.test.mean
    
    X1 <- cbind(rep(1, nrow(X)), X)
    p1 <- predict(object = object$model, newdata = X1)
    Yhat1 <- p1$prob.test.mean
    
    CATE <- Yhat1 - Yhat0
    ATE <- mean(CATE)
    
    causal_result <- data.frame(Yhat1 = Yhat1, Yhat0 = Yhat0, CATE = CATE)
    return(list(causal_result = causal_result, CATE = CATE, ATE = ATE))
  }
  
  
  if(object$PS_adjusted == TRUE){
    PS_pred <- predict(object = object$PS_model, newdata = X)
    PS <- PS_pred$prob.test.mean
    
    
    X0 <- cbind(treatment = rep(0, nrow(X)), X, PS = PS)
    p0 <- predict(object = object$model, newdata = X0)
    Yhat0 <- p0$prob.test.mean
    
    X1_new <- cbind(treatment = rep(1, nrow(X)), X, PS = PS)
    p1 <- predict(object = object$model, newdata = X1_new)
    Yhat1 <- p1$prob.test.mean
    
    CATE <- Yhat1 - Yhat0
    ATE <- mean(CATE)
    
    causal_result <- data.frame(Yhat1 = Yhat1, Yhat0 = Yhat0, PS = PS, CATE = CATE)
    
    return(list(causal_result = causal_result, CATE = CATE, ATE = ATE))
  }
}
