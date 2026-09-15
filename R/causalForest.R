#' Fit a Causal Forest and Estimate Treatment Effects
#'
#' @param X A matrix or data frame of covariates.
#' @param Y A numeric vector of outcomes.
#' @param treatment A numeric vector of treatment assignments (typically
#'   binary, 0/1).
#' @param no_of_tree Integer. Number of trees to grow in the forest.
#'   Defaults to \code{500}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{model}{An object of class \code{causal_forest} returned by
#'   \code{\link[grf]{causal_forest}}. Can be used directly with
#'   prediction and treatment-effect estimation methods provided by the
#'   \pkg{grf} package.}
#'   \item{CATE}{A numeric vector of out-of-bag conditional average
#'   treatment effect (CATE) estimates, one per observation in
#'   \code{X}.}
#'   \item{ATE}{A named numeric vector with the average treatment
#'   effect (ATE) estimate and its standard error, as returned by
#'   \code{\link[grf]{average_treatment_effect}}.}
#' }
#'
#' @details
#' This function fits a causal forest using the
#' \code{\link[grf]{causal_forest}} function from the \pkg{grf} package.
#' Causal forests are designed to estimate heterogeneous treatment
#' effects by allowing the treatment effect to vary across observations
#' according to their covariate values.
#'
#' In addition to the fitted model, this function returns out-of-bag
#' conditional average treatment effect (CATE) estimates for each
#' observation, and a doubly robust estimate of the overall average
#' treatment effect (ATE) across the sample.
#'
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' X <- model.matrix(~ x1 + x2)[, -1]
#' Y <- rbinom(100, 1, 0.5)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- causalForest(X = X, Y = Y, treatment = treatment, no_of_tree = 500)
#'
#' # Fitted causal_forest object
#' result$model
#'
#' # Conditional average treatment effects (CATE)
#' result$CATE
#'
#' # Average treatment effect (ATE)
#' result$ATE
#' }
#'
#' @importFrom grf causal_forest average_treatment_effect
#' @export
causalForest <- function(X, Y, treatment, no_of_tree = 500){

  # Fit the causal forest
  model <- grf::causal_forest(X = X, Y = Y, W = treatment, num.trees = no_of_tree)

  # CATE: out-of-bag individual treatment effect estimates
  cate <- predict(model)$predictions

  # ATE: average treatment effect (doubly robust estimate)
  ate <- grf::average_treatment_effect(model, target.sample = "all")

  return(list(
    model = model,
    CATE  = cate,
    ATE   = ate
  ))
}
