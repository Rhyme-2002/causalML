
#' Causal Forest Estimation
#'
#' Estimates heterogeneous treatment effects using a causal forest.
#'
#' @param X A numeric matrix or data frame containing the covariates.
#' @param Y Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param treatment Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param no_of_tree Number of trees to grow in the causal forest.
#'   The default is 500.
#'
#' @return An object of class \code{causal_forest} returned by
#'   \code{\link[grf]{causal_forest}}. The returned object contains
#'   the fitted causal forest and can be used with prediction and
#'   treatment-effect estimation methods provided by the \pkg{grf}
#'   package.
#'
#' @details
#' This function fits a causal forest using the
#' \code{\link[grf]{causal_forest}} function from the \pkg{grf} package.
#' Causal forests are designed to estimate heterogeneous treatment
#' effects by allowing the treatment effect to vary across observations
#' according to their covariate values.
#'
#' The fitted model can be used to estimate conditional average
#' treatment effects (CATE) using the prediction methods provided by
#' \pkg{grf}.
#'
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' X <- model.matrix(~ x1 + x2)
#' Y <- rbinom(100, 1, 0.5)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- causalForest(
#'   X = X,
#'   Y = Y,
#'   treatment = treatment,
#'   no_of_tree = 500
#' )
#'
#' # Estimate CATE
#' cate <- predict(result)$predictions
#' }
#'
#' @importFrom grf causal_forest
#' @export
causalForest <- function(X, Y, treatment, no_of_tree = 500){
  grf::causal_forest(X = X, Y = Y, W = treatment, num.trees = no_of_tree) 
}
