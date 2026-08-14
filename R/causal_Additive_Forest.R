#' Additive Causal Forest
#'
#' Estimates conditional average treatment effects (CATE)
#' and average treatment effect (ATE) using Bayesian Causal Forests.
#'
#' @param X A model matrix containing the confounder variables.
#' @param Y Numeric outcome variable. Must be numeric and not a factor.
#' @param treatment Binary treatment variable coded as numeric 0 and 1.
#' Must be numeric and not a factor.
#' @param nburn Number of burn-in iterations.
#' @param nsim Number of posterior simulation iterations.
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted Bayesian Causal Forest model.}
#'   \item{CATE}{Estimated conditional average treatment effects.}
#'   \item{ATE}{Estimated average treatment effect.}
#' }
#'
#' @examples
#' \dontrun{
#' x1 = rnorm(100)
#' x2 = rnorm(100) 
#' X = model.matrix(~ x1 + x2)[, -1]
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
causal_Additive_Forest <- function(X, Y, treatment, nburn = 20, nsim = 100){
  ps_model <- glm(treatment ~ X, family = binomial())
  pihat <- predict(ps_model, type = "response")
  fit <- bcf::bcf(y = Y, z = treatment, x_control = X, pihat = pihat,
                  nburn = nburn, nsim = nsim, save_tree_directory = FALSE)
  CATE <- colMeans(fit$tau)
  ATE <- mean(CATE)
  
  return(list(model = fit, CATE = CATE, ATE = ATE))
}
