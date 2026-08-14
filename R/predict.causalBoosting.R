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
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' X <- model.matrix(~ x1 + x2)
#' Y <- rnorm(100)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' result <- causal_boosting(
#'   X = X,
#'   y = Y,
#'   treatment = treatment
#' )
#'
#' prediction <- predict.causalBoosting(
#'   object = result,
#'   newdata = X
#' )
#'
#' prediction$CATE_results
#' prediction$ATE
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
