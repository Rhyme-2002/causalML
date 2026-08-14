#' Predict from a Causal MARS Model
#'
#' Generates potential outcome predictions, conditional average
#' treatment effects (CATE), and the average treatment effect (ATE)
#' from a fitted causal MARS model.
#'
#' @param object A fitted object of class `causalMARS`.
#' @param newdata A model matrix containing the same covariates used to
#'  fit the causal MARS model. If `NULL`, the original training
#' data are used.
#' @param type Character string specifying the type of prediction.
#' Must be one of `"CATE"`, `"mu0"`, `"mu1"`, or `"response"`.
#' The default is `"response"`.
#'
#' @return If `type = "CATE"`, a numeric vector containing the estimated
#' CATE values.
#'
#' If `type = "mu0"`, a numeric vector containing the predicted control
#' potential outcomes Y(0).
#'
#' If `type = "mu1"`, a numeric vector containing the predicted treated
#' potential outcomes Y(1).
#'
#' If `type = "response"`, a list containing:
#' \describe{
#'   \item{CATE_result}{A data frame containing predicted Y(0), Y(1),
#'   and CATE.}
#'   \item{ATE}{The average treatment effect calculated from the
#'   predicted CATE values.}
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
#' prediction <- predict.causalMARS(object = result, newdata = X)
#'
#' prediction$CATE_result
#' prediction$ATE
#'
#' cate <- predict.causalMARS(object = result, newdata = X, type = "CATE")
#' }
#'
#' @export
predict.causalMARS <- function(object, newdata = NULL, type = "response"){
predict.causalMARS <- function(object, newdata = NULL, type = "response"){
  if(!inherits(object, "causalMARS")){
    stop("object must be a causalMARS model.")
  }
  
  if(is.null(newdata)){
    newdata <- object$x
  }
  newdata <- as.matrix(newdata)
  
  if(ncol(newdata) != ncol(object$x)){
    stop("newdata must have ", ncol(object$x), " columns. ", "Received ", ncol(newdata), ".")
  }
  
  if(any(!is.finite(newdata))){
    stop("newdata contains NA/NaN/Inf.")
  }
  newdata_scaled <- sweep(newdata, 2, object$x_center, "-")
  B <- matrix(1, nrow = nrow(newdata_scaled), ncol = 1)
  n_children <- length(object$parent)
  if(n_children > 0){
    if(n_children %% 2 != 0){
      stop("Invalid MARS metadata: number of children must be even.")
    }
    nsplits <- n_children / 2
    for(s in seq_len(nsplits)){
      idx1 <- 2 * s - 1
      idx2 <- 2 * s
      parent_id <- object$parent[idx1]
      variable_id <- object$variable[idx1]
      knot_value <- object$knot[idx1]
      if(parent_id > ncol(B)){
        stop("Invalid parent index during prediction.")
      }
      
      parent_basis <- B[, parent_id]
      h1 <- parent_basis * truncpow(newdata_scaled[, variable_id], knot_value, direction = 1)
      h2 <- parent_basis * truncpow(newdata_scaled[, variable_id], knot_value, direction = 2)
      B <- cbind(B, h1, h2)
    }
  }
  
  if(ncol(B) != object$nterms){
    stop("Prediction basis has ", ncol(B), " columns but model expects ", object$nterms, ".")
  }
  
  mu0 <- as.vector(B %*% object$fit0$coef)
  mu1 <- as.vector(B %*% object$fit1$coef)
  CATE <- mu1 - mu0
  ATE <- mean(CATE)
  if(any(!is.finite(mu0)) || any(!is.finite(mu1)) || any(!is.finite(CATE))){
    stop("Prediction generated NA/NaN/Inf values.")
  }
  
  if(type == "CATE"){
    return(CATE)
  }
  
  if(type == "mu0"){
    return(mu0)
  }
  
  if(type == "mu1"){
    return(mu1)
  }
  
  if(type == "response"){return(list(CATE_result = data.frame(
    Y0_hat = mu0, Y1_hat = mu1, CATE = CATE), ATE = ATE))
  }
  
  
  stop("type must be one of: ","'CATE', 'mu0', 'mu1', or 'response'.")
}
