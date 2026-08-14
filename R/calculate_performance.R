#' Calculate Performance Measures for Causal Effect Estimation
#'
#' Calculates the root mean squared error (RMSE) for the estimated
#' conditional average treatment effects (CATE) and the absolute
#' relative bias (ARB) for the estimated average treatment effect (ATE).
#'
#' @param CATE_hat A numeric vector containing the estimated CATE values.
#' @param ATE_hat A numeric value containing the estimated average
#'   treatment effect (ATE).
#' @param TRUE_CATE A numeric vector containing the true CATE values
#'   used in the simulation study.
#' @param TRUE_ATE A numeric value containing the true average treatment
#'   effect (ATE) used in the simulation study.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{\code{CATE_RMSE}}{
#'     The root mean squared error between the estimated and true CATE
#'     values.
#'   }
#'   \item{\code{ARB}}{
#'     The absolute relative bias of the estimated ATE relative to the
#'     true ATE.
#'   }
#' }
#'
#' @details
#' The CATE root mean squared error is calculated as
#' \deqn{
#' RMSE =
#' \sqrt{
#' \frac{1}{n}\sum_{i=1}^{n}
#' \left(\widehat{\tau}_i-\tau_i\right)^2
#' }.
#' }
#'
#' The absolute relative bias of the ATE is calculated as
#' \deqn{
#' ARB =
#' \left|
#' \frac{\widehat{ATE}-TRUE\_ATE}{TRUE\_ATE}
#' \right|.
#' }
#'
#' Smaller values of both CATE RMSE and ARB indicate better estimation
#' performance.
#'
#' @examples
#' \dontrun{
#' TRUE_CATE <- c(0.10, 0.20, 0.15, 0.30, 0.25)
#' CATE_hat <- c(0.12, 0.18, 0.17, 0.28, 0.24)
#'
#' TRUE_ATE <- mean(TRUE_CATE)
#' ATE_hat <- mean(CATE_hat)
#'
#' result <- calculate_performance(CATE_hat = CATE_hat, ATE_hat = ATE_hat,
#'                                 TRUE_CATE = TRUE_CATE, TRUE_ATE = TRUE_ATE)
#'
#' result
#' }
#'
#' @export
calculate_performance <- function(CATE_hat, ATE_hat, TRUE_CATE, TRUE_ATE){
  
  CATE_RMSE <- sqrt(mean((CATE_hat - TRUE_CATE)^2))
  ARB <- abs((ATE_hat - TRUE_ATE) / TRUE_ATE)
  
  return(data.frame(CATE_RMSE = CATE_RMSE, ARB = ARB))
}
