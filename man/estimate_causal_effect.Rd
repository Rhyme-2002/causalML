
#' Unified Causal Effect Estimation
#'
#' Estimate causal effects using a selected causal inference method
#' available in the causalML package. The function provides a unified
#' interface for fitting causal models and extracting the estimated
#' average treatment effect (ATE) and conditional average treatment
#' effect (CATE), when available.
#'
#' @param Y Numeric outcome vector coded as 0 and 1. Must not be a factor
#'   and must contain no missing values.
#' @param X A numeric model matrix or data frame containing the
#'   covariates used for causal effect estimation. Must contain no
#'   missing values.
#' @param treatment Numeric binary treatment vector coded as 0 and 1.
#'   Both treatment groups must be present. Must not contain missing
#'   values or be a factor.
#' @param method Character string specifying the causal estimation
#'   method. Available options are \code{"causalForest"},
#'   \code{"causal_Additive_Forest"}, \code{"causal_boosting"},
#'   \code{"causalMARS"}, \code{"causalBART"}, and
#'   \code{"parametric_standardization"}. The default is
#'   \code{"causalForest"}.
#' @param ... Additional arguments passed to the selected causal
#'   estimation method.
#'
#' @return A list of class \code{"causalML_effect"} containing:
#'   \describe{
#'     \item{method}{The selected causal estimation method.}
#'     \item{fit}{The fitted model or returned model object.}
#'     \item{ATE}{The estimated average treatment effect, when available.}
#'     \item{CATE}{The estimated conditional average treatment effects,
#'       when available.}
#'   }
#'
#' @details
#' The function validates the dimensions of the outcome, treatment,
#' and covariate data and checks for missing values. It then fits the
#' selected causal model and attempts to extract the ATE and CATE
#' estimates from the returned object.
#'
#' The availability and interpretation of ATE and CATE estimates depend
#' on the selected method and its output structure. Users should inspect
#' the fitted object when estimates are not extracted automatically.
#'
#' @examples
#' \dontrun{
#' library(causalML)
#'
#' result <- estimate_causal_effect(
#'   Y = Y,
#'   X = X,
#'   treatment = treatment,
#'   method = "causalForest"
#' )
#'
#' # Estimated average treatment effect
#' result$ATE
#'
#' # Estimated conditional treatment effects
#' result$CATE
#'
#' # Fitted model object
#' result$fit
#' }
#'
#' @seealso
#' \code{\link{causalForest}},
#' \code{\link{causal_boosting}},
#' \code{\link{causalMARS}},
#' \code{\link{causalBART}},
#' \code{\link{parametric_standardization}}
#'
#' @export
