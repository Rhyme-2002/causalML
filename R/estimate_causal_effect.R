
#' Unified Causal Effect Estimation
#'
#' Estimate causal effects using a selected method from causalML.
#'
#' @param Y Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param X A model matrix containing the covariates.
#' @param treatment Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param compare_model Character vector specifying the causal models
#'   to be compared. Available models are \code{"causalMARS"},
#'   \code{"causal_boosting"}, \code{"causal_forest"},
#'   \code{"causal_BART"}, and \code{"Bayesian_Causal_Forest"}.
#'   The default compares all available causal machine learning models.
#' @param ... Additional arguments passed to the selected method.
#'
#' @return A list containing the selected method, fitted model,
#'   estimated ATE, and estimated CATE when available.
#' @export

estimate_causal_effect <- function(
    Y,
    X,
    treatment,
    method = c(
      "causalForest",
      "causal_Additive_Forest",
      "causal_boosting",
      "causalMARS",
      "causalBART",
      "parametric_standardization"
    ),
    ...
) {

  method <- match.arg(method)
  extra_args <- list(...)

  # Validate inputs
  if (length(Y) != length(treatment) ||
      nrow(as.data.frame(X)) != length(Y)) {
    stop("Y, X, and treatment must have the same number of observations.")
  }

  if (anyNA(Y) || anyNA(treatment) || anyNA(X)) {
    stop("Missing values are not allowed. Handle them before fitting.")
  }

  if (!all(unique(treatment) %in% c(0, 1)) ||
      length(unique(treatment)) != 2) {
    stop("treatment must contain both 0 and 1.")
  }

  # Fit the selected causal model
  fit <- switch(
    method,

    causalForest = do.call(
      causalML::causalForest,
      c(list(Y = Y, X = X, treatment = treatment), extra_args)
    ),

    causal_Additive_Forest = do.call(
      causalML::causal_Additive_Forest,
      c(list(X = X, Y = Y, treatment = treatment), extra_args)
    ),

    causal_boosting = do.call(
      causalML::causal_boosting,
      c(list(X = X, y = Y, treatment = treatment), extra_args)
    ),

    causalMARS = do.call(
      causalML::causalMARS,
      c(list(x = X, y = Y, treatment = treatment), extra_args)
    ),

    causalBART = do.call(
      causalML::causalBART,
      c(list(X = X, Y = Y, treatment = treatment), extra_args)
    ),

    parametric_standardization = do.call(
      causalML::parametric_standardization,
      c(list(X = X, Y = Y, treatment = treatment), extra_args)
    )
  )

  # Obtain predictions for causal forest
  if (method == "causalForest") {
    fit <- causalML::predict_causal(
      object = fit,
      newdata = X
    )
  }

  # Extract estimates using supported output field names
  get_estimate <- function(object, candidates) {
    if (is.null(object)) return(NULL)

    for (nm in candidates) {
      value <- object[[nm]]
      if (!is.null(value)) return(value)
    }

    NULL
  }

  ate <- get_estimate(
    fit,
    c("ATE", "ate", "ATE_estimate", "estimate")
  )

  cate <- get_estimate(
    fit,
    c("CATE", "cate", "CATE_estimate")
  )

  if (is.null(ate) && !is.null(cate)) {
    ate <- mean(cate, na.rm = TRUE)
  }

  if (is.null(ate)) {
    warning(
      "ATE was not found automatically. Inspect the fitted model in $fit."
    )
  }

  result <- list(
    method = method,
    fit = fit,
    ATE = ate,
    CATE = cate
  )

  class(result) <- "causalML_effect"
  return(result)
}
