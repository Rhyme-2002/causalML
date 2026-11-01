\name{estimate_causal_effect}
\alias{estimate_causal_effect}
\title{Unified Causal Effect Estimation}
\description{
  Estimates causal effects using a selected causal inference method
  available in the causalML package. Provides a unified interface for
  fitting causal models and extracting the estimated average treatment
  effect (ATE) and conditional average treatment effect (CATE), when
  available.
}
\usage{
estimate_causal_effect(Y, X, treatment,
                       method = "causalForest", ...)
}
\arguments{
  \item{Y}{
    Numeric binary outcome vector coded as 0 and 1. Must not be a
    factor and must contain no missing values.
  }

  \item{X}{
    A numeric model matrix or data frame containing the covariates
    used for causal effect estimation. Must contain no missing values.
  }

  \item{treatment}{
    Numeric binary treatment vector coded as 0 and 1. Both treatment
    groups must be present. Must not be a factor or contain missing
    values.
  }

  \item{method}{
    Character string specifying the causal estimation method.
    Available options are \code{"causalForest"},
    \code{"causal_Additive_Forest"}, \code{"causal_boosting"},
    \code{"causalMARS"}, \code{"causalBART"}, and
    \code{"parametric_standardization"}.
    The default is \code{"causalForest"}.
  }

  \item{...}{
    Additional arguments passed to the selected causal estimation
    method.
  }
}
\value{
  A list of class \code{causalML_effect} containing:

  \describe{
    \item{method}{
      The selected causal estimation method.
    }

    \item{fit}{
      The fitted model or returned model object.
    }

    \item{ATE}{
      The estimated average treatment effect (ATE), when available.
    }

    \item{CATE}{
      The estimated conditional average treatment effects (CATE),
      when available.
    }
  }
}
\details{
  The function validates the dimensions of the outcome, treatment,
  and covariate data and checks for missing values. It then fits the
  selected causal model and attempts to extract the ATE and CATE
  estimates from the returned object.

  The availability and interpretation of ATE and CATE estimates
  depend on the selected method and its output structure. Users
  should inspect the fitted object when estimates are not extracted
  automatically.
}
\examples{
  \dontrun{
    library(causalML)

    result <- estimate_causal_effect(
      Y = Y,
      X = X,
      treatment = treatment,
      method = "causalForest"
    )

    # Estimated average treatment effect
    result$ATE

    # Estimated conditional treatment effects
    result$CATE

    # Fitted model object
    result$fit
  }
}
\seealso{
  \code{\link{causalForest}},
  \code{\link{causal_Additive_Forest}},
  \code{\link{causal_boosting}},
  \code{\link{causalMARS}},
  \code{\link{causalBART}},
  \code{\link{parametric_standardization}}
}
\keyword{causal}
\keyword{models}
\keyword{nonparametric}

