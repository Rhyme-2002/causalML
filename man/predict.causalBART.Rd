\name{predict.causalBART}
\alias{predict.causalBART}
\title{Predict from a Causal BART Model}
\description{
  Generates estimated potential outcomes, conditional average treatment
  effects (CATE), and the average treatment effect (ATE) from a fitted
  causal BART model.
}
\usage{
predict.causalBART(object, newdata)
}
\arguments{
  \item{object}{
    A fitted object returned by \code{\link{causalBART}}.
  }
  
  \item{newdata}{
    A numeric matrix or data frame containing the same covariates used
    to fit the causal BART model.
  }
}
\value{
  A list containing:
    \describe{
      \item{causal_result}{
        A data frame containing the estimated treated potential outcome
        (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
        CATE. When the fitted model uses propensity-score adjustment,
        the estimated propensity score (\code{PS}) is also included.
      }
      
      \item{CATE}{
        A numeric vector containing the estimated conditional average
        treatment effects.
      }
      
      \item{ATE}{
        The estimated average treatment effect, calculated as the mean
        of the estimated CATE values.
      }
    }
}
\details{
  For a model fitted with \code{PS_adjusted = FALSE}, the function
  predicts the potential outcomes by setting treatment to 0 and 1
  while keeping the covariates fixed.
  
  For a model fitted with \code{PS_adjusted = TRUE}, the propensity
  score is first predicted for \code{newdata} using the fitted
  propensity-score BART model. The predicted propensity score is then
  included in the outcome BART model when predicting the potential
  outcomes.
  
  The CATE is calculated as
  \deqn{
    CATE(X) = E[Y(1) \mid X] - E[Y(0) \mid X].
  }
  
  The ATE is calculated as
  \deqn{
    ATE = \frac{1}{n}\sum_{i=1}^{n} CATE(X_i).
  }
}
\examples{
  \dontrun{
    x1 <- rnorm(100)
    x2 <- rnorm(100)
    X <- model.matrix(~ x1 + x2)
    Y <- rbinom(100, 1, 0.5)
    treatment <- rbinom(100, 1, 0.5)
    
    result <- causalBART(X = X, Y = Y, treatment = treatment)
    
    prediction <- predict.causalBART(object = result, newdata = X)
    
    prediction$causal_result
    prediction$CATE
    prediction$ATE
  }
}
\keyword{causal}
\keyword{prediction}
