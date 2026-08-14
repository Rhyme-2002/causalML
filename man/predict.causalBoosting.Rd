\name{predict.causalBoosting}
\alias{predict.causalBoosting}
\title{Predict from a Causal Boosting Model}

\description{
  Generates potential outcome predictions, conditional average
  treatment effect (CATE) estimates, and the average treatment
  effect (ATE) from a fitted causal boosting model.
}

\usage{
predict.causalBoosting(object, newdata)
}

\arguments{
  \item{object}{A fitted object of class \code{causal_boosting}.}
  
  \item{newdata}{A model matrix containing covariates
    for prediction.}
}

\value{
  A list containing:
    \itemize{
      \item \code{CATE_results}: A data frame containing estimated
      potential outcomes under control (\code{Y0_hat}), treatment
      (\code{Y1_hat}), and CATE.
      \item \code{ATE}: Estimated average treatment effect.
    }
}

\examples{
  \dontrun{
    x1 <- rnorm(100)
    x2 <- rnorm(100)
    X <- model.matrix(~ x1 + x2)
    Y <- rnorm(100)
    treatment <- rbinom(100, 1, 0.5)

    result <- causal_boosting(X = X, y = Y, treatment = treatment)
    prediction <- predict.causalBoosting(object = result, newdata = X)
    
    prediction$CATE_results
    prediction$ATE
  }
}

\keyword{causal}
