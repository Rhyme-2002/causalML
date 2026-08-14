\name{causalBART}
\alias{causalBART}
\title{Causal BART Estimation}
\description{
  Estimates heterogeneous treatment effects using Bayesian Additive
  Regression Trees (BART).
  
  The function estimates the potential outcomes under treatment and
  control and calculates the conditional average treatment effect (CATE)
  and average treatment effect (ATE). Optionally, the propensity score
  can be estimated using BART and included in the outcome model.
}

\usage{
causalBART(X, Y, treatment, PS_adjusted = FALSE)
}

\arguments{
  \item{X}{
    A model matrix containing the covariates.
  }
  
  \item{Y}{
    Binary outcome variable coded as numeric 0 and 1.
    Must be numeric and not a factor.
  }
  
  \item{treatment}{
    Binary treatment variable coded as numeric 0 and 1.
    Must be numeric and not a factor.
  }
  
  \item{PS_adjusted}{
    Logical indicating whether propensity-score adjustment
    should be used. The default is \code{FALSE}.
  }
}

\value{
  A list containing:
    \describe{
      \item{causal_result}{
        A data frame containing the estimated treated potential outcome
        (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
        CATE. When \code{PS_adjusted = TRUE}, the estimated propensity
        score (\code{PS}) is also included.
      }
      
      \item{CATE}{
        A numeric vector containing the estimated conditional average
        treatment effects.
      }
      
      \item{ATE}{
        The estimated average treatment effect, calculated as the mean
        of the estimated CATE values.
      }
      
      \item{PS_adjusted}{
        Logical value indicating whether propensity-score adjustment
        was used.
      }
      
      \item{model}{
        The fitted BART outcome model.
      }
      
      \item{PS_model}{
        The fitted BART propensity-score model. Returned only when
        \code{PS_adjusted = TRUE}.
      }
    }
}

\details{
  When \code{PS_adjusted = FALSE}, BART is fitted using treatment and
  covariates as predictors of the outcome.
  
  When \code{PS_adjusted = TRUE}, a separate BART model is first fitted
  to estimate the propensity score from the covariates. The estimated
  propensity score is then included together with treatment and the
  covariates in the outcome BART model.
  
  The potential outcomes are estimated by setting treatment to 0 and 1
  for every observation. The CATE is calculated as
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
    
    # BART without propensity-score adjustment
    result <- causalBART(X = X, Y = Y, treatment = treatment)
    
    result$causal_result
    result$CATE
    result$ATE
    
    # BART with propensity-score adjustment
    result_ps <- causalBART(X = X, Y = Y, treatment = treatment, PS_adjusted = TRUE)
    
    result_ps$causal_result
  }
}

\keyword{causal}
\keyword{prediction}
