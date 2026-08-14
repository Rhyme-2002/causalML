\name{parametric_standardization}
\alias{parametric_standardization}
\title{Parametric Standardization for Causal Effect Estimation}
\description{
  Estimates individual potential outcome probabilities, conditional
  average treatment effects (CATE), and the average treatment effect
  (ATE) using parametric standardization with a logistic regression model.
}
\usage{
parametric_standardization(X, Y, treatment)
}
\arguments{
  \item{X}{
    A model matrix containing the covariates used in the
    outcome model.
  }
  
  \item{Y}{
    Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
  }
  
  \item{treatment}{
    Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
  }
}
\value{
  An object of class \code{parametric_standardization} containing:
    
    \describe{
      \item{\code{model}}{
        The fitted logistic regression model used to estimate the potential
        outcome probabilities.
      }
      
      \item{\code{causal_result}}{
        A data frame containing the estimated treated potential outcome
        (\code{Yhat1}), control potential outcome (\code{Yhat0}), and
        conditional average treatment effect (\code{CATE}) for each
        observation.
      }
      
      \item{\code{CATE}}{
        A numeric vector containing the estimated conditional average
        treatment effects.
      }
      
      \item{\code{ATE}}{
        The estimated average treatment effect, calculated as the mean
        of the estimated CATE values.
      }
    }
}
\details{
  The function fits a logistic regression model for the binary outcome
  conditional on treatment and the observed covariates. The fitted model
  is then used to estimate the potential outcome probabilities under
  treatment (\code{Y(1)}) and control (\code{Y(0)}) for every observation.
  
  The conditional average treatment effect is calculated as
  
  \deqn{
    CATE(X) = P\{Y(1)=1 \mid X\} - P\{Y(0)=1 \mid X\}.
  }
  
  The average treatment effect is calculated as
  
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
    
    result <- parametric_standardization(X = X, Y = Y, treatment = treatment)
    
    result$causal_result
    result$CATE
    result$ATE
  }
}
\keyword{causal}
