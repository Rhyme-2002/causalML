\name{causalMARS}
\alias{causalMARS}
\title{Causal Multivariate Adaptive Regression Splines}
\description{
  Estimates heterogeneous treatment effects using a causal
  multivariate adaptive regression splines (MARS) algorithm.
}
\usage{
causalMARS(x, treatment, y, maxterms = 11, nquant = 5,
           degree = 2, eps = 1, lambda = 1e-6, verbose = TRUE)
}
\arguments{
  \item{x}{
    A model matrix containing the confounder variables.
  }
  
  \item{treatment}{
    Binary treatment variable coded as numeric 0 and 1. Must not be a factor.
  }
  
  \item{y}{
    Binary outcome variable coded as numeric 0 and 1. Must not be a factor.
  }
  
  \item{maxterms}{
    Maximum number of basis functions allowed in the model.
  }
  
  \item{nquant}{
    Number of quantile-based candidate knot points considered
    for each covariate.
  }
  
  \item{degree}{
    Maximum interaction degree of the basis functions.
  }
  
  \item{eps}{
    Shrinkage parameter controlling the update of the residuals.
  }
  
  \item{lambda}{
    Ridge penalty parameter used in the fitting procedure.
  }
  
  \item{verbose}{
    Logical; if \code{TRUE}, displays progress messages during
    model fitting.
  }
}
\value{
  An object of class \code{causalMARS} containing the fitted
  causal MARS model and treatment effect estimates.
  
  \describe{
    \item{basis}{
      Final MARS basis-function matrix.
    }
    
    \item{parent}{
      Parent basis-function indices.
    }
    
    \item{variable}{
      Covariate indices used to construct new basis functions.
    }
    
    \item{knot}{
      Knot values used in the basis functions.
    }
    
    \item{direction}{
      Direction of the truncated power basis functions.
    }
    
    \item{basis_degree}{
      Degree of each basis function.
    }
    
    \item{basis_variables}{
      Variables used in each basis function.
    }
    
    \item{x}{
      Processed covariate matrix used for fitting.
    }
    
    \item{x_scaled}{
      Centered covariate matrix.
    }
    
    \item{y}{
      Processed outcome variable.
    }
    
    \item{treatment}{
      Processed binary treatment variable.
    }
    
    \item{x_center}{
      Column means used to center the covariates.
    }
    
    \item{variable_names}{
      Names of the covariates.
    }
    
    \item{quantiles}{
      Candidate knot values for each covariate.
    }
    
    \item{fit0}{
      Fitted model for the control potential outcome Y(0).
    }
    
    \item{fit1}{
      Fitted model for the treated potential outcome Y(1).
    }
    
    \item{y0_hat}{
      Estimated control potential outcomes.
    }
    
    \item{y1_hat}{
      Estimated treated potential outcomes.
    }
    
    \item{cate}{
      Estimated conditional average treatment effects (CATE).
    }
    
    \item{ate}{
      Estimated average treatment effect (ATE).
    }
    
    \item{cate_history}{
      CATE estimates recorded during model building.
    }
    
    \item{mse}{
      Factual prediction mean squared error at each step.
    }
    
    \item{nterms}{
      Number of basis functions in the final model.
    }
  }
}
\examples{
  \dontrun{
    x1 <- rnorm(100)
    x2 <- rnorm(100)
    X <- model.matrix(~ x1 + x2)[, -1]
    Y <- rnorm(100)
    treatment <- rbinom(100, 1, 0.5)
    
    result <- causalMARS(x = X, treatment = treatment, y = Y)
    
    result$cate
    result$ate
  }
}
\keyword{causal}
\keyword{models}
\keyword{nonparametric}
