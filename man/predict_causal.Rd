\name{predict_causal}
\alias{predict_causal}
\title{Predict Causal Effects from Supported Causal Models}

\description{
  Provides a unified prediction interface for several causal machine
  learning models supported by the \pkg{causalML} package.
  
  The function identifies the class of the fitted causal model and
  automatically calls the corresponding prediction method. Supported
  models include causal MARS, causal boosting, causal forest, and
  causal BART.
}

\usage{
  predict_causal(object, newdata = NULL, ...)
}

\arguments{
  \item{object}{
    A fitted causal model object. The object must have one of the
    supported classes: \code{causalMARS}, \code{causal_boosting},
    \code{causal_forest}, or \code{causal_BART}.
  }
  
  \item{newdata}{
    A numeric matrix or data frame containing the covariates for which
    causal effects are to be predicted. If \code{NULL}, the
    corresponding prediction method determines how the original data
    are handled.
  }
  
  \item{...}{
    Additional arguments passed to the corresponding prediction method.
  }
}

\value{
  The returned object depends on the class of \code{object}.
  
  \describe{
    \item{\code{causalMARS}}{
      The result returned by \code{\link{predict.causalMARS}}.
    }
    
    \item{\code{causal_boosting}}{
      The result returned by \code{\link{predict.causalBoosting}}.
    }
    
    \item{\code{causal_forest}}{
      A list containing:
        
        \describe{
          \item{\code{CATE}}{
            A data frame containing the predicted conditional average
            treatment effects.
          }
          
          \item{\code{ATE}}{
            The average treatment effect calculated as the mean of the
            predicted CATE values.
          }
        }
    }
    
    \item{\code{BART}}{
      The result returned by \code{\link{predict.causalBART}}.
    }
  }
}

\details{
  The function provides a common interface for obtaining causal effect
  predictions from different causal machine learning models.
  
  For a causal MARS model, the function calls
  \code{\link{predict.causalMARS}}.
  
  For a causal boosting model, the function calls
  \code{\link{predict.causalBoosting}}.
  
  For a causal forest model, the function uses the prediction method
  provided by the \pkg{grf} package and extracts the predicted treatment
  effects. The ATE is then calculated as the mean of the predicted CATE
  values.
  
  For a causal BART model, the function calls
  \code{\link{predict.causalBART}}.
  
  If the supplied object does not belong to one of the supported model
  classes, the function returns an error.
}

\examples{
  \dontrun{
    x1 <- rnorm(100)
    x2 <- rnorm(100)
    x3 <- rnorm(100)
    x4 <- rnorm(100)
    x5 <- rnorm(100)
    
    X <- model.matrix(~ x1 + x2 + x3 + x4 + x5)
    Y <- rbinom(100, 1, 0.5)
    treatment <- rbinom(100, 1, 0.5)
    
    model <- causalBART(X = X, Y = Y, treatment = treatment)
    
    result <- predict_causal(object = model, newdata = X)
    
    result$CATE
    result$ATE
  }
}

\keyword{causal}
\keyword{prediction}
