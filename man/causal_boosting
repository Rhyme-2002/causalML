\name{causal_boosting}
\alias{causal_boosting}
\title{Causal Boosting}

\description{
  Estimates conditional average treatment effects (CATE)
  and average treatment effect (ATE) using a causal boosting
  algorithm based on causal trees.
}

\usage{
  causal_boosting(X, y, treatment, n_trees = 100, max_leaves = 4,
                  learning_rate = 0.01, split_spread = 0.1,
                  verbose = TRUE)
}

\arguments{
  \item{X}{A numeric matrix of covariates.}
  
  \item{y}{Numeric outcome variable. Must be numeric, not a factor.}
  
  \item{treatment}{Binary treatment variable. Must be numeric and
    coded as 0 and 1, not a factor.}
  
  \item{n_trees}{Number of causal trees to be fitted.}
  
  \item{max_leaves}{Maximum number of terminal leaves allowed in
    each causal tree.}
  
  \item{learning_rate}{Learning rate used to update the fitted
    treatment-specific predictions at each boosting iteration.}
  
  \item{split_spread}{Proportion used to control the number of
    candidate split points evaluated when constructing each tree.}
  
  \item{verbose}{Logical value indicating whether progress messages
    should be displayed during model fitting.}
}

\value{
  A list containing:
    \itemize{
      \item \code{trees}: A list of fitted causal trees.
      \item \code{cate}: Estimated conditional average treatment effects.
      \item \code{ate}: Estimated average treatment effect.
      \item \code{y0_hat}: Estimated potential outcomes under control.
      \item \code{y1_hat}: Estimated potential outcomes under treatment.
      \item \code{tauhat}: CATE estimates from each boosting iteration.
      \item \code{residuals}: Final model residuals.
      \item \code{mse}: Mean squared error for each boosting iteration.
      \item \code{n_trees}: Number of trees used.
      \item \code{max_leaves}: Maximum number of leaves per tree.
      \item \code{learning_rate}: Learning rate used for boosting.
      \item \code{split_spread}: Split spread parameter used during
      tree construction.
    }
}

\examples{
  \dontrun{
    x1 <- rnorm(100)
    x2 <- rnorm(100)
    X <- model.matrix(~ x1 + x2)
    Y <- rnorm(100)
    treatment <- rbinom(100, 1, 0.5)
    
    result <- causal_boosting(
      X = X,
      y = Y,
      treatment = treatment
    )
    
    result$cate
    result$ate
  }
}

\keyword{causal}
