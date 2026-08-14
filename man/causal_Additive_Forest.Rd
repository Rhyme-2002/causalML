\name{causal_Additive_Forest}
\alias{causal_Additive_Forest}
\title{Causal Additive Forest}

\description{
Estimates conditional average treatment effects (CATE)
and average treatment effect (ATE) using Bayesian Causal Forests.
}

\usage{
causal_Additive_Forest(X, Y, treatment, nburn = 20, nsim = 100)
}

\arguments{
\item{X}{A model matrix of confounders. The input must be a model matrix.}

\item{Y}{Binary outcome variable. Must be numeric and coded as 0 and 1, not a factor.}

\item{treatment}{Binary treatment variable. Must be numeric and coded as 0 and 1, not a factor.}

\item{nburn}{Number of burn-in iterations.}

\item{nsim}{Number of posterior simulation iterations.}
}

\value{A list containing the fitted model, CATE, and ATE.}

\examples{
\dontrun{
x1 <- rnorm(100)
x2 <- rnorm(100)
X <- model.matrix(~ x1 + x2)[, -1]
Y <- rnorm(100)
treatment <- rbinom(100, 1, 0.5)

result <- causal_Additive_Forest(X = X, Y = Y, treatment = treatment)

result$model
result$CATE
result$ATE
}
}

\keyword{causal}
