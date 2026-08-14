\name{causal_Additive_Forest}
\alias{causal_Additive_Forest}
\title{Additive Causal Forest}
\description{Estimates conditional average treatment effects (CATE) and average treatment effect (ATE) using Bayesian Causal Forests.}
\usage{causal_Additive_Forest(X, Y, treatment, nburn = 20, nsim = 100)}
\arguments{
\item{X}{A model matrix of confounders. The input must be a model matrix.}
\item{Y}{Numeric outcome variable. Must be numeric, not a factor.}
\item{treatment}{Binary treatment variable. Must be numeric and coded as 0 and 1, not a factor.}
\item{nburn}{Number of burn-in iterations.}
\item{nsim}{Number of posterior simulation iterations.}
}
\value{A list containing the fitted model, CATE, and ATE.}
\keyword{causal}
