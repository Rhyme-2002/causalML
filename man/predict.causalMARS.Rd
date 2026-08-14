\name{predict.causalMARS}
\alias{predict.causalMARS}
\title{Predict from a Causal MARS Model}
\description{
Generates potential outcome predictions, conditional average
treatment effects (CATE), and the average treatment effect (ATE)
from a fitted causal MARS model.
}
\usage{
causalMARS(object, newdata = NULL, type = "response")
}
\arguments{
\item{object}{
A fitted object of class \code{causalMARS}.
}

\item{newdata}{
A numeric matrix containing the same covariates used to fit
the causal MARS model. If \code{NULL}, the original training
data are used.
}

\item{type}{
Character string specifying the type of prediction.
Must be one of \code{"CATE"}, \code{"mu0"}, \code{"mu1"},
or \code{"response"}. The default is \code{"response"}.
}
}
\value{
A value determined by the \code{type} argument.

\describe{
\item{\code{type = "CATE"}}{
A numeric vector containing the estimated conditional
average treatment effects (CATE).
}

\item{\code{type = "mu0"}}{
A numeric vector containing the predicted control potential
outcomes Y(0).
}

\item{\code{type = "mu1"}}{
A numeric vector containing the predicted treated potential
outcomes Y(1).
}

\item{\code{type = "response"}}{
A list containing:

\describe{
\item{\code{CATE_result}}{
A data frame containing predicted Y(0), Y(1), and CATE.
}

\item{\code{ATE}}{
The average treatment effect calculated from the
predicted CATE values.
}
}
}
}
}
\examples{
\dontrun{
x1 <- rnorm(100)
x2 <- rnorm(100)
X <- model.matrix(~ x1 + x2)
Y <- rnorm(100)
treatment <- rbinom(100, 1, 0.5)

result <- causalMARS(x = X, treatment = treatment, y = Y)

prediction <- predict.causalMARS(object = result, newdata = X)

prediction$CATE_result
prediction$ATE

cate <- predict.causalMARS(object = result, newdata = X, type = "CATE")

mu0 <- predict.causalMARS(object = result, newdata = X, type = "mu0")
mu1 <- predict.causalMARS(object = result, newdata = X, type = "mu1")
}
}
\keyword{causal}
\keyword{prediction}
