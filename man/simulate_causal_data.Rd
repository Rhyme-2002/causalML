\name{simulate_causal_data}
\alias{simulate_causal_data}
\title{Simulate Data for Causal Machine-Learning Methods}
\description{
Generates potential outcome probabilities and a simulated binary
outcome for evaluating causal machine-learning methods. The
simulation framework follows the general approach of Wendling et al.
(2018), adapted to estimate heterogeneous treatment effects for a
binary treatment and binary outcome.
}
\usage{
simulate_causal_data(Y, X, treatment, model = "logit", k = 5,
                     OR = NULL, seed = NULL)
}
\arguments{
  \item{Y}{
    Binary outcome variable coded as numeric 0 and 1. Must be numeric
    and not a factor.
  }

  \item{X}{
    A model matrix or data frame containing the observed covariates.
  }

  \item{treatment}{
    Binary treatment variable coded as numeric 0 and 1. Must be numeric
    and not a factor.
  }

  \item{model}{
    Character string specifying the model used to estimate the
    conditional outcome probabilities. Available methods are
    \code{"logit"}, \code{"BART"}, \code{"GBM"},
    \code{"Random_Forest"}, and \code{"KNN"}. The default is
    \code{"logit"}.
  }

  \item{k}{
    Number of nearest neighbors used when
    \code{model = "KNN"}. The default is 5.
  }

  \item{OR}{
    The prespecified odds ratio used to calibrate the treatment effect.
    If \code{NULL}, the odds ratio is estimated from a logistic
    regression model fitted to the observed outcome, treatment, and
    covariates.
  }

  \item{seed}{
    Optional integer used to set the random-number seed.
  }
}
\value{
A list containing:

\describe{
  \item{simulated_data}{
    A data frame containing the estimated untreated probability
    (\code{pi0}), the original treated probability (\code{pi1}),
    the adjusted treated probability (\code{pi1_adjusted}),
    the factual outcome probability (\code{pi_factual}), and the
    true conditional treatment effect (\code{tau_x}).
  }

  \item{X}{
    The covariate matrix used in the simulation.
  }

  \item{treatment}{
    The binary treatment variable.
  }

  \item{pi1}{
    Estimated conditional outcome probability under treatment
    before odds-ratio calibration.
  }

  \item{pi0}{
    Estimated conditional outcome probability under control.
  }

  \item{pi1_adjusted}{
    Adjusted conditional outcome probability under treatment after
    calibration using the specified odds ratio.
  }

  \item{pi_factual}{
    The conditional probability of the observed factual outcome,
    determined by the observed treatment assignment.
  }

  \item{TRUE_CATE}{
    A numeric vector containing the true conditional average
    treatment effects.
  }

  \item{TRUE_ATE}{
    The true average treatment effect, calculated as the mean of
    the true CATE values.
  }

  \item{model}{
    The fitted model used to estimate the conditional outcome
    probabilities.
  }

  \item{method}{
    The name of the model used for simulation.
  }

  \item{OR}{
    The odds ratio used to calibrate the treatment effect.
  }
}
}
\details{
This function implements a simulation framework for evaluating
heterogeneous treatment-effect estimators. The framework is
adapted from Wendling et al. (2018), who considered simulation
settings with binary treatment and binary outcomes and used the
conditional risk difference as the treatment-effect estimand.

First, the conditional outcome probabilities under treatment and
control are estimated from the observed data using one of five
models: logistic regression, Bayesian Additive Regression Trees
(BART), gradient boosting machines (GBM), random forest, or
K-nearest neighbors (KNN).

Let \eqn{X} denote the observed covariates and \eqn{Z} the binary
treatment indicator. The fitted model is used to obtain

\deqn{
\pi^1(x) = P(Y = 1 \mid Z = 1, X = x)
}

and

\deqn{
\pi^0(x) = P(Y = 1 \mid Z = 0, X = x).
}

The estimated potential outcome probabilities are then used to
construct the treatment-effect function on the logit scale.
Specifically,

\deqn{
\gamma(x) =
\operatorname{logit}\{\pi^1(x)\}
-
\operatorname{logit}\{\pi^0(x)\}
+ c.
}

Here, \eqn{c} is a calibration constant. The constant is selected
so that the mean of the treatment-effect function on the logit
scale corresponds to the logarithm of the specified odds ratio:

\deqn{
\frac{1}{n}\sum_{i=1}^{n}\gamma(X_i) = \log(OR).
}

The adjusted treated potential outcome probability is then
obtained using the inverse-logit transformation:

\deqn{
\pi^1_{\mathrm{adjusted}}(x) =
\operatorname{logit}^{-1}
\left[
\operatorname{logit}\{\pi^0(x)\} + \gamma(x)
\right].
}

The conditional treatment effect is defined on the probability
scale as the conditional risk difference:

\deqn{
\tau(x) =
\pi^1_{\mathrm{adjusted}}(x) - \pi^0(x).
}

The true average treatment effect is calculated as

\deqn{
ATE =
\frac{1}{n}\sum_{i=1}^{n}\tau(X_i).
}

The factual outcome probability is determined by the observed
treatment assignment:

\deqn{
\pi_{\mathrm{factual},i} =
Z_i\pi^1_{\mathrm{adjusted}}(X_i)
+
(1-Z_i)\pi^0(X_i).
}

The binary outcome is then generated from a Bernoulli distribution:

\deqn{
Y_i \sim
\operatorname{Bernoulli}
\left(\pi_{\mathrm{factual},i}\right).
}

Therefore, the simulation provides known values of the true CATE
and ATE, which can be used to evaluate the performance of causal
machine-learning estimators using measures such as RMSE and
absolute relative bias.

If \code{OR = NULL}, the odds ratio is estimated from a logistic
regression model using the observed outcome, treatment, and
covariates.

Small numerical bounds are applied to the estimated probabilities
before calculating the logit transformation to avoid numerical
problems when probabilities are exactly 0 or 1.
}
\references{
Wendling, T., Jung, K., Callahan, A., Schuler, A., Shah, N. H.,
& Gallego, B. (2018).
Comparing methods for estimation of heterogeneous treatment effects
using observational data from health care databases.
\emph{Statistics in Medicine}, 37(23), 3309--3324.
\doi{10.1002/sim.7820}.
}
\examples{
\dontrun{
set.seed(123)

X <- matrix(rnorm(500), nrow = 100, ncol = 5)
treatment <- rbinom(100, 1, 0.5)
Y <- rbinom(100, 1, 0.5)

# Logistic regression
result <- simulate_causal_data(Y = Y, X = X, treatment = treatment, model = "logit", OR = 1.5, seed = 123)

result$simulated_data
result$TRUE_CATE
result$TRUE_ATE

# BART
result_bart <- simulate_causal_data(Y = Y, X = X, treatment = treatment, model = "BART", OR = 1.5, seed = 123)
}
}
\keyword{simulation}
\keyword{causal}
\keyword{causal inference}
\keyword{machine learning}
