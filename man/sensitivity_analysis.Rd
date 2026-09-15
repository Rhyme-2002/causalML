\name{sensitivity_analysis}
\alias{sensitivity_analysis}
\title{Sensitivity Analysis for an Unmeasured Binary Confounder}

\description{
Performs a quantitative sensitivity analysis for an unmeasured binary
confounder by varying its prevalence and its associations with treatment
and outcome. The causal effect is re-estimated for each combination of
sensitivity parameters using the selected causal inference method.
}

\usage{
sensitivity_analysis(
  Y,
  treatment,
  X,
  P_C = seq(0.1, 1, 0.1),
  RD_CZ = seq(-1, 1, 0.1),
  RD_CY = seq(-1, 1, 0.1),
  causal_method,
  Simulation = 1,
  ci_level = 0.95,
  digits = 3
)
}

\arguments{
  \item{Y}{
    A binary outcome vector containing only 0 and 1.
  }

  \item{treatment}{
    A binary treatment or exposure vector containing only 0 and 1.
  }

  \item{X}{
    A matrix or data frame of observed covariates used in the causal
    estimation.
  }

  \item{P_C}{
    A numeric vector specifying the assumed prevalence of the unmeasured
    binary confounder C. Values must be between 0 and 1.
    The default is \code{seq(0.1, 1, 0.1)}.
  }

  \item{RD_CZ}{
    A numeric vector specifying the assumed risk difference between the
    confounder C and treatment Z. Values range from -1 to 1.
  }

  \item{RD_CY}{
    A numeric vector specifying the assumed risk difference between the
    confounder C and outcome Y. Values range from -1 to 1.
  }

  \item{causal_method}{
    A character string specifying the causal inference method used to
    estimate the average treatment effect. Supported methods are
    \code{"causalForest"},
    \code{"causal_Additive_Forest"},
    \code{"causal_boosting"},
    \code{"causalMARS"},
    \code{"causalBART"}, and
    \code{"parametric_standardization"}.
  }

  \item{Simulation}{
    An integer specifying the number of simulations performed for each
    combination of sensitivity parameters. The default is 1.
  }

  \item{ci_level}{
    A numeric value between 0 and 1 specifying the confidence level used
    to calculate empirical confidence intervals across simulations.
    The default is 0.95.
  }

  \item{digits}{
    An integer specifying the number of decimal places used when reporting
    numerical results. The default is 3.
  }
}

\details{
The function first estimates the observed average treatment effect (ATE)
using the selected causal inference method and the observed covariates.

For each combination of the sensitivity parameters, a three-way frequency
table for the binary outcome, treatment, and hypothetical unmeasured
confounder is constructed. The confounder is then distributed among the
observations according to the specified prevalence and risk differences.

The causal effect is subsequently re-estimated after including the
simulated confounder as an additional covariate.

For a single simulation, the absolute percentage change in the ATE is
calculated as

\deqn{
100 \times \frac{|\mathrm{ATE}_{adjusted}-\mathrm{ATE}_{observed}|}
{|\mathrm{ATE}_{observed}|}
}

The directional percentage change is calculated as

\deqn{
100 \times \frac{\mathrm{ATE}_{adjusted}-\mathrm{ATE}_{observed}}
{|\mathrm{ATE}_{observed}|}
}

A positive directional percentage change indicates that adjustment for
the assumed unmeasured confounder increases the estimated ATE, whereas a
negative value indicates that it decreases the estimated ATE.

When \code{Simulation > 1}, the function repeats the confounder
generation and causal-effect estimation for each parameter combination.
The mean adjusted ATE, standard deviation, and empirical confidence
interval are then reported across simulations.

Some combinations of prevalence and risk-difference parameters may be
infeasible because they result in invalid or negative joint cell
frequencies. Such combinations are omitted from the results.

The graphical output includes text plots and heat maps showing the
adjusted ATE, absolute percentage change, and directional percentage
change across the sensitivity parameter space.
}

\value{
A list containing the following components:

\item{Observed_ATE}{
  The observed ATE estimated without the simulated unmeasured confounder.
}

\item{N_Simulations}{
  The number of simulations performed for each parameter combination.
}

\item{CI_Level}{
  The confidence level used for empirical confidence intervals when
  \code{Simulation > 1}.
}

\item{Results}{
  A data frame containing the sensitivity-analysis results for each
  feasible combination of \code{P_C}, \code{RD_CZ}, and \code{RD_CY}.
  For a single simulation, the data frame contains the adjusted ATE,
  absolute percentage change, and directional percentage change.
  For multiple simulations, it additionally contains the standard
  deviation, empirical confidence interval, and number of valid
  simulations.
}

\item{Raw_Effect_Matrix}{
  A matrix containing the estimated adjusted ATE from each simulation
  for every feasible sensitivity-parameter combination. Returned when
  \code{Simulation > 1}.
}

\item{graph_adjusted_ATE}{
  A list of text plots showing the adjusted ATE for each assumed
  confounder prevalence.
}

\item{graph_percentage_change_ATE}{
  A list of text plots showing the absolute percentage change in the ATE.
}

\item{graph_percentage_change_with_direction_ATE}{
  A list of text plots showing the directional percentage change in the
  ATE.
}

\item{graph_adjusted_ATE_heat_map}{
  A list of heat maps showing the adjusted ATE across values of
  \code{RD_CZ} and \code{RD_CY}.
}

\item{graph_percentage_change_ATE_heat_map}{
  A list of heat maps showing the absolute percentage change in the ATE.
}

\item{graph_percentage_change_with_direction_ATE_heat_map}{
  A list of heat maps showing the directional percentage change in the
  ATE.
}
}

\examples{
\dontrun{

# Single sensitivity-analysis combination
result <- sensitivity_analysis(
  Y = Y,
  treatment = treatment,
  X = X,
  P_C = 0.3,
  RD_CZ = 0.2,
  RD_CY = 0.4,
  causal_method = "causalForest",
  Simulation = 1
)

result$Observed_ATE
result$Results


# Multiple sensitivity-parameter combinations
result <- sensitivity_analysis(
  Y = Y,
  treatment = treatment,
  X = X,
  P_C = seq(0.1, 0.9, 0.1),
  RD_CZ = seq(-1, 1, 0.1),
  RD_CY = seq(-1, 1, 0.1),
  causal_method = "causalForest",
  Simulation = 10,
  ci_level = 0.95
)

head(result$Results)

}
}

\author{
Abu Sufiun Rhyme
}

\seealso{
\code{\link{causalForest}},
\code{\link{causal_Additive_Forest}},
\code{\link{causal_boosting}},
\code{\link{causalMARS}},
\code{\link{causalBART}},
\code{\link{parametric_standardization}}
}

\keyword{causal}
\keyword{sensitivity}
\keyword{confounding}
\keyword{methods}
