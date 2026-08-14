\name{compare_model}
\alias{compare_model}
\title{Comparison of Causal Machine Learning Models}
\description{
  Compares several causal inference and causal machine learning methods
  using a simulation-based study. The function generates potential outcome
  probabilities using a specified data-generating model, repeatedly
  simulates binary outcomes, fits the selected causal models, and evaluates
  their ability to estimate the conditional average treatment effect (CATE)
  and average treatment effect (ATE).
  
  The simulation procedure follows the general framework proposed by
  Wendling et al. (2018).
}
\usage{
compare_model(Y, X, treatment, model_data_simulation = "logit", sim = 100,
    compare_model = c("causalMARS", "causal_boosting", "causal_forest", "causal_BART", "causal_Additive_Forest"),
    ...)
}
\arguments{
  \item{Y}{
    A binary numeric outcome variable coded as 0 and 1. This outcome
    is used to estimate the data-generating mechanism.
  }
  
  \item{X}{
    A numeric matrix or data frame containing the observed covariates.
    The number of rows must be equal to the length of \code{Y} and
    \code{treatment}.
  }
  
  \item{treatment}{
    A binary numeric treatment variable coded as 0 and 1.
  }
  
  \item{model_data_simulation}{
    Character string specifying the model used to estimate the
    potential outcome probabilities in the data-generating process.
    Available options are \code{"logit"}, \code{"BART"}, \code{"GBM"},
    \code{"Random_Forest"}, and \code{"KNN"}.
    The default is \code{"logit"}.
  }
  
  \item{sim}{
    Number of simulation replications. The default is 100.
  }
  
  \item{compare_model}{
    Character vector specifying the causal models to be compared.
    Available options are \code{"causalMARS"},
    \code{"causal_boosting"}, \code{"causal_forest"},
    \code{"causal_BART"}, and \code{"causal_Additive_Forest"}.
    The default compares all five methods.
    
    Parametric standardization is automatically included as a
    benchmark method.
  }
  
  \item{...}{
    Additional arguments passed to the underlying causal modeling
    functions.
  }
}
\value{
  A list containing the following components:
    
    \describe{
      \item{\code{RMSE_df}}{
        A data frame containing the CATE RMSE for every simulation
        replication and each compared model.
      }
      
      \item{\code{RMSE_long}}{
        A long-format data frame containing the CATE RMSE values,
        with columns for simulation number, model, and RMSE.
      }
      
      \item{\code{ARB_df}}{
        A data frame containing the absolute relative bias (ARB) of
        the ATE for every simulation replication and each compared model.
      }
      
      \item{\code{ARB_long}}{
        A long-format data frame containing the ARB values,
        with columns for simulation number, model, and ARB.
      }
      
      \item{\code{RMSE_plot}}{
        A \code{ggplot2} boxplot showing the distribution of CATE RMSE
        across simulation replications for each model.
      }
      
      \item{\code{ARB_plot}}{
        A \code{ggplot2} boxplot showing the distribution of ATE ARB
        across simulation replications for each model.
      }
      
      \item{\code{percentage_of_exposure_per_iteration}}{
        The percentage of observations receiving the simulated outcome
        equal to 1 in the final simulation iteration.
      }
    }
}
\details{
  \subsection{Simulation procedure}{
    The simulation study follows the general data-generating framework
    proposed by Wendling et al. (2018).
    
    First, the observed covariates \code{X}, binary outcome \code{Y}, and
    binary treatment variable are supplied to the function. The function
    calls \code{\link{simulate_causal_data}} to estimate the conditional
    potential outcome probabilities.
    
    The potential outcome probabilities are estimated using the model
    specified by \code{model_data_simulation}. The available models are
    logistic regression, BART, gradient boosting, random forest, and
    K-nearest neighbors.
    
    For every observation, the fitted model is used to estimate
    
    \deqn{
      \pi_1(X) = P\{Y(1)=1 \mid X\}
    }
    
    and
    
    \deqn{
      \pi_0(X) = P\{Y(0)=1 \mid X\}.
    }
    
    The simulation procedure then constructs the treatment effect
    heterogeneity through
    
    \deqn{
      \tau(X) =
        \pi_1(X) - \pi_0(X),
    }
    
    after applying the specified odds-ratio calibration. The resulting
    \code{tau_x} values are treated as the true CATE values, and their
    mean is treated as the true ATE.
    
    The resulting true quantities are therefore
    
    \deqn{
      TRUE\_CATE_i = \tau(X_i)
    }
    
    and
    
    \deqn{
      TRUE\_ATE =
        \frac{1}{n}\sum_{i=1}^{n}\tau(X_i).
    }
  }
  
  \subsection{Outcome simulation}{
    After obtaining the true potential outcome probabilities, a binary
    outcome is generated for each simulation replication according to
    
    \deqn{
      Y_i^{obs} \sim
      Bernoulli\left(\pi_i^{factual}\right),
    }
    
    where the factual probability is determined by the observed treatment
    status.
    
    For treated observations, the simulated outcome is generated using
    the treated potential outcome probability, whereas for control
    observations it is generated using the control potential outcome
    probability.
    
    This process is repeated \code{sim} times.
  }
  
  \subsection{Causal model comparison}{
    For each simulation replication, the selected causal models are fitted
    to the newly simulated outcome and the original covariates and
    treatment assignment.
    
    The following methods can be compared:
      
      \itemize{
        \item Causal MARS (\code{causalMARS})
        \item Causal boosting (\code{causal_boosting})
        \item Causal forest (\code{causal_forest})
        \item Causal BART (\code{causal_BART})
        \item Additive causal forest (\code{causal_Additive_Forest})
        \item Parametric standardization
      }
    
    For causal BART, both an unadjusted BART model and a propensity-score
    adjusted BART model are evaluated.
  }
  
  \subsection{Performance measures}{
    Two performance measures are calculated.
    
    The root mean squared error (RMSE) is used to evaluate the estimation
    of the CATE:
      
      \deqn{
        RMSE =
          \sqrt{
            \frac{1}{n}
            \sum_{i=1}^{n}
            \left(
              \widehat{\tau}(X_i)-\tau(X_i)
              \right)^2
          }.
      }
    
    A smaller RMSE indicates better estimation of the heterogeneous
    treatment effects.
    
    The absolute relative bias (ARB) is used to evaluate the estimation
    of the ATE:
      
      \deqn{
  ARB =
  \left|
  \frac{\widehat{ATE}-TRUE\_ATE}
       {TRUE\_ATE}
  \right| \times 100
}
    
    A smaller ARB indicates that the estimated ATE is closer to the true
    ATE.
  }
  
  \subsection{Graphical comparison}{
    After all simulation replications are completed, the RMSE and ARB
    results are converted to long format and summarized using boxplots.
    The RMSE plot compares the distribution of CATE estimation errors,
    whereas the ARB plot compares the distribution of relative ATE bias
    across the competing methods.
  }
  
}
\references{
The simulation study follows the framework described by:
Wendling, T., Jung, K., Callahan, A., Schuler, A., Shah, N. H., Gallego, B. (2018).
Comparing methods for estimation of heterogeneous treatment effects using observational data from health care databases.
\emph{Statistics in Medicine}.
\url{https://doi.org/10.1002/sim.7820}
}
\seealso{
  \code{\link{simulate_causal_data}},
  \code{\link{causalBART}},
  \code{\link{causalMARS}},
  \code{\link{causal_boosting}},
  \code{\link{causalForest}},
  \code{\link{causal_Additive_Forest}},
  \code{\link{parametric_standardization}}
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
    
    result <- compare_model(Y = Y, X = X, treatment = treatment, model_data_simulation = "logit", sim = 100,
      compare_model = c("causalMARS", "causal_boosting", "causal_forest", "causal_BART", "causal_Additive_Forest"))
    
    result$RMSE_df
    result$ARB_df
    
    result$RMSE_plot
    result$ARB_plot
  }
}
\keyword{simulation}
\keyword{causal}
\keyword{model comparison}
\keyword{performance}
