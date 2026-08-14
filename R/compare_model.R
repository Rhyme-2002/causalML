#' Compare Causal Machine Learning Models
#'
#' Compares the performance of multiple causal inference and causal
#' machine learning methods for estimating heterogeneous treatment
#' effects through a simulation study.
#'
#' The function first uses \code{\link{simulate_causal_data}} to construct
#' potential outcome probabilities and the true conditional average
#' treatment effect (CATE). Binary outcomes are then repeatedly generated
#' from the simulated factual outcome probabilities. Each selected causal
#' model is fitted to the simulated data, and its estimated CATE and
#' average treatment effect (ATE) are compared with the corresponding
#' true values.
#'
#' @param Y Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param X A model matrix containing the covariates.
#' @param treatment Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param model_data_simulation Character string specifying the method
#'   used to estimate the conditional outcome probabilities in the
#'   data-generation procedure. Available methods are
#'   \code{"logit"}, \code{"BART"}, \code{"GBM"},
#'   \code{"Random_Forest"}, and \code{"KNN"}. The default is
#'   \code{"logit"}.
#' @param sim Number of simulation replications. The default is 100.
#' @param compare_model Character vector specifying the causal models
#'   to be compared. Available models are \code{"causalMARS"},
#'   \code{"causal_boosting"}, \code{"causal_forest"},
#'   \code{"causal_BART"}, and \code{"causal_Additive_Forest"}.
#'   The default compares all available causal machine learning models.
#' @param ... Additional arguments passed to the corresponding model
#'   fitting functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{RMSE_df}}{
#'     A data frame containing the CATE root mean squared error (RMSE)
#'     for each model across all simulation replications.
#'   }
#'   \item{\code{RMSE_long}}{
#'     A long-format data frame containing the CATE RMSE values with
#'     simulation number and model name.
#'   }
#'   \item{\code{ARB_df}}{
#'     A data frame containing the absolute relative bias (ARB) of the
#'     estimated ATE for each model across all simulation replications.
#'   }
#'   \item{\code{ARB_long}}{
#'     A long-format data frame containing the ARB values with simulation
#'     number and model name.
#'   }
#'   \item{\code{RMSE_plot}}{
#'     A boxplot showing the distribution of CATE RMSE across simulation
#'     replications for each model.
#'   }
#'   \item{\code{ARB_plot}}{
#'     A boxplot showing the distribution of ATE ARB across simulation
#'     replications for each model.
#'   }
#'   \item{\code{percentage_of_exposure_per_iteration}}{
#'     The percentage of observations receiving the simulated outcome
#'     value 1 in the final simulation iteration.
#'   }
#' }
#'
#' @details
#' The simulation-based model comparison consists of the following
#' steps.
#'
#' \strong{Step 1: Data-generation procedure}
#'
#' The function first calls \code{\link{simulate_causal_data}} using
#' the observed outcome, covariates, and treatment assignment. The
#' selected \code{model_data_simulation} method is used to estimate
#' the conditional outcome probabilities under treatment and control.
#' These probabilities are used to construct the simulated potential
#' outcomes and the underlying treatment-effect surface.
#'
#' The simulation procedure is motivated by the framework proposed by
#' Wendling et al. (2018), in which observed covariate and treatment
#' assignment structures are retained while outcomes are generated from
#' estimated conditional outcome models.
#'
#' The resulting true conditional treatment effect is
#'
#' \deqn{
#' \tau(X_i) =
#' \pi_1(X_i) - \pi_0(X_i),
#' }
#'
#' where \eqn{\pi_1(X_i)} and \eqn{\pi_0(X_i)} represent the simulated
#' potential outcome probabilities under treatment and control,
#' respectively.
#'
#' The true average treatment effect is calculated as
#'
#' \deqn{
#' ATE =
#' \frac{1}{n}\sum_{i=1}^{n}\tau(X_i).
#' }
#'
#' \strong{Step 2: Generate simulated outcomes}
#'
#' For each simulation replication, a new binary outcome is generated
#' from the factual outcome probability returned by
#' \code{\link{simulate_causal_data}}.
#'
#' Specifically, for observation \eqn{i}, the simulated outcome is
#' generated according to
#'
#' \deqn{
#' Y_i^{*} \sim \operatorname{Bernoulli}
#' \left(\pi_{\mathrm{factual},i}\right).
#' }
#'
#' The factual probability is determined by the observed treatment
#' assignment, such that treated observations use the treated potential
#' outcome probability and control observations use the control
#' potential outcome probability.
#'
#' \strong{Step 3: Fit causal models}
#'
#' Each selected method in \code{compare_model} is fitted separately
#' to the simulated outcome, while the covariates and treatment
#' assignment are kept fixed.
#'
#' The available causal models are:
#'
#' \itemize{
#'   \item \code{causalMARS}: Causal Multivariate Adaptive Regression
#'   Splines.
#'   \item \code{causal_boosting}: Causal boosting.
#'   \item \code{causal_forest}: Causal forest.
#'   \item \code{causal_BART}: Bayesian Additive Regression Trees.
#'   \item \code{causal_Additive_Forest}: Additive causal forest.
#'   \item \code{parametric_standardization}: Logistic-regression-based
#'   parametric standardization, which is included as a benchmark
#'   method.
#' }
#'
#' When \code{causal_BART} is selected, both BART without propensity
#' score adjustment and BART with propensity score adjustment are
#' evaluated.
#'
#' Parametric standardization is automatically included as a benchmark
#' regardless of whether it is specified in \code{compare_model}.
#'
#' \strong{Step 4: Estimate CATE and ATE}
#'
#' For each fitted model, the estimated conditional treatment effects
#' are obtained and compared with the known simulated CATE values.
#' The estimated ATE is calculated from the estimated CATE values.
#'
#' \strong{Step 5: Evaluate CATE performance}
#'
#' The accuracy of CATE estimation is evaluated using the root mean
#' squared error (RMSE):
#'
#' \deqn{
#' RMSE =
#' \sqrt{
#' \frac{1}{n}
#' \sum_{i=1}^{n}
#' \left\{
#' \widehat{\tau}(X_i)-\tau(X_i)
#' \right\}^{2}
#' }.
#' }
#'
#' Smaller RMSE values indicate better estimation of the heterogeneous
#' treatment-effect function.
#'
#' \strong{Step 6: Evaluate ATE performance}
#'
#' The accuracy of ATE estimation is evaluated using absolute relative
#' bias (ARB):
#'
#' \deqn{
#' ARB =
#' \left|
#' \frac{\widehat{ATE}-ATE}{ATE}
#' \right|.
#' }
#'
#' Smaller ARB values indicate that the estimated ATE is closer to the
#' true simulated ATE.
#'
#' \strong{Step 7: Repeat the simulation}
#'
#' Steps 2--6 are repeated \code{sim} times. The RMSE and ARB from
#' every simulation replication are stored in separate data frames.
#' These results are subsequently converted to long format for
#' visualization.
#'
#' Boxplots are produced to compare the distribution of CATE RMSE and
#' ATE ARB across the competing methods.
#'
#' @references
#' Wendling, T., Jung, K., Callahan, A., Schuler, A., Shah, N. H.,
#' and Gallego, B. (2018).
#' Comparing methods for estimation of heterogeneous treatment effects
#' using observational data from health care databases.
#' \emph{Statistics in Medicine}, 37(23), 3309--3324.
#' \doi{10.1002/sim.7820}.
#'
#' @seealso
#' \code{\link{simulate_causal_data}},
#' \code{\link{causalMARS}},
#' \code{\link{causal_boosting}},
#' \code{\link{causalForest}},
#' \code{\link{causalBART}},
#' \code{\link{causal_Additive_Forest}},
#' \code{\link{parametric_standardization}}
#'
#' @examples
#' \dontrun{
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' x4 <- rnorm(100)
#' x5 <- rnorm(100)
#'
#' X <- model.matrix(~ x1 + x2 + x3 + x4 + x5)
#' Y <- rbinom(100, 1, 0.5)
#' treatment <- rbinom(100, 1, 0.5)
#'
#' # Compare causal models using a logistic outcome model
#' result <- compare_model(Y = Y, X = X, treatment = treatment, model_data_simulation = "logit", sim = 10,
#' compare_model = c("causalMARS", "causal_boosting", "causal_forest", "causal_BART", "causal_Additive_Forest"))
#'
#' # RMSE results
#' result$RMSE_df
#'
#' # ARB results
#' result$ARB_df
#'
#' # RMSE boxplot
#' result$RMSE_plot
#'
#' # ARB boxplot
#' result$ARB_plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_classic
#' @importFrom tidyr pivot_longer
#' @export
compare_model <- function(Y, X, treatment,
                          model_data_simulation = "logit",
                          sim = 100, compare_model = c("causalMARS", "causal_boosting",
                                                       "causal_forest", "causal_BART", "causal_Additive_Forest"), ...){
  
  allowed_sim_model <- c("logit", "BART", "GBM", "Random_Forest", "KNN")
  if(!(model_data_simulation %in% allowed_sim_model)) {
    stop("Invalid model in model_data_simulation.")
  }
  compare_model <- unique(compare_model)
  allowed_models <- c("causalMARS", "causal_boosting", "causal_forest", "causal_BART", "causal_Additive_Forest")
  if(!all(compare_model %in% allowed_models)) {
    stop("Invalid model in compare_model.")
  }
  sim_data <- simulate_causal_data(Y = Y, X = X, treatment = treatment,
                                   model = model_data_simulation)$simulated_data
  true_CATE <- sim_data$tau_x
  true_ATE <- mean(true_CATE)
  model_names <- setdiff(compare_model, "causal_BART")
  model_names <- c(model_names, "parametric_standardization")
  
  RMSE_df <- data.frame(Simulation = 1:sim, matrix(nrow = sim,
                                                   ncol = length(model_names),
                                                   dimnames = list(NULL, model_names)))
  ARB_df <- data.frame(Simulation = 1:sim, matrix(nrow = sim,
                                                  ncol = length(model_names),
                                                  dimnames = list(NULL, model_names)))
  percentage_of_exposure_per_iteration <- c()
  for(i in 1:sim){
    Y1 <- rbinom(n = nrow(X), size = 1, prob = sim_data$pi_factual)
    percentage_of_exposure_per_iteration[i] <- mean(Y1) * 100
    
    if("causal_BART" %in% compare_model){
      model_bart <- causalBART(X = X, Y = Y1, treatment = treatment, PS_adjusted = FALSE)
      CATE <- model_bart$CATE
      ATE <- model_bart$ATE
      RMSE_df[i, "BART without PS"] <- sqrt(mean((CATE - true_CATE)^2))
      ARB_df[i, "BART without PS"] <- abs((ATE - true_ATE) / true_ATE)
      
      model_bart_ps <- causalBART(X = X, Y = Y1, treatment = treatment, PS_adjusted = TRUE)
      CATE <- model_bart_ps$CATE
      ATE <- model_bart_ps$ATE
      RMSE_df[i, "BART with PS"] <- sqrt(mean((CATE - true_CATE)^2))
      ARB_df[i, "BART with PS"] <- abs((ATE - true_ATE) / true_ATE)
      
      print("BART END")
    }
    
    if("causalMARS" %in% compare_model){
      model_mars <- causalMARS(x = X, y = Y1, treatment = treatment)
      CATE <- model_mars$cate
      ATE <- model_mars$ate
      RMSE_df[i, "causalMARS"] <- sqrt(mean((CATE - true_CATE)^2))
      ARB_df[i, "causalMARS"] <- abs((ATE - true_ATE) / true_ATE)
      
      print("MARS END")
    }
    
    if("causal_boosting" %in% compare_model){
      model_boosting <- causal_boosting(X = X, y = Y1, treatment = treatment)
      CATE <- model_boosting$cate
      ATE <- model_boosting$ate
      RMSE_df[i, "causal_boosting"] <- sqrt(mean((CATE - true_CATE)^2))
      ARB_df[i, "causal_boosting"] <- round((abs((ATE - true_ATE) / true_ATE)) * 100, 2)
      print("CAUSAL BOOSTING END")
    }
    
    if("causal_forest" %in% compare_model){
      model_forest <- causalForest(X = X, Y = Y1, treatment = treatment)
      result_forest <- predict_causal(object = model_forest, newdata = X)
      CATE <- result_forest$CATE
      ATE <- result_forest$ATE
      RMSE_df[i, "causal_forest"] <- sqrt(mean((CATE$CATE - true_CATE)^2))
      ARB_df[i, "causal_forest"] <- abs((ATE - true_ATE) / true_ATE)
      print("CAUSAL FOREST END")
    }
    
    if("causal_Additive_Forest" %in% compare_model){
      model_addative_forest <- causal_Additive_Forest(X = X, Y = Y1, treatment = treatment)
      CATE <- model_addative_forest$CATE
      ATE <- model_addative_forest$ATE
      RMSE_df[i, "causal_Additive_Forest"] <- sqrt(mean((CATE - true_CATE)^2))
      ARB_df[i, "causal_Additive_Forest"] <- abs((ATE - true_ATE) / true_ATE)
      print("ADDITIVE CAUSAL FOREST END")
    }
    
    model_parametric_standardization <- parametric_standardization(X = X, Y = Y1, treatment = treatment)
    CATE <- model_parametric_standardization$CATE
    ATE <- model_parametric_standardization$ATE
    RMSE_df[i, "parametric_standardization"] <- sqrt(mean((CATE - true_CATE)^2))
    ARB_df[i, "parametric_standardization"] <- abs((ATE - true_ATE) / true_ATE)
  }
  
  RMSE_long <- RMSE_df %>%
    tidyr::pivot_longer(cols = -Simulation, names_to = "Model", values_to = "RMSE")
  
  RMSE_plot <- ggplot2::ggplot(RMSE_long, ggplot2::aes(x = Model, y = RMSE)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x = "Model", y = "RMSE", title = "RMSE Across Simulations") +
    ggplot2::theme_classic()
  
  ARB_long <- ARB_df %>%
    tidyr::pivot_longer(cols = -Simulation, names_to = "Model", values_to = "ARB")
  
  ARB_plot <- ggplot2::ggplot(ARB_long, ggplot2::aes(x = Model, y = ARB)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x = "Model", y = "ARB", title = "ARB Across Simulations") +
    ggplot2::theme_classic()
  
  result <- list(RMSE_df = RMSE_df, RMSE_long = RMSE_long,
                 ARB_df = ARB_df, ARB_long = ARB_long,
                 RMSE_plot = RMSE_plot, ARB_plot = ARB_plot,
                 percentage_of_exposure_per_iteration = percentage_of_exposure_per_iteration)
  
  return(result)
}
