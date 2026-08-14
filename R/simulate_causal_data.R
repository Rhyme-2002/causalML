#' Simulate Data for Causal Machine-Learning Methods
#'
#' Generates potential outcome probabilities and a simulated binary
#' outcome for evaluating causal machine-learning methods. The
#' simulation framework follows the general approach of Wendling et al.
#' (2018), adapted to estimate heterogeneous treatment effects for a
#' binary treatment and binary outcome.
#'
#' @param Y Binary outcome variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param X A model matrix or data frame containing the observed covariates.
#' @param treatment Binary treatment variable coded as numeric 0 and 1. Must be numeric and not a factor.
#' @param model Character string specifying the model used to estimate
#'   the conditional outcome probabilities. Available methods are
#'   \code{"logit"}, \code{"BART"}, \code{"GBM"},
#'   \code{"Random_Forest"}, and \code{"KNN"}. The default is
#'   \code{"logit"}.
#' @param k Number of nearest neighbors used when
#'   \code{model = "KNN"}. The default is 5.
#' @param OR The prespecified odds ratio used to calibrate the
#'   treatment effect. If \code{NULL}, the odds ratio is estimated
#'   from a logistic regression model fitted to the observed outcome,
#'   treatment, and covariates.
#' @param seed Optional integer used to set the random-number seed.
#'
#' @return A list containing:
#' \describe{
#'   \item{simulated_data}{
#'     A data frame containing the estimated untreated probability
#'     (\code{pi0}), the original treated probability (\code{pi1}),
#'     the adjusted treated probability (\code{pi1_adjusted}),
#'     the factual outcome probability (\code{pi_factual}), and the
#'     true conditional treatment effect (\code{tau_x}).
#'   }
#'
#'   \item{X}{
#'     The covariate matrix used in the simulation.
#'   }
#'
#'   \item{treatment}{
#'     The binary treatment variable.
#'   }
#'
#'   \item{pi1}{
#'     Estimated conditional outcome probability under treatment
#'     before odds-ratio calibration.
#'   }
#'
#'   \item{pi0}{
#'     Estimated conditional outcome probability under control.
#'   }
#'
#'   \item{pi1_adjusted}{
#'     Adjusted conditional outcome probability under treatment after
#'     calibration using the specified odds ratio.
#'   }
#'
#'   \item{pi_factual}{
#'     The conditional probability of the observed factual outcome,
#'     determined by the observed treatment assignment.
#'   }
#'
#'   \item{TRUE_CATE}{
#'     A numeric vector containing the true conditional average
#'     treatment effects.
#'   }
#'
#'   \item{TRUE_ATE}{
#'     The true average treatment effect, calculated as the mean of
#'     the true CATE values.
#'   }
#'
#'   \item{model}{
#'     The fitted model used to estimate the conditional outcome
#'     probabilities.
#'   }
#'
#'   \item{method}{
#'     The name of the model used for simulation.
#'   }
#'
#'   \item{OR}{
#'     The odds ratio used to calibrate the treatment effect.
#'   }
#' }
#'
#' @details
#' This function implements a simulation framework for evaluating
#' heterogeneous treatment-effect estimators. The framework is
#' adapted from Wendling et al. (2018), who considered simulation
#' settings with binary treatment and binary outcomes and used the
#' conditional risk difference as the treatment-effect estimand.
#'
#' First, the conditional outcome probabilities under treatment and
#' control are estimated from the observed data using one of five
#' models: logistic regression, Bayesian Additive Regression Trees
#' (BART), gradient boosting machines (GBM), random forest, or
#' K-nearest neighbors (KNN).
#'
#' Let \eqn{X} denote the observed covariates and \eqn{Z} the binary
#' treatment indicator. The fitted model is used to obtain
#'
#' \deqn{
#' \pi^1(x) = P(Y = 1 \mid Z = 1, X = x)
#' }
#'
#' and
#'
#' \deqn{
#' \pi^0(x) = P(Y = 1 \mid Z = 0, X = x).
#' }
#'
#' The estimated potential outcome probabilities are then used to
#' construct the treatment-effect function on the logit scale.
#' Specifically,
#'
#' \deqn{
#' \gamma(x) =
#' \operatorname{logit}\{\pi^1(x)\}
#' -
#' \operatorname{logit}\{\pi^0(x)\}
#' + c,
#' }
#'
#' where \eqn{c} is a calibration constant. The constant is selected
#' so that the mean of the treatment-effect function on the logit
#' scale corresponds to the logarithm of the specified odds ratio:
#'
#' \deqn{
#' \frac{1}{n}\sum_{i=1}^{n}\gamma(X_i) = \log(OR).
#' }
#'
#' The adjusted treated potential outcome probability is then
#' obtained using the inverse-logit transformation:
#'
#' \deqn{
#' \pi^1_{\mathrm{adjusted}}(x) =
#' \operatorname{logit}^{-1}
#' \left[
#' \operatorname{logit}\{\pi^0(x)\} + \gamma(x)
#' \right].
#' }
#'
#' The conditional treatment effect is defined on the probability
#' scale as the conditional risk difference:
#'
#' \deqn{
#' \tau(x) =
#' \pi^1_{\mathrm{adjusted}}(x) - \pi^0(x).
#' }
#'
#' The true average treatment effect is calculated as
#'
#' \deqn{
#' ATE =
#' \frac{1}{n}\sum_{i=1}^{n}\tau(X_i).
#' }
#'
#' The factual outcome probability is determined by the observed
#' treatment assignment:
#'
#' \deqn{
#' \pi_{\mathrm{factual},i} =
#' Z_i\pi^1_{\mathrm{adjusted}}(X_i)
#' +
#' (1-Z_i)\pi^0(X_i).
#' }
#'
#' The binary outcome is then generated from a Bernoulli distribution:
#'
#' \deqn{
#' Y_i \sim
#' \operatorname{Bernoulli}
#' \left(\pi_{\mathrm{factual},i}\right).
#' }
#'
#' Therefore, the simulation provides known values of the true CATE
#' and ATE, which can be used to evaluate the performance of causal
#' machine-learning estimators using measures such as RMSE and
#' absolute relative bias.
#'
#' If \code{OR = NULL}, the odds ratio is estimated from a logistic
#' regression model using the observed outcome, treatment, and
#' covariates.
#'
#' Small numerical bounds are applied to the estimated probabilities
#' before calculating the logit transformation to avoid numerical
#' problems when probabilities are exactly 0 or 1.
#'
#' @references
#' Wendling, T., Jung, K., Callahan, A., Schuler, A., Shah, N. H.,
#' & Gallego, B. (2018).
#' Comparing methods for estimation of heterogeneous treatment effects
#' using observational data from health care databases.
#' \emph{Statistics in Medicine}, 37(23), 3309--3324.
#' \doi{10.1002/sim.7820}.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' X <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' treatment <- rbinom(100, 1, 0.5)
#' Y <- rbinom(100, 1, 0.5)
#'
#' # Logistic regression
#' result <- simulate_causal_data(Y = Y, X = X, treatment = treatment, model = "logit", OR = 1.5, seed = 123)
#'
#' result$simulated_data
#' result$TRUE_CATE
#' result$TRUE_ATE
#'
#' # BART
#' result_bart <- simulate_causal_data(Y = Y, X = X, treatment = treatment, model = "BART", OR = 1.5, seed = 123)
#' }
#'
#' @importFrom stats glm predict qlogis plogis
#' @export
simulate_causal_data <- function(Y, X, treatment, model = "logit", k = 5, OR = NULL, seed = NULL){
  
  allowed_sim_model <- c("logit", "BART", "GBM", "Random_Forest", "KNN")
  if(!(model %in% allowed_sim_model)){
    stop("Invalid model in model_data_simulation.")
  }
  X1 <- X
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  model <- match.arg(model, choices = allowed_sim_model)
  
  X1 <- as.matrix(X1)
  Y <- as.numeric(Y)
  treatment <- as.numeric(treatment)
  
  if(length(Y) != nrow(X1)){
    stop("Y and X1 have different numbers of observations.")
  }
  
  if(length(treatment) != nrow(X1)){
    stop("treatment and X1 have different numbers of observations.")
  }
  
  if(!all(treatment %in% c(0, 1))){
    stop("treatment must contain only 0 and 1.")
  }
  
  if(!all(Y %in% c(0, 1))){
    stop("Y must contain only 0 and 1.")
  }
  
  
  X_obs <- data.frame(treatment = treatment, X1)
  X_treat <- data.frame(treatment = 1, X1)
  X_control <- data.frame(treatment = 0, X1)
  
  
  ## Logistic Regression
  
  if(model == "logit"){
    dat <- data.frame(Y = Y, X_obs)
    fit <- glm(Y ~ ., data = dat, family = binomial())
    
    pi1 <- predict(fit, newdata = X_treat, type = "response")
    pi0 <- predict(fit, newdata = X_control,type = "response")
  }
  
  ## BART
  
  if(model == "BART"){
    if(!requireNamespace("BART", quietly = TRUE)){
      stop("Package 'BART' is required. ",
           "Install it using install.packages('BART').")
    }
    
    fit <- BART::pbart(x.train = as.matrix(X_obs), y.train = Y)
    pi1 <- predict(fit, newdata = as.matrix(X_treat))$prob.test.mean
    pi0 <- predict(fit, newdata = as.matrix(X_control))$prob.test.mean
  }
  
  ## Gradient Boosting Machine
  
  
  if(model == "GBM"){
    if (!requireNamespace("gbm", quietly = TRUE)){
      stop("Package 'gbm' is required. ",
           "Install it using install.packages('gbm').")
    }
    dat <- data.frame(Y = Y, X_obs)
    fit <- gbm::gbm(Y ~ ., data = dat, distribution = "bernoulli", n.trees = 500,
                    interaction.depth = 3, shrinkage = 0.01, verbose = FALSE)
    pi1 <- predict(fit, newdata = X_treat, n.trees = 500, type = "response")
    pi0 <- predict(fit, newdata = X_control, n.trees = 500, type = "response")
  }
  ## Random Forest
  if(model == "Random_Forest"){
    
    if(!requireNamespace("randomForest", quietly = TRUE)){
      stop("Package 'randomForest' is required. ",
           "Install it using install.packages('randomForest').")
    }
    
    dat <- data.frame(Y = factor(Y, levels = c(0, 1)), X_obs)
    fit <- randomForest::randomForest(Y ~ ., data = dat, ntree = 500, importance = TRUE)
    
    pi1 <- predict(fit, newdata = X_treat, type = "prob")[, "1"]
    pi0 <- predict(fit, newdata = X_control, type = "prob")[, "1"]
  }
  
  ## KNN
  if(model == "KNN"){
    if(!requireNamespace("class", quietly = TRUE)){
      stop("Package 'class' is required. ",
           "Install it using install.packages('class').")
    }
    
    dat <- data.frame(Y = factor(Y, levels = c(0, 1)), X_obs)
    fit <- list(method = "KNN", k = k, X = X_obs, Y = Y)
    
    knn_probability <- function(newdata){
      distances <- as.matrix(dist(rbind(as.matrix(X_obs), as.matrix(newdata))))
      n <- nrow(X_obs)
      prob <- numeric(nrow(newdata))
      
      for(i in seq_len(nrow(newdata))){
        d <- distances[n + i, seq_len(n)]
        nearest <- order(d)[seq_len(min(k, n))]
        prob[i] <- mean(Y[nearest])
      }
      
      return(prob)
    }
    
    pi1 <- knn_probability(X_treat)
    pi0 <- knn_probability(X_control)
  }
  
  if(is.null(OR)){
    dat <- data.frame(Y = Y, X_obs)
    fit <- glm(Y ~ ., data = dat, family = binomial())
    OR <- as.numeric(exp(coef(fit)[["treatment"]]))
  }
  
  eps <- 1e-6
  pi1 <- ifelse(pi1 <= 0, eps, ifelse(pi1 >= 1, 1 - eps, pi1))
  pi0 <- ifelse(pi0 <= 0, eps, ifelse(pi0 >= 1, 1 - eps, pi0)) 
  
  
  c_value <- log(OR) * 2 - mean(qlogis(pi1) - qlogis(pi0))
  gamma_x <- qlogis(pi1) - qlogis(pi0) + c_value
  pi1_adjusted <- plogis(qlogis(pi0) + gamma_x)
  tau_x <- pi1_adjusted - pi0
  pi_factual <- ifelse(treatment == 1, pi1_adjusted, pi0)
  TRUE_ATE <- mean(tau_x)
  simulated_data <- data.frame(pi1 = pi1, pi1_adjusted = pi1_adjusted, pi0 = pi0,
                               pi_factual = pi_factual, tau_x = tau_x)
  
  return(list(simulated_data = simulated_data, X = X1, treatment = treatment, pi1 = pi1, 
              pi0 = pi0, pi1_adjusted = pi1_adjusted, pi_factual = pi_factual, TRUE_CATE = tau_x,
              TRUE_ATE = TRUE_ATE, model = fit, method = model, OR = OR))
}
