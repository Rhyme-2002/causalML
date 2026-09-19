construct_3way_table <- function(Y, Z, P_C, RD_CZ, RD_CY) {
  
  keep <- !is.na(Y) & !is.na(Z)
  Y <- Y[keep]
  Z <- Z[keep]
  
  if (!all(Y %in% c(0, 1))) {
    stop("Y must contain only 0 and 1.")
  }
  
  if (!all(Z %in% c(0, 1))) {
    stop("Z must contain only 0 and 1.")
  }
  
  if (P_C < 0 || P_C > 1) {
    stop("P_C must be between 0 and 1.")
  }
  
  N <- length(Y)
  P_Z <- mean(Z == 1)
  P_Y_Z0 <- mean(Y[Z == 0] == 1)
  P_Y_Z1 <- mean(Y[Z == 1] == 1)
  
  b11 <- N * (P_Z * P_C + RD_CZ * P_C * (1 - P_C))
  b10 <- N * P_C - b11
  b01 <- N * P_Z - b11
  b00 <- N - b11 - b10 - b01
  
  if ((b10 + b00) == 0 || (b11 + b01) == 0) {
    return(NULL)
  }
  
  P_C_Z0 <- b10 / (b10 + b00)
  P_C_Z1 <- b11 / (b11 + b01)
  n011 <- b10 * (P_Y_Z0 + RD_CY * (1 - P_C_Z0))
  n001 <- b10 - n011
  n010 <- b00 * (P_Y_Z0 - RD_CY * P_C_Z0)
  n000 <- b00 - n010
  n111 <- b11 * (P_Y_Z1 + RD_CY * (1 - P_C_Z1))
  n101 <- b11 - n111
  n110 <- b01 * (P_Y_Z1 - RD_CY * P_C_Z1)
  n100 <- b01 - n110
  
  frequency <- data.frame(
    Z = c(0, 0, 0, 0, 1, 1, 1, 1),
    Y = c(0, 0, 1, 1, 0, 0, 1, 1),
    C = c(0, 1, 0, 1, 0, 1, 0, 1),
    Frequency = c(n000, n001, n010, n011, n100, n101, n110, n111))
  
  if (any(!is.finite(frequency$Frequency))) {
    return(NULL)
  }
  
  if (any(frequency$Frequency < 0)) {
    return(NULL)
  }
  
  x <- frequency$Frequency
  target <- round(sum(x))
  result <- floor(x)
  n <- target - sum(result)
  
  if (n > 0) {
    decimal <- x - floor(x)
    index <- order(decimal, decreasing = TRUE)[seq_len(n)]
    result[index] <- result[index] + 1
  }
  
  frequency$Frequency <- result
  return(frequency)
}

risk_difference <- function(Y, Z) {
  
  outcome <- Y
  treatment  <- Z
  risk_treated <- mean(outcome[treatment == 1], na.rm = TRUE)
  risk_control <- mean(outcome[treatment == 0], na.rm = TRUE)
  rd <- risk_treated - risk_control
  
  return(rd)
}

distribute_C <- function(Y, Z, frequency_table, set_seed = 102) {
  
  set.seed(set_seed)
  C <- rep(NA_integer_, length(Y))
  for (z_value in 0:1) {
    for (y_value in 0:1) {
      n_c1 <- frequency_table$Frequency[frequency_table$Z == z_value & frequency_table$Y == y_value & frequency_table$C == 1]
      n_c1 <- as.integer(n_c1)
      
      index <- which(!is.na(Y) & !is.na(Z) & Y == y_value & Z == z_value)
      if (n_c1 > length(index)) {
        stop(paste0("Required C = 1 frequency (", n_c1,") exceeds available observations (",length(index), ") for Z = ", z_value, ", Y = ", y_value))
      }
      
      C[index] <- 0
      if (n_c1 > 0) {
        selected <- sample(index, n_c1, replace = FALSE)
        C[selected] <- 1
      }
    }
  }
  
  return(C)
}

#' Sensitivity Analysis for an Unmeasured Binary Confounder
#'
#' Performs a quantitative sensitivity analysis to assess how an unmeasured
#' binary confounder may affect the estimated average treatment effect (ATE).
#' The analysis varies the prevalence of the unmeasured confounder and its
#' associations with the treatment and outcome, represented by risk
#' differences.
#'
#' For each combination of confounder prevalence, treatment-confounder risk
#' difference, and outcome-confounder risk difference, a three-way frequency
#' table is constructed and a binary confounder is assigned to the observed
#' data. The selected causal inference method is then refitted after adding
#' the simulated confounder to the observed covariates.
#'
#' When \code{Simulation = 1}, the function returns the adjusted ATE and
#' percentage changes relative to the observed ATE. When
#' \code{Simulation > 1}, the function summarizes the adjusted ATE across
#' simulations using the mean, standard deviation, and empirical confidence
#' interval.
#'
#' The function also generates text plots and heatmaps for visualizing the
#' adjusted ATE and its percentage change across the sensitivity-analysis
#' parameter space.
#'
#' @param Y A numeric binary outcome vector containing only 0 and 1.
#' @param treatment A numeric binary treatment vector containing only 0 and 1.
#' @param X A matrix or data frame containing the observed covariates used
#'   for adjustment in the causal effect estimation.
#' @param P_C A numeric vector specifying the assumed prevalence of the
#'   unmeasured binary confounder. Values must lie between 0 and 1.
#'   Default is \code{seq(0.1, 0.9, 0.1)}.
#' @param RD_CZ A numeric vector specifying the assumed risk difference
#'   between the unmeasured confounder (\code{C}) and treatment (\code{Z}).
#'   Values are typically specified between -1 and 1.
#'   Default is \code{seq(-1, 1, 0.1)}.
#' @param RD_CY A numeric vector specifying the assumed risk difference
#'   between the unmeasured confounder (\code{C}) and outcome (\code{Y}).
#'   Values are typically specified between -1 and 1.
#'   Default is \code{seq(-1, 1, 0.1)}.
#' @param causal_method Character string specifying the causal inference
#'   method used to estimate the ATE. Supported methods are
#'   \code{"causalForest"}, \code{"causal_Additive_Forest"},
#'   \code{"causal_boosting"}, \code{"causalMARS"}, \code{"causalBART"},
#'   and \code{"parametric_standardization"}.
#' @param Simulation Integer specifying the number of simulations performed
#'   for each sensitivity-analysis parameter combination. If
#'   \code{Simulation = 1}, the adjusted ATE from the single simulation is
#'   returned. If \code{Simulation > 1}, the mean, standard deviation, and
#'   empirical confidence interval of the adjusted ATE are calculated.
#'   Default is \code{1}.
#' @param ci_level Numeric value between 0 and 1 specifying the confidence
#'   level for the empirical confidence interval when \code{Simulation > 1}.
#'   Default is \code{0.95}.
#' @param digits Integer specifying the number of decimal places used for
#'   reported sensitivity-analysis results. Default is \code{3}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{Observed_ATE}{
#'     The ATE estimated using the original observed covariates without the
#'     simulated unmeasured confounder.
#'   }
#'   \item{N_Simulations}{
#'     The number of simulations performed for each parameter combination.
#'   }
#'   \item{CI_Level}{
#'     The confidence level used for the empirical confidence intervals.
#'     Returned when \code{Simulation > 1}.
#'   }
#'   \item{Results}{
#'     A data frame containing the sensitivity-analysis results. For a single
#'     simulation, this includes the assumed confounder prevalence,
#'     treatment-confounder risk difference, outcome-confounder risk
#'     difference, adjusted ATE, absolute percentage change, and directional
#'     percentage change. For multiple simulations, it additionally includes
#'     the standard deviation, lower and upper confidence limits, and number
#'     of valid simulations.
#'   }
#'   \item{Raw_Effect_Matrix}{
#'     A matrix containing the adjusted ATE obtained from each simulation for
#'     every sensitivity-analysis parameter combination. Returned when
#'     \code{Simulation > 1}.
#'   }
#'   \item{graph_adjusted_ATE}{
#'     A list of ggplot2 text plots displaying the adjusted ATE for each
#'     assumed confounder prevalence. For multiple simulations, the labels
#'     display the mean adjusted ATE and its confidence interval.
#'   }
#'   \item{graph_percentage_change_ATE}{
#'     A list of ggplot2 text plots displaying the absolute percentage change
#'     in the estimated ATE for each assumed confounder prevalence.
#'   }
#'   \item{graph_percentage_change_with_direction_ATE}{
#'     A list of ggplot2 text plots displaying the directional percentage
#'     change in the estimated ATE for each assumed confounder prevalence.
#'   }
#'   \item{graph_adjusted_ATE_heat_map}{
#'     A list of ggplot2 heatmaps displaying the adjusted ATE across the
#'     treatment-confounder and outcome-confounder risk differences.
#'   }
#'   \item{graph_percentage_change_ATE_heat_map}{
#'     A list of ggplot2 heatmaps displaying the absolute percentage change
#'     in the ATE.
#'   }
#'   \item{graph_percentage_change_with_direction_ATE_heat_map}{
#'     A list of ggplot2 heatmaps displaying the directional percentage
#'     change in the ATE. A diverging colour scale is used to distinguish
#'     increases and decreases in the estimated effect.
#'   }
#' }
#'
#' @details
#' The observed ATE is first estimated using the selected causal inference
#' method and the observed covariates. For each sensitivity-analysis
#' parameter combination, \code{construct_3way_table()} is used to construct
#' a feasible joint distribution of the binary outcome, treatment, and
#' unmeasured confounder.
#'
#' Parameter combinations producing invalid or negative cell frequencies are
#' treated as infeasible and are excluded from the analysis.
#'
#' The simulated confounder is assigned to the observed observations using
#' \code{distribute_C()}, after which the causal inference method is refitted
#' using the original covariates together with the simulated confounder.
#'
#' The absolute percentage change is calculated as
#'
#' \deqn{
#' \frac{|\mathrm{ATE}*{adjusted} - \mathrm{ATE}*{observed}|}
#' {|\mathrm{ATE}*{observed}|} \times 100.
#' }
#'
#' The directional percentage change is calculated as
#'
#' \deqn{
#' \frac{\mathrm{ATE}*{adjusted} - \mathrm{ATE}*{observed}}
#' {|\mathrm{ATE}*{observed}|} \times 100.
#' }
#'
#' For multiple simulations, the confidence interval is obtained from the
#' empirical quantiles of the simulated adjusted ATEs using the confidence
#' level specified by \code{ci_level}.
#'
#' @seealso
#' \code{\link{construct_3way_table}},
#' \code{\link{distribute_C}},
#' \code{\link{risk_difference}}
#'
#' @examples
#' \dontrun{
#' # Single sensitivity-analysis scenario
#' result <- sensitivity_analysis(
#'   Y = Y,
#'   treatment = treatment,
#'   X = X,
#'   P_C = 0.5,
#'   RD_CZ = 0.4,
#'   RD_CY = 0.5,
#'   causal_method = "causalForest",
#'   Simulation = 1
#' )
#'
#' result$Observed_ATE
#' result$Results
#'
#' # Sensitivity analysis over multiple parameter combinations
#' result <- sensitivity_analysis(
#'   Y = Y,
#'   treatment = treatment,
#'   X = X,
#'   P_C = seq(0.1, 0.9, 0.1),
#'   RD_CZ = seq(-1, 1, 0.1),
#'   RD_CY = seq(-1, 1, 0.1),
#'   causal_method = "causalForest",
#'   Simulation = 3,
#'   ci_level = 0.95
#' )
#'
#' # View the adjusted ATE results
#' head(result$Results)
#'
#' # Display the heatmap for the first prevalence level
#' result$graph_adjusted_ATE_heat_map[[1]]
#' }
#'
#' @import ggplot2
#' @importFrom stats quantile
#' @export  
sensitivity_analysis <- function(Y, treatment, X,
                                 P_C = seq(0.1, 0.9, 0.1),
                                 RD_CZ = seq(-1, 1, 0.1),
                                 RD_CY  = seq(-1, 1, 0.1),
                                 causal_method,
                                 Simulation = 1,
                                 ci_level = 0.95,
                                 digits = 3){
  
  run_method <- function(Yv, Xv, treatv) {
    if (causal_method == "causalForest") {
      causalForest(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else if (causal_method == "causal_Additive_Forest") {
      causal_Additive_Forest(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else if (causal_method == "causal_boosting") {
      causal_boosting(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else if (causal_method == "causalMARS") {
      causalMARS(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else if (causal_method == "causalBART") {
      causalBART(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else if (causal_method == "parametric_standardization") {
      parametric_standardization(Y = Yv, X = Xv, treatment = treatv)$ATE
    } else {
      stop("Invalid causal method.")
    }
  }
  
  obs_ATE <- run_method(Y, X, treatment)
  
  grid   <- expand.grid(P_C = P_C, RD_CZ = RD_CZ, RD_CY = RD_CY, KEEP.OUT.ATTRS = FALSE)
  n_grid <- nrow(grid)
  
  single_point <- (length(P_C) == 1 && length(RD_CZ) == 1 && length(RD_CY) == 1)
  
  # ---------------------------------------------------------------------
  # Feasibility check (done ONCE - it does not depend on the simulation)
  # ---------------------------------------------------------------------
  tables   <- vector("list", n_grid)
  feasible <- logical(n_grid)
  for (row in seq_len(n_grid)) {
    tab <- construct_3way_table(Y = Y, Z = treatment,
                                P_C = grid$P_C[row],
                                RD_CZ = grid$RD_CZ[row],
                                RD_CY = grid$RD_CY[row])
    if (!is.null(tab)) {
      tables[[row]] <- tab
      feasible[row] <- TRUE
    }
  }
  
  if (single_point && !feasible[1]) {
    cat("\n====================================\n")
    cat("This combination is NOT feasible:",
        "P_C =", grid$P_C[1],
        "| RD_CZ =", grid$RD_CZ[1],
        "| RD_CY =", grid$RD_CY[1], "\n")
    cat("====================================\n")
    names(grid) <- c("Confounder_Prevalence", "RD_Confounder_Treatment", "RD_Confounder_Outcome")
    return(invisible(list(
      Observed_ATE  = obs_ATE,
      N_Simulations = Simulation,
      Results       = grid[0, ]
    )))
  }
  
  if (!single_point && any(!feasible)) {
    cat("\n", sum(!feasible), " of ", n_grid,
        " combinations are NOT feasible and will be skipped.\n", sep = "")
  }
  
  effect_mat <- matrix(NA_real_, nrow = n_grid, ncol = Simulation)
  
  # ---------------------------------------------------------------------
  # Timing helper (seconds)
  # ---------------------------------------------------------------------
  fmt_sec <- function(secs) paste0(round(secs), " seconds")
  start_time <- Sys.time()
  iter_times <- numeric(Simulation)
  
  for (i in seq_len(Simulation)) {
    
    iter_start <- Sys.time()
    
    cat("\n====================================\n")
    cat("Simulation", i, "of", Simulation, "\n")
    cat("====================================\n")
    
    for (row in seq_len(n_grid)) {
      
      # Progress percentage
      if (row %% 100 == 0 || row == 1 || row == n_grid) {
        percentage <- round(row / n_grid * 100, 1)
        cat(
          "\rSimulation ", i, "/", Simulation,
          " | Combination ", row, "/", n_grid,
          " | Progress: ", percentage, "%",
          sep = ""
        )
        flush.console()
      }
      
      if (!feasible[row]) next
      
      conf_column <- distribute_C(Y = Y, Z = treatment,
                                  frequency_table = tables[[row]],
                                  set_seed = i)
      new_X <- cbind(X, C = conf_column)
      effect_mat[row, i] <- run_method(Y, new_X, treatment)
    }
    
    cat("\nSimulation", i, "completed.\n")
    
    # Time report (only when Simulation > 1)
    iter_times[i] <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
    if (Simulation > 1) {
      elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time  <- mean(iter_times[1:i])
      remaining <- avg_time * (Simulation - i)
      total_est <- elapsed + remaining
      
      cat("Current simulation time   :", fmt_sec(iter_times[i]), "\n")
      cat("Average simulation time   :", fmt_sec(avg_time), "\n")
      cat("Elapsed time              :", fmt_sec(elapsed), "\n")
      cat("Estimated remaining time  :", fmt_sec(remaining), "\n")
      cat("Estimated total time      :", fmt_sec(total_est), "\n")
      flush.console()
    }
  }
  
  cat("\n====================================\n")
  cat("ALL SIMULATIONS COMPLETED\n")
  if (Simulation > 1) {
    cat("Total time:", fmt_sec(as.numeric(difftime(Sys.time(), start_time, units = "secs"))), "\n")
  }
  cat("====================================\n")
  
  names(grid) <- c("Confounder_Prevalence", "RD_Confounder_Treatment", "RD_Confounder_Outcome")
  
  make_text_plot <- function(df, label_col, title_prefix) {
    plots <- list()
    for (pc_val in P_C) {
      sub <- df[df$Confounder_Prevalence == pc_val, ]
      if (nrow(sub) == 0) next
      plots[[as.character(pc_val)]] <- sub |>
        ggplot2::ggplot(ggplot2::aes(x = RD_Confounder_Outcome, y = RD_Confounder_Treatment)) +
        ggplot2::geom_text(ggplot2::aes(label = .data[[label_col]]), vjust = -0.8, size = 3) +
        ggplot2::labs(x = "Risk Difference of response given confounder",
                      y = "Risk Difference of treatment given confounder",
                      title = paste0(title_prefix, pc_val)) +
        ggplot2::theme_minimal()
    }
    plots
  }
  
  make_heatmap <- function(df, fill_col, label_col, title_prefix, diverging = FALSE) {
    plots <- list()
    for (pc_val in P_C) {
      sub <- df[df$Confounder_Prevalence == pc_val, ]
      if (nrow(sub) == 0) next
      p <- sub |>
        ggplot2::ggplot(ggplot2::aes(x = factor(RD_Confounder_Outcome),
                                     y = factor(RD_Confounder_Treatment),
                                     fill = .data[[fill_col]])) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = .data[[label_col]]), size = 3) +
        ggplot2::labs(x = "Risk Difference of response given confounder",
                      y = "Risk Difference of treatment given confounder",
                      title = paste0(title_prefix, pc_val)) +
        ggplot2::theme_minimal()
      
      p <- if (diverging) {
        p + ggplot2::scale_fill_gradient2(name = fill_col, low = "blue", mid = "white", high = "red", midpoint = 0)
      } else {
        p + ggplot2::scale_fill_viridis_c(name = fill_col)
      }
      plots[[as.character(pc_val)]] <- p
    }
    plots
  }
  
  exclude_from_plot <- function(df) {
    rd_sum <- df$RD_Confounder_Outcome + df$RD_Confounder_Treatment
    (df$Confounder_Prevalence == rd_sum) & (rd_sum == 1)
  }
  
  if (Simulation == 1) {
    adjust_ATE     <- effect_mat[, 1]
    per_change     <- abs(adjust_ATE - obs_ATE) / abs(obs_ATE) * 100
    dir_per_change <- (adjust_ATE - obs_ATE) / abs(obs_ATE) * 100
    
    grid$Adjusted_ATE       <- round(adjust_ATE, digits)
    grid$Abs_Percent_Change <- round(per_change, digits)
    grid$Percent_Change     <- round(dir_per_change, digits)
    grid <- na.omit(grid)
    
    if (single_point) {
      return(list(
        Observed_ATE  = obs_ATE,
        N_Simulations = 1L,
        Results       = grid
      ))
    }
    
    grid_labeled <- grid
    grid_labeled$Label_Adjusted_ATE  <- as.character(round(grid_labeled$Adjusted_ATE, 4))
    grid_labeled$Label_Abs_Change    <- as.character(round(grid_labeled$Abs_Percent_Change, 4))
    grid_labeled$Label_Signed_Change <- as.character(round(grid_labeled$Percent_Change, 4))
    
    plot_data <- grid_labeled[!exclude_from_plot(grid_labeled), ]
    
    return(list(
      Observed_ATE  = obs_ATE,
      N_Simulations = 1L,
      Results       = grid,
      
      graph_adjusted_ATE = make_text_plot(plot_data, "Label_Adjusted_ATE",
                                          "Adjusted ATE | Prevalence of confounder: "),
      graph_percentage_change_ATE = make_text_plot(plot_data, "Label_Abs_Change",
                                                   "Percentage Change of ATE | Prevalence of confounder: "),
      graph_percentage_change_with_direction_ATE = make_text_plot(plot_data, "Label_Signed_Change",
                                                                  "Percentage Change of ATE | Prevalence of confounder: "),
      
      graph_adjusted_ATE_heat_map = make_heatmap(plot_data, "Adjusted_ATE", "Label_Adjusted_ATE",
                                                 "Adjusted ATE | Prevalence of confounder: "),
      graph_percentage_change_ATE_heat_map = make_heatmap(plot_data, "Abs_Percent_Change", "Label_Abs_Change",
                                                          "Percentage Change of ATE | Prevalence of confounder: "),
      graph_percentage_change_with_direction_ATE_heat_map = make_heatmap(plot_data, "Percent_Change",
                                                                         "Label_Signed_Change",
                                                                         "Percentage Change of ATE | Prevalence of confounder: ", diverging = TRUE)
    ))
  }
  
  # -----------------------------------------------------------------------
  # MULTIPLE REPS (Simulation > 1): mean + SD + CI across reps, per grid point
  # -----------------------------------------------------------------------
  alpha <- 1 - ci_level
  
  mean_ATE <- apply(effect_mat, 1, mean, na.rm = TRUE)
  sd_ATE   <- apply(effect_mat, 1, sd,   na.rm = TRUE)
  n_valid  <- apply(effect_mat, 1, function(r) sum(!is.na(r)))
  
  ci_lo <- apply(effect_mat, 1, function(r)
    if (all(is.na(r))) NA_real_ else quantile(r, alpha / 2, na.rm = TRUE))
  ci_hi <- apply(effect_mat, 1, function(r)
    if (all(is.na(r))) NA_real_ else quantile(r, 1 - alpha / 2, na.rm = TRUE))
  
  change         <- mean_ATE - obs_ATE
  per_change     <- abs(change) / abs(obs_ATE) * 100
  dir_per_change <- change / abs(obs_ATE) * 100
  
  grid$Adjusted_ATE                  <- round(mean_ATE, digits)
  grid$Adjusted_ATE_SD               <- round(sd_ATE, digits)
  grid$Adjusted_ATE_CI_Lower         <- round(ci_lo, digits)
  grid$Adjusted_ATE_CI_Upper         <- round(ci_hi, digits)
  grid$Abs_Percent_Change            <- round(per_change, digits)
  grid$Percent_Change_with_direction <- round(dir_per_change, digits)
  grid$N_Valid_Simulations           <- n_valid
  
  row_ok     <- !is.na(mean_ATE)
  grid       <- grid[row_ok, ]
  effect_mat <- effect_mat[row_ok, , drop = FALSE]
  
  if (single_point) {
    return(list(
      Observed_ATE      = obs_ATE,
      N_Simulations     = Simulation,
      CI_Level          = ci_level,
      Results           = grid,
      Raw_Effect_Matrix = effect_mat
    ))
  }
  
  grid_labeled <- grid
  grid_labeled$Label_Adjusted_ATE <- paste0(
    grid_labeled$Adjusted_ATE, " [", grid_labeled$Adjusted_ATE_CI_Lower, ", ", grid_labeled$Adjusted_ATE_CI_Upper, "]"
  )
  grid_labeled$Label_Abs_Change    <- as.character(grid_labeled$Abs_Percent_Change)
  grid_labeled$Label_Signed_Change <- as.character(grid_labeled$Percent_Change_with_direction)
  
  plot_data <- grid_labeled[!exclude_from_plot(grid_labeled), ]
  
  list(
    Observed_ATE      = obs_ATE,
    N_Simulations     = Simulation,
    CI_Level          = ci_level,
    Results           = grid,
    Raw_Effect_Matrix = effect_mat,
    
    graph_adjusted_ATE = make_text_plot(plot_data, "Label_Adjusted_ATE",
                                        "Adjusted ATE (with CI) | Prevalence of confounder: "),
    graph_percentage_change_ATE = make_text_plot(plot_data, "Label_Abs_Change",
                                                 "Percentage Change of ATE | Prevalence of confounder: "),
    graph_percentage_change_with_direction_ATE = make_text_plot(plot_data, "Label_Signed_Change",
                                                                "Percentage Change of ATE | Prevalence of confounder: "),
    
    graph_adjusted_ATE_heat_map = make_heatmap(plot_data, "Adjusted_ATE", "Label_Adjusted_ATE",
                                               "Adjusted ATE (with CI) | Prevalence of confounder: "),
    graph_percentage_change_ATE_heat_map = make_heatmap(plot_data, "Abs_Percent_Change", "Label_Abs_Change",
                                                        "Percentage Change of ATE | Prevalence of confounder: "),
    graph_percentage_change_with_direction_ATE_heat_map = make_heatmap(plot_data, "Percent_Change_with_direction",
                                                                       "Label_Signed_Change",
                                                                       "Percentage Change of ATE | Prevalence of confounder: ", diverging = TRUE)
  )
}
