#' Generate a covariance matrix from a correlation matrix and corresponding
#' standard deviations
#'
#' @description
#' Generates a covariance matrix based on a user-given correlation matrix.
#' This version assumes a single standard deviation for all pairs of
#' measurements.
#'
#' adapted from: https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
#'
#' @param cor_mat Correlation matrix.
#'
#' @param sd_var Vector of standard deviations
#'
#' @returns The corresponding covariance matrix.
#'
correlation2covariance <- function(cor_mat, sd_var = 50) {
  sweep(sweep(cor_mat, 1L, sd_var, "*"), 2L, sd_var, "*")
}

#' Calculates a covariance matrix from a set of correlation coefficients and
#' corresponding standard deviations
#'
#' @description
#' Calculates a covariance matrix based on a series of correlation coefficients
#' and standard deviations.

#' This version assumes a single correlation coefficient and a single standard
#' deviation for all pairs of measurements.
#'
#' @param rho Correlation coefficient, or a number between -1 and 1.
#'
#' @param std_dev Vector of standard deviations.
#'
#' @param ncond How many columns should the covariance matrix have. Default = 2.
#'
#' @returns The corresponding covariance matrix.
#'
calc_cov_mat <- function(rho, std_dev, ncond = 2) {
  corr_matrix <- matrix(rep(rho, times = ncond^2), ncol = ncond)
  diag(corr_matrix) <- 1
  cov_matrix <- correlation2covariance(cor_mat = corr_matrix, sd_var = std_dev)
  cov_matrix
}

#' Fast calculation of a paired t-test
#'
#' @description
#' Calculates a paired t-test more rapidly than the base R version, making it
#' suitable for simulations where the same operation needs to be repeated a very
#' large number of times.
#'
#' This version only calculates the paired t-test, based on a matrix with two
#' data columns.
#'
#' @param data_matrix A matrix with two columns upon which the paired t-test
#'   will be computed.
#'
#' @returns A list with the same statistics as the base R version, namely:
#'   t statistic, degrees of freedom, p-value, 95% confidence interval, mean
#'   difference and its standard error.
#'
fast_t <- function(data_matrix) {
  diff <- data_matrix[, 1] - data_matrix[, 2]
  n <- length(diff)
  deg_fre <- n - 1
  sd_diff <- stats::sd(diff)
  se_diff <- sd_diff / sqrt(n)
  mean_diff <- mean(diff)
  t_val <- mean_diff / se_diff
  t_crit <- stats::qt(.975, df = deg_fre)
  conf.int <- c(mean_diff - (t_crit * se_diff), mean_diff + (t_crit * se_diff))
  p_val <- stats::pt(-abs(t_val), df = deg_fre, lower.tail = TRUE) * 2 # change this for one-tail
  list(statistic = t_val, parameter = deg_fre, p.value = p_val, conf.int = conf.int,
       estimate = mean_diff, stderr = se_diff)
}

#' Convenient wrapper for fast paired t-test calculation
#'
#' @description
#' Convenient wrapper for fast paired t-test that outputs a sensible default
#' subset of information, and deals with the names of the columns.
#'
#' The use case for this function is in large simulations, to automatically
#' take care of handling column names and subsetting output.
#'
#' @param paired_data A matrix with two columns upon which the paired t-test
#'   will be computed.
#'
#' @param fast Boolean. Should the fast paired t-test function be used? Defaults
#'   to TRUE.
#'
#' @param use_names Boolean. Use existing column names in the output? Defaults
#'   to TRUE. If FALSE, column names are set to NA.
#'
#' @param info_out Vector of strings. Names the fields from the t-test function
#'   used in the output. Defaults to "estimate", "statistic" and "p.value".
#'
#' @returns A vector with the summary information described in "info_out".
#'
t_test_paired_data <- function(paired_data, fast = TRUE, use_names = TRUE,
                               info_out = c("estimate", "statistic", "p.value")) {
  if (fast) {
    test_result <- fast_t(paired_data)[info_out] |>
      unlist()
  } else {
    test_result <- stats::t.test(paired_data[, 1], paired_data[, 2],
                                 paired = TRUE)[info_out] |>
      unlist()
  }
  if(use_names) {
    true_means <- as.numeric(sub("C[0-9]+_", "", colnames(paired_data)))
    true_es <- true_means[1] - true_means[2]
    test_result <- c(true_es, test_result)
  } else {
    test_result <- c(NA, test_result)
  }
  names(test_result) <- c("true_es", "estimate", "statistic", "p.value")
  test_result
}

#' Calculates the main effects and the interaction in a 2x2 within-subjects
#' factorial design
#'
#' @description
#' Calculates the main effects and the interaction in a 2x2 within-subjects
#' factorial design by means of three paired t-tests. The interaction is
#' equivalent as a paired t-test of the difference between the main effects.
#'
#' It automatically calculates the adjusted p-values for the three tests, in
#' case it is needed.
#'
#' The use case for this function is in large simulations, to automatically
#' take care of handling column names and subsetting output.
#'
#' @param sim_results A matrix with four columns upon which the three
#'   paired t-tests will be computed.
#'
#' @param rename_cols Boolean. Automatically rename the columns? Defaults
#'   to TRUE, which should give sensible default names.
#'
#' @returns A vector with the summary information of the three paired t-tests.
#'
test_2x2_paired_data <- function(sim_results, rename_cols = TRUE) {
  main1 <- t_test_paired_data(sim_results[, 1:2])
  main2 <- t_test_paired_data(sim_results[, 3:4])
  diff1 <- sim_results[, 1] - sim_results[, 2]
  diff2 <- sim_results[, 3] - sim_results[, 4]
  diff_diff <- cbind(diff1, diff2)
  diff_es <- c(main1["true_es"], main2["true_es"])
  colnames(diff_diff) <- paste0("C", seq(diff_es), "_", diff_es)
  interaction <- t_test_paired_data(diff_diff)
  output <-  rbind(main1, main2, interaction)
  output[["p.adjust"]] <- stats::p.adjust(p = output$p.value, "holm")
  output[["condition"]] <- rownames(output)
  output
}

#' Calculates the type-S and type-M errors in a given power simulation.
#'
#' @description
#' Calculates the type-S (sign) and type-M (exaggeration ratio) errors for a
#' particular design on a power simulation, for both unadjusted and adjusted
#' p-values.
#'
#' @param sim_results A matrix with all the simulation results in the columns.
#'
#' @param alpha A double between 0 and 1, defining the significance level used in the
#'   simulations Defaults to .05.
#'
#' @param rename_cols A boolean. Whether or not to rename the columns. Defaults
#'   to TRUE. At the moment, this function requires this to be set to true.
#'
#' @returns A vector with the relevant information about the simulations (the
#'   simulated effect size, power (unadjusted and adjusted), type S and M error
#'   rates (unadjusted and adjusted)).
#'
calc_error_design <- function(sim_results, alpha = .05, rename_cols = TRUE) {
  stopifnot(rename_cols)
  col_labels <- colnames(sim_results)
  colnames(sim_results) <- sub("cond[0-9]+_", "", col_labels)

  sim_res_df <- as.data.frame(sim_results)

  sig_subset_unadjusted <- sim_res_df[which(sim_res_df$p.value < alpha), ]
  sig_subset_adjusted <- sim_res_df[which(sim_res_df$p.adjust < alpha), ]

  #sig_subset_unadjusted <- subset(as.data.frame(sim_results), p.value < alpha)
  #sig_subset_adjusted <- subset(as.data.frame(sim_results), p.adjust < alpha)
  sim_true_es <- mean(sig_subset_unadjusted$true_es)
  nsim <- nrow(sim_results)
  nsig_unadjusted <- nrow(sig_subset_unadjusted)
  nsig_adjusted <- nrow(sig_subset_adjusted)
  power_unadjusted <- nsig_unadjusted / nsim
  power_adjusted <- nsig_adjusted / nsim

  type_s_cases_unadjusted <- sig_subset_unadjusted[which(sig_subset_unadjusted$estimate < 0), ]
  type_s_cases_adjusted <- sig_subset_adjusted[which(sig_subset_adjusted$estimate < 0), ]
  type_s_unadjusted <- nrow(type_s_cases_unadjusted) / nsig_unadjusted
  type_s_adjusted <- nrow(type_s_cases_adjusted) / nsig_adjusted

  exaggeration_ratio_unadjusted <- mean(abs(sig_subset_unadjusted$estimate)) /
    mean(sig_subset_unadjusted$true_es)
  exaggeration_ratio_adjusted <- mean(abs(sig_subset_unadjusted$estimate)) /
    mean(sig_subset_adjusted$true_es)
  #exaggeration_ratio_unadjusted <- with(sig_subset_unadjusted,
  #                                      mean(abs(estimate)) / mean(true_es))
  #exaggeration_ratio_adjusted <- with(sig_subset_adjusted,
  #                                    mean(abs(estimate)) / mean(true_es))

  errors <- c(ES = sim_true_es,
              power_unadjusted = power_unadjusted,
              type_s_unadjusted = type_s_unadjusted,
              exaggeration_ratio_unadjusted = exaggeration_ratio_unadjusted,
              power_adjusted = power_adjusted,
              type_s_adjusted = type_s_adjusted,
              exaggeration_ratio_adjusted = exaggeration_ratio_adjusted)
  errors
}

#' Adjust p-values stored in a data.frame using the Holm-Bonferroni correction
#'
#' @description
#' Adjust p-values stored in a data.frame using the Holm-Bonferroni correction
#'
#' @param df A data.frame containing a column of p-values to be adjusted
#'
#' @param p_col A string defining the name of the column in the input data.frame
#'   that contains the p-values to be adjusted.
#'
#' @param method A string defining the name of the p-value adjustment method to
#'   be used. Defautls to Holm's method.
#'
#' @returns A data.frame identical to the original, with an added column
#'   containing the adjusted p-values.
#'
adjust_p <- function(df, p_col = "p.value", method = "holm") {
  df[["p.adjust"]] <- stats::p.adjust(p = df[[p_col]], method = method)
  df
}
