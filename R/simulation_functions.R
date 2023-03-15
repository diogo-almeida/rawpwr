#' Function factory to generated a simulation function based on a user-specified
#' set of effect sizes
#'
#' @description
#' This is a function factory, i.e., a closure that takes a set of effect sizes
#' specified by the user and outputs a function that the user can parametrize
#' in terms of correlations, standard deviations and sample sizes.
#'
#' The use case in in large scale power simulations.
#'
#' @param es Vector representing the effects sizes to be simulated in the
#'   generated function
#'
#' @param base_mu A double. If one desires to use more realistic numbers in the
#'   simulation. Defaults to 600, a sensible mean reaction time for lexical
#'   decision experiments.
#'
#' @returns A function that closes over the effects sizes and can be
#'   further parametrized by the user.
#' @export
simulate_data_from_es <- function(es, base_mu = 600) {
  function(rho, std_dev, nsubj) {
    condition_means <- as.vector(es + base_mu)
    list(rho, std_dev)
    cov_mat <- calc_cov_mat(rho, std_dev, ncond = length(es))  # one rho all conds
    simulation <- MASS::mvrnorm(n = nsubj, mu = condition_means, Sigma = cov_mat,
                                empirical = FALSE)
    colnames(simulation) <- paste0("C", seq(condition_means), "_", condition_means)
    simulation
  }
}

#' Performs a power simulation based on user-specified parameters
#'
#' @description
#' Simulates data using user-specified criteria and design, and evaluate the
#' rate of detection (power) and types of errors (S and M) of a given
#' experimental design.
#'
#' The use case in in large scale power simulations.
#'
#' @param sim_params A list containing three fields: the correlations (rho),
#'   the standard deviations (std_dev) and sample size (nsubj).
#'
#' @param sim_func A function. This function produces the simulated datasets
#'   incorporating the user-specified effect sizes. This function is produced by
#'   the `simulate_data_from_es()` function factory.
#'
#' @param Nsim An integer. The number of simulations per set of parameters to be
#'   conducted.
#'
#' @param alpha A double between 0 and 1, defining the significance level used in the
#'   simulations Defaults to .05.
#'
#' @returns A data.frame with all the summary information for each set of
#'   parameters specified by the user.
#' @export
do_sim <- function(sim_params, sim_func, Nsim = 5000, alpha = 0.05){
  replicate(Nsim, do.call(sim_func, sim_params), simplify = FALSE) |>
    lapply(FUN = matcols2lists, group_by = 2) |>
    rapply(t_test_paired_data, how = "list") |>
    lapply(FUN = rbind_list) |>
    lapply(FUN = adjust_p) |>
    organize_output() |>
    lapply(FUN = calc_error_design, alpha = alpha) |>
    rbind_list() |>
    data.frame(sim_params)
}

#' Performs a power simulation based on user-specified parameters for a 2x2
#' within-subject factorial design.
#'
#' @description
#' Simulates data using user-specified criteria and design, and evaluate the
#' rate of detection (power) and types of errors (S and M) of a given
#' experimental design.
#'
#' The use case in in large scale power simulations.
#'
#' @param sim_params A list containing three fields: the correlations (rho),
#'   the standard deviations (std_dev) and sample size (nsubj).
#'
#' @param sim_func A function. This function produces the simulated datasets
#'   incorporating the user-specified effect sizes. This function is produced by
#'   the `simulate_data_from_es()` function factory.
#'
#' @param Nsim An integer. The number of simulations per set of parameters to be
#'   conducted.
#'
#' @param alpha A double between 0 and 1, defining the significance level used in the
#'   simulations Defaults to .05.
#'
#' @returns A data.frame with all the summary information for each set of
#'   parameters specified by the user.
#' @export
do_sim_2x2 <- function(sim_params, sim_func, Nsim = 5000, alpha = 0.05) {
  replicate(Nsim, do.call(sim_func, sim_params), simplify = FALSE) |>
    lapply(FUN = matcols2lists, group_by = 4) |>
    rapply(test_2x2_paired_data, how = "list" )  #|>
    #unlist(recursive = FALSE) |>
    #rbind_list() |>
    #split(~condition) |>
    #lapply(FUN = calc_error_design, alpha = alpha) |>
    #rbind_list() |>
    #data.frame(sim_params)
}

