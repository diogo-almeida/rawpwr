% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{do_sim_2x2}
\alias{do_sim_2x2}
\title{Performs a power simulation based on user-specified parameters for a 2x2
within-subject factorial design.}
\usage{
do_sim_2x2(sim_params, sim_func, Nsim = 5000, alpha = 0.05)
}
\arguments{
\item{sim_params}{A list containing three fields: the correlations (rho),
the standard deviations (std_dev) and sample size (nsubj).}

\item{sim_func}{A function. This function produces the simulated datasets
incorporating the user-specified effect sizes. This function is produced by
the `simulate_data_from_es()` function factory.}

\item{Nsim}{An integer. The number of simulations per set of parameters to be
conducted.}

\item{alpha}{A double between 0 and 1, defining the significance level used in the
simulations Defaults to .05.}
}
\value{
A data.frame with all the summary information for each set of
  parameters specified by the user.
}
\description{
Simulates data using user-specified criteria and design, and evaluate the
rate of detection (power) and types of errors (S and M) of a given
experimental design.

The use case in in large scale power simulations.
}
