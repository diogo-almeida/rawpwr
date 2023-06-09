% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{merge_with_sim_params}
\alias{merge_with_sim_params}
\title{Utility function to combine output of power simulations with the simulation
parameters.}
\usage{
merge_with_sim_params(sim_results, sim_params)
}
\arguments{
\item{sim_results}{A list. The output of a power simulation.}

\item{sim_params}{A list. The parameters of a given a power simulation.}
}
\value{
A data.frame with all the simulation parameters added to each
  simulation result
}
\description{
This function takes a simulation list output and adds the simulation
parameters and outputs a proper data.frame.
}
