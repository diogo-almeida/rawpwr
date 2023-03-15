#' Applies rbind() to a list
#'
#' @description
#' `rbind_list()` takes a list as an argument and applies `rbind()` to it to
#' create a data.table or matrix element.
#'
#' With the introduction of native pipes in base R, this functions allows one
#' one to pipe a list that is being processed in a pipe context and pass it to
#' rbind using to `do.call()`
#'
#' @param input_list Input list, generally the product of a pipe context.
#'
#' @returns A new data.table or matrix that pulls together all the items from
#'  the input list.
#'
#' @seealso [rbind()], which this function wraps.
#'
#' @examples
#' list(a = rnorm(10), b = rnorm(10)) |>
#'   lapply(mean) |>
#'   rbind_list()
#' @export
rbind_list <- function(input_list) {
  do.call(rbind, input_list)
}


#' Breaks down a matrix into a list where each item is an arbitrary number of
#' columns from the original matrix.
#'
#' @description
#' This function takes a matrix and breaks it into a user-specified number of
#' columns, where each group becomes a list item.
#'
#' For the time being, this is a simple helper function, and does not perform
#' any checks on the input, so the user needs to make sure the input is correct,
#' otherwise unpredictable behavior can occur.
#'
#' @param data_matrix Input matrix.
#'
#' @param group_by How many columns are supposed to be in each group that will
#'   become a list item. Defaults to 2
#'
#' @returns A list with all the columns from the input matrix rearranged as
#'   items, each with the user-specified number of columns
#'
matcols2lists <- function(data_matrix, group_by = 2) {
  ncols <- ncol(data_matrix)
  stopifnot(!(ncols%%group_by))
  start_indices <- seq(from = 1, to = ncols, by = group_by)
  n_conditions <- length(start_indices)
  grouped_conditions <- vector(mode = "list", length = n_conditions)
  for (first_item in start_indices) {
    grouped_conditions[[which(start_indices == first_item)]] <- data_matrix[, seq(first_item, first_item + group_by - 1)]
  }
  grouped_conditions
}


#' Utility function to organize output of power simulations.
#'
#' @description
#' This function takes a simulation list output and organizes it for further
#' processing. It is mostly a convenient way of shortening the data processing
#' pipeline.
#'
#' @param simdata_list The output of a power simulation..
#'
#' @returns A list with all the columns adequately labeled.
organize_output <- function(simdata_list) {
  step_1 <- lapply(simdata_list, FUN = t)
  n_conds <- ncol(step_1[[1]])
  n_measures <- nrow(step_1[[1]])
  measure_labels <- rownames(step_1[[1]])
  c_lbl <- paste0("cond", rep(1:n_conds, each = n_measures), "_", rep(measure_labels, n_conds))
  step_2 <- step_1 |>
    rapply(as.vector, how = "list") |>
    rbind_list()
  colnames(step_2) <- c_lbl
  step_2 |>
    matcols2lists(group_by = n_measures)
}

#' Utility function to combine output of power simulations with multiple
#' parameters.
#'
#' @description
#' This function takes a simulation list output and adds simulation and
#' condition numbers and outputs a proper data.frame for data exploration
#'
#' @param sim_list The output of a power simulation.
#'
#' @returns A list with all the columns adequately labeled.
combine_sims <- function(sim_list) {
  n_conds <- nrow(sim_list[[1]])
  n_sims <- length(sim_list)
  sim_df <- rbind_list(sim_list)
  sim_df[["simulation"]]  <- paste0("simulation_", rep(1:n_sims, each = n_conds))
  sim_df[["condition"]] <- paste0("condition_", rep(1:n_conds, times = n_sims))
  sim_df
}
