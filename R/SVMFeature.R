utils::globalVariables(c("FRONT"))
utils::globalVariables(c("DIST"))
utils::globalVariables(c("EPS"))
utils::globalVariables(c("MCPOS"))
utils::globalVariables(c("MCNEG"))

library(dplyr)

#' SVMFeature Class (S3 Implementation)
#'
#' Implements a feature selection algorithm using Support Vector Machines (SVM).
#'
#' @param data DataFrame. The input data.
#' @param inputs Character vector. The names of the input features.
#' @param output Character. The name of the output variable.
#' @param costs Numeric vector. The costs associated with the features.
#' @param pop_size Numeric. The population size.
#' @param num_fea Numeric. The number of features.
#' @param n_iter Numeric. The number of iterations (default: 10).
#' @param max_time Numeric. The maximum execution time (default: 300 seconds).
#' @param mode Character. "iters" for iteration-based execution, "time" for time-based execution.
#' @param objective Character. "distance-epsilon" or "confusion-matrix".
#' @param num_obj Numeric. Number of objectives (default: 2).
#' @param clones Numeric. Number of clones (default: 0).
#' @param p_mutation Numeric. Mutation probability (default: 0.7).
#' @param p_mut_ind Numeric. Individual mutation probability (default: 0.4).
#' @param p_mut_fea Numeric. Feature mutation probability (default: 0.4).
#' @param p_mut_coord Numeric. Mutation coordinates probability (default: 0.2).
#' @param mut_coord Numeric. Mutation coordinates (default: 0).
#' @return An S3 object of class "SVMFeature".
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 1),
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#'
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize n
SVMFeature <- function(data, inputs, output, costs, pop_size, num_fea,
                       n_iter = 10, max_time = 300, mode = "iters", objective="distance-epsilon",
                       num_obj = 2, clones = 0, p_mutation = 0.7, p_mut_ind = 0.4,
                       p_mut_fea = 0.4, p_mut_coord = 0.2, mut_coord = 0) {

  if (!is.data.frame(data)) stop('Error: "data" must be a data frame')
  if (!mode %in% c("iters", "time")) stop('Error: "mode" must be either "iters" or "time"')
  if (!is.character(inputs) || !all(inputs %in% colnames(data))) stop('Error: "inputs" must be a character vector of existing column names in "data"')
  if (!is.character(output) || length(output) != 1 || !(output %in% colnames(data))) stop('Error: "output" must be a single column name in "data"')
  if (!is.numeric(max_time) || max_time <= 0) stop('Error: "max_time" must be a positive number')
  if (!is.numeric(n_iter) || n_iter <= 0) stop('Error: "n_iter" must be a positive number')
  if (!is.numeric(pop_size) || pop_size <= 0) stop('Error: "pop_size" must be a positive number')
  if (!is.numeric(num_fea) || num_fea <= 0) stop('Error: "num_fea" must be a positive number')
  if (!objective %in% c("distance-epsilon", "distance-epsilon-costs",
                        "confusion-matrix", "confusion-matrix-costs"))
    stop('Error: "objective" must be either "distance-epsilon",
         "distance-epsilon-costs", "confusion-matrix" or
         "confusion-matrix-costs"')

  # Normalization function
  scaler <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(0.5, length(x)))
    return((x - rng[1]) / (rng[2] - rng[1]))
  }

  # Normalize input data
  norm_data <- as.data.frame(lapply(data[, inputs], scaler))
  norm_data[[output]] <- data[[output]]

  # Create the S3 object (a list with attributes)
  object <- list(
    data = norm_data,
    inputs = inputs,
    output = output,
    costs = costs,
    pop_size = pop_size,
    num_fea = num_fea,
    n_iter = n_iter,
    max_time = max_time,
    mode = mode,
    objective = objective,
    population = Population(pop_size = pop_size, num_features = num_fea,
                            num_obj = num_obj, clones = clones, p_mutation = p_mutation,
                            p_mut_ind = p_mut_ind, p_mut_fea = p_mut_fea,
                            p_mut_coord = p_mut_coord, mut_coord = mut_coord,
                            objective = objective)
  )

  class(object) <- "SVMFeature"
  return(object)
}

#' Run the SVMFeature optimization process
#' @param object An object of class "SVMFeature".
#' @return The final population after the optimization process.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 1),
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.SVMFeature(obj)
run.SVMFeature <- function(object) {

  object$population <- generate_initial_population(object$population, object$data,
                                                   object$inputs, object$output,
                                                   object$costs)

  object$population <- fnds(object$population)
  object <- update_df_solutions.SVMFeature(object)

  init_time <- Sys.time()

  run_iteration <- function() {
    object$population <- new_population(object$population, object$data,
                                        object$inputs, object$output, object$costs)
    object$population <- fnds(object$population)

    reduced_population <- reduce_population(object$population)
    reduced_population <- fnds(reduced_population)

    object$population <- reduced_population
    object <- update_df_solutions.SVMFeature(object)
    fronts_count <- object$population$df_solutions %>% dplyr::group_by(FRONT) %>% dplyr::summarize(count = dplyr::n())
    return(object)
  }

  if (object$mode == "iters") {
    for (iter in seq_len(object$n_iter)) {
      message(sprintf("ITERATION %d: %f seconds", iter, as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
      object <- run_iteration()
    }
  } else if (object$mode == "time") {
    i <- 0
    while (as.numeric(difftime(Sys.time(), init_time, units = "secs")) < object$max_time) {
      message(sprintf("ITERATION %d: %f seconds", i, as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
      object <- run_iteration()
      i <- i + 1
    }
  } else {
    stop("Invalid execution mode. Use 'iters' or 'time'.")
  }

  print_population(object$population)

  return(object)
}

#' Draw the Pareto Front of an SVMFeature object
#' @param x An object of class "SVMFeature".
#' @param ... Additional arguments (not used).
#' @returns An object of class `ggplot` representing the selected solution.
#'    The plot summarizes the strcture of characteristics of the solution
#'    generated by the optimization procedure.
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal element_line
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 1),
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.SVMFeature(obj)
#' draw_solution(res)
draw_solution <- function(x, ...) {
  if (is.null(x$population$df_solutions)) {
    message("No solutions to plot.")
    return(NULL)
  }

  front_population1 <- dplyr::filter(x$population$df_solutions, FRONT == 1)

  if (nrow(front_population1) == 0) {
    message("There is no solution in front 1 to plot.")
    return(NULL)
  }

  if (x$objective == "distance-epsilon") {
    front_population1$DIST <- as.numeric(front_population1$DIST)
    front_population1$EPS <- as.numeric(front_population1$EPS)

    xaxis_label <- 'Distance'
    yaxis_label <- 'Epsilon'
    xvalues <- front_population1$DIST
    yvalues <- front_population1$EPS

  } else if (x$objective == "confusion-matrix") {
    front_population1$MCPOS <- as.numeric(front_population1$MCPOS)
    front_population1$MCNEG <- as.numeric(front_population1$MCNEG)

    xaxis_label <- 'False Positives'
    yaxis_label <- 'False Negatives'
    xvalues <- front_population1$MCPOS
    yvalues <- front_population1$MCNEG
  } else {
    message("Objective not supported for plotting.")
    return(NULL)
  }

  p <- ggplot2::ggplot(front_population1, ggplot2::aes(x = xvalues, y = yvalues)) +
    ggplot2::geom_point(color = 'blue') +
    ggplot2::labs(x = xaxis_label, y = yaxis_label, title = 'Non-Dominated Solutions') +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid', colour = "gray"))

  print(p)
  return(invisible(p))
}

#' Save the Pareto front results to a CSV file
#' @param object An object of class "SVMFeature".
#' @param output_dir Character. The directory to save the results.
#'    This argument must be supplied explicitly.
#' @param dataset_name Character. The name of the dataset (default: "dataset").
#' @returns A character string with the path to the saved file.
#'    This value indicates where the results were written.
#' @importFrom dplyr filter mutate across
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 1),
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.SVMFeature(obj)
#' temp_dir <- tempdir()
#' save_results(res, output_dir = temp_dir, dataset_name = "test")
save_results <- function(object, output_dir, dataset_name = "dataset") {
  if (missing(output_dir) || is.null(output_dir) || !nzchar(output_dir)) {
    stop("Please provide 'output_dir' explicitly.")
  }

  if (is.null(object$population$df_solutions)) {
    stop("No solutions found in the object.")
  }

  # Extract results from the Pareto front
  pareto_front <- object$population$df_solutions %>% dplyr::filter(FRONT == 1)

  if (nrow(pareto_front) == 0) {
    message("No solutions in front 1 to save.")
    return(invisible(NULL))
  }

  # Convert list-type columns to text
  pareto_front <- pareto_front %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.list), ~ sapply(., function(x) paste(x, collapse = ","))))

  final_output_dir <- file.path(output_dir, dataset_name, object$objective)
  if (!dir.exists(final_output_dir)) dir.create(final_output_dir, recursive = TRUE)

  iter_val <- if (object$mode == "iters") object$n_iter else object$max_time
  base_filename <- paste0(object$pop_size, "_", object$mode, "_", iter_val, ".csv")
  output_file <- file.path(final_output_dir, base_filename)

  counter <- 1
  while (file.exists(output_file)) {
    output_file <- file.path(final_output_dir, paste0(object$pop_size, "_", object$mode, "_", iter_val, "_", counter, ".csv"))
    counter <- counter + 1
  }

  utils::write.csv(pareto_front, output_file, row.names = FALSE)
  message(sprintf("Pareto front saved in: %s", output_file))
  return(invisible(output_file))
}

#' Update df_solutions for an SVMFeature object
#' @param object An object of class "SVMFeature".
#' @return The updated SVMFeature object.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 1),
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' obj$population <- generate_initial_population(obj$population, obj$data, obj$inputs, obj$output, obj$costs)
#' obj$population <- fnds(obj$population)
#' obj <- update_df_solutions.SVMFeature(obj)
update_df_solutions.SVMFeature <- function(object) {
  df_solutions <- do.call(rbind, lapply(object$population$solution_list, function(solution) {
    solution_dict <- to_dict(solution)
    solution_dict <- lapply(solution_dict, function(item) if (is.list(item)) unlist(item) else item)
    as.data.frame(t(solution_dict), stringsAsFactors = FALSE)
  }))

  df_solutions$FRONT <- as.numeric(as.character(df_solutions$FRONT))
  df_solutions <- df_solutions[order(df_solutions$FRONT), ]
  object$population$df_solutions <- df_solutions

  return(object)
}
