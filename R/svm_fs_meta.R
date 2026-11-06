utils::globalVariables(c("FRONT"))
utils::globalVariables(c("DIST"))
utils::globalVariables(c("EPS"))
utils::globalVariables(c("MCPOS"))
utils::globalVariables(c("MCNEG"))

#' svm_fs_meta Class (S3 Implementation)
#'
#' Configures multi-objective support vector machine (SVM) optimization for
#' binary classification with a fixed number of selected features.
#'
#' @param data A non-empty data frame with unique, non-empty column names.
#' @param inputs A non-empty character vector of unique predictor column names.
#'   Selected predictors must be real numeric vectors with finite values and
#'   must not include the response column.
#' @param output A single response column name. The response must be numeric,
#'   contain only -1 and 1, and have at least two observations of each class.
#' @param pop_size A single finite integer of at least 3. The population size.
#' @param num_fea A single positive integer, no larger than the number of inputs.
#'   Each candidate selects exactly this many features throughout optimization.
#' @param n_iter A single positive finite integer. Number of iterations
#'   (default: 10).
#' @param max_time A single positive finite number. Approximate optimization
#'   time budget in seconds when \code{mode = "time"} (default: 300).
#'   Initialization is excluded. The budget is checked before each iteration;
#'   an iteration already started finishes, so total runtime may exceed it.
#' @param mode Character. "iters" for iteration-based execution, "time" for time-based execution.
#' @param objective Character. "DIST-EPS" maximizes the distance between margin
#'   planes and minimizes the sum of margin-violation distances. "MC" minimizes
#'   false-negative and false-positive counts from the intermediate classifier.
#' @param clones A single 0/1 flag: 0 forbids clone solutions and 1 allows them
#'   (default: 0).
#'   Logical values are also accepted.
#' @param p_mutation A single finite mutation probability in \code{[0, 1]}
#'   (default: 0.7).
#' @param p_mut_ind A single finite individual mutation probability in \code{[0, 1]}
#'   (default: 0.4).
#' @param p_mut_fea A single finite feature mutation probability in \code{[0, 1]}
#'   (default: 0.4).
#' @param p_mut_coord A single finite coordinate mutation probability in \code{[0, 1]}
#'   (default: 0.2).
#' @param mut_coord A single 0/1 flag for the coordinate mutation strategy
#'   (default: 0). Logical values are also accepted.
#' @details
#'   The current evolutionary operators require at least three individuals and
#'   two observations in each response class. Data-frame subclasses are converted
#'   locally to a base data frame. Columns outside the selected inputs and output
#'   are not used and need not be numeric or complete.
#'
#'   The feature count is fixed by \code{num_fea}; candidate solutions differ in
#'   which features they select, rather than in how many they select.
#'
#'   In "MC" mode, class 1 is positive. The third objective is the false-negative
#'   count and the fourth is the false-positive count. These use the same
#'   intermediate-plane predictions as the reported confusion matrix; scores
#'   of zero are assigned to class 1.
#' @return An S3 object of class "svm_fs_meta" containing the input data,
#'   predictor and response names, optimization settings, and a configured
#'   \code{population} element. Call \code{run.svm_fs_meta()} to generate and
#'   optimize its candidate solutions.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#'
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize n
svm_fs_meta <- function(data, inputs, output, pop_size, num_fea,
                       n_iter = 10, max_time = 300, mode = "iters", objective="DIST-EPS",
                       clones = 0, p_mutation = 0.7, p_mut_ind = 0.4,
                       p_mut_fea = 0.4, p_mut_coord = 0.2, mut_coord = 0) {

  is_scalar_number <- function(value) {
    is.numeric(value) && !is.complex(value) && is.null(dim(value)) &&
      length(value) == 1L && is.finite(value)
  }
  is_scalar_choice <- function(value, choices) {
    is.character(value) && is.null(dim(value)) && length(value) == 1L &&
      !is.na(value) && value %in% choices
  }
  is_finite_numeric_vector <- function(value) {
    is.numeric(value) && !is.complex(value) && is.null(dim(value)) &&
      all(is.finite(value))
  }

  if (!is.data.frame(data)) {
    stop('"data" must be a data frame.', call. = FALSE)
  }
  data <- as.data.frame(data)
  if (nrow(data) == 0L || ncol(data) == 0L) {
    stop('"data" must contain at least one row and one column.', call. = FALSE)
  }
  if (anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop('"data" must have unique, non-empty column names.', call. = FALSE)
  }
  if (!is_scalar_choice(mode, c("iters", "time"))) {
    stop('"mode" must be either "iters" or "time".', call. = FALSE)
  }
  if (!is.character(inputs) || !is.null(dim(inputs)) || length(inputs) == 0L ||
      anyNA(inputs) || any(!nzchar(inputs)) || anyDuplicated(inputs) ||
      !all(inputs %in% names(data))) {
    stop('"inputs" must be a non-empty character vector of unique existing column names in "data".',
         call. = FALSE)
  }
  if (!is_scalar_choice(output, names(data))) {
    stop('"output" must name a single column in "data".', call. = FALSE)
  }
  if (output %in% inputs) {
    stop('"inputs" must not include "output".', call. = FALSE)
  }

  counts <- list(pop_size = pop_size, num_fea = num_fea, n_iter = n_iter)
  minimums <- c(pop_size = 3, num_fea = 1, n_iter = 1)
  for (parameter in names(counts)) {
    value <- counts[[parameter]]
    if (!is_scalar_number(value) || value != floor(value) || value < minimums[[parameter]]) {
      stop(sprintf('"%s" must be a finite integer greater than or equal to %d.',
                   parameter, minimums[[parameter]]), call. = FALSE)
    }
  }
  if (num_fea > length(inputs)) {
    stop('"num_fea" must not exceed the number of "inputs".', call. = FALSE)
  }
  if (!is_scalar_number(max_time) || max_time <= 0) {
    stop('"max_time" must be a positive finite number.', call. = FALSE)
  }
  if (!is_scalar_choice(objective, c("DIST-EPS", "MC"))) {
    stop('"objective" must be either "DIST-EPS" or "MC".', call. = FALSE)
  }

  probabilities <- list(p_mutation = p_mutation, p_mut_ind = p_mut_ind,
                        p_mut_fea = p_mut_fea, p_mut_coord = p_mut_coord)
  for (parameter in names(probabilities)) {
    value <- probabilities[[parameter]]
    if (!is_scalar_number(value) || value < 0 || value > 1) {
      stop(sprintf('"%s" must be a finite number between 0 and 1.', parameter),
           call. = FALSE)
    }
  }
  flags <- list(clones = clones, mut_coord = mut_coord)
  for (parameter in names(flags)) {
    value <- flags[[parameter]]
    is_logical_scalar <- is.logical(value) && is.null(dim(value)) &&
      length(value) == 1L && !is.na(value)
    if (!(is_scalar_number(value) || is_logical_scalar) || !value %in% c(0, 1)) {
      stop(sprintf('"%s" must be a single 0/1 or logical flag.', parameter),
           call. = FALSE)
    }
  }

  if (!all(vapply(data[inputs], is_finite_numeric_vector, logical(1)))) {
    stop('"inputs" must select real numeric columns containing only finite values.',
         call. = FALSE)
  }
  response <- data[[output]]
  if (!is_finite_numeric_vector(response) || !all(response %in% c(-1, 1)) ||
      !all(c(-1, 1) %in% response)) {
    stop('"output" must contain both numeric classes -1 and 1, with no other values.',
         call. = FALSE)
  }
  if (sum(response == -1) < 2L || sum(response == 1) < 2L) {
    stop('"output" must contain at least two observations in each class.', call. = FALSE)
  }

  # Create the S3 object (a list with attributes)
  object <- list(
    data = data,
    inputs = inputs,
    output = output,
    pop_size = pop_size,
    num_fea = num_fea,
    n_iter = n_iter,
    max_time = max_time,
    mode = mode,
    objective = objective,
    population = population(pop_size = pop_size, num_features = num_fea,
                            clones = clones, p_mutation = p_mutation,
                            p_mut_ind = p_mut_ind, p_mut_fea = p_mut_fea,
                            p_mut_coord = p_mut_coord, mut_coord = mut_coord,
                            objective = objective)
  )

  class(object) <- "svm_fs_meta"
  return(object)
}

#' Run the svm_fs_meta optimization process
#' @param object An object of class "svm_fs_meta".
#' @return An updated S3 object of class "svm_fs_meta", retaining the input data
#'   and optimization settings. Its \code{population} element is an object of
#'   class "population" containing:
#'   \describe{
#'     \item{solution_list}{A list of the retained candidate solutions. Each
#'       solution is a list containing its selected feature names in
#'       \code{features}, normal coefficients in \code{plane_coord}, and three
#'       intercepts in \code{plane_term_b}. The coefficients follow the order of
#'       \code{features}; the intercepts correspond to the negative-class margin
#'       plane, positive-class margin plane, and intermediate classifier plane.
#'       Its \code{objective} vector contains distance, epsilon, false negatives,
#'       and false positives, in that order.}
#'     \item{df_solutions}{A data frame with one row per retained solution,
#'       ordered by front. The numeric \code{FRONT} column gives the front number;
#'       the other columns currently contain lists of solution values.}
#'     \item{fronts}{The front membership information produced by non-dominated
#'       sorting.}
#'   }
#' @details
#'   In time mode, \code{object$max_time} is an approximate budget for the
#'   optimization loop. Timing starts after the initial population has been
#'   generated, sorted into fronts, and summarized. The budget is checked before
#'   each iteration, and an iteration already started runs to completion. Total
#'   runtime therefore includes initialization and may exceed the requested
#'   budget; this is not a strict timeout.
#'
#'   Front 1 contains the non-dominated solutions in the final population and
#'   provides an approximation to the problem's Pareto front. Both
#'   \code{solution_list} and \code{df_solutions} include all retained fronts.
#'
#'   Each solution also stores \code{accuracy}, \code{precision}, \code{recall},
#'   \code{f1}, \code{kappa}, \code{auc}, and \code{conf_matrix}. These metrics
#'   are calculated on the data used for optimization, with class 1 treated as
#'   positive. The confusion matrix has true classes in rows and predicted
#'   classes in columns, ordered as -1, 1.
#'   The historical summary-table and CSV columns \code{MCPOS} and \code{MCNEG}
#'   contain false-negative and false-positive counts, respectively.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.svm_fs_meta(obj)
#' front <- res$population$df_solutions
#' front <- front[front$FRONT == 1, , drop = FALSE]
#' front_solutions <- Filter(function(s) s$front == 1,
#'                           res$population$solution_list)
#' selected <- front_solutions[[1]]
#' selected$features
#' selected$plane_coord
#' selected$plane_term_b
#' selected[c("accuracy", "precision", "recall", "f1", "kappa", "auc",
#'            "conf_matrix")]
run.svm_fs_meta <- function(object) {

  object$population <- generate_initial_population(object$population, object$data,
                                                   object$inputs, object$output)

  object$population <- fnds(object$population)
  object <- update_df_solutions.svm_fs_meta(object)

  init_time <- Sys.time()

  run_iteration <- function() {
    object$population <- new_population(object$population, object$data,
                                        object$inputs, object$output)
    object$population <- fnds(object$population)

    reduced_population <- reduce_population(object$population)
    reduced_population <- fnds(reduced_population)

    object$population <- reduced_population
    object <- update_df_solutions.svm_fs_meta(object)
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

#' Draw the Pareto Front of an svm_fs_meta object
#' @param x An object of class "svm_fs_meta".
#' @param ... Additional arguments (not used).
#' @returns A \code{ggplot} object displaying front 1, returned invisibly after
#'   printing the plot. Returns \code{NULL} if there are no solutions to plot
#'   or the objective is unsupported.
#' @details
#'   For "MC", the horizontal axis shows false negatives (\code{MCPOS}) and
#'   the vertical axis shows false positives (\code{MCNEG}).
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal element_line
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.svm_fs_meta(obj)
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

  if (x$objective == "DIST-EPS") {
    front_population1$DIST <- as.numeric(front_population1$DIST)
    front_population1$EPS <- as.numeric(front_population1$EPS)

    xaxis_label <- 'Distance'
    yaxis_label <- 'Epsilon'
    xvalues <- front_population1$DIST
    yvalues <- front_population1$EPS

  } else if (x$objective == "MC") {
    front_population1$MCPOS <- as.numeric(front_population1$MCPOS)
    front_population1$MCNEG <- as.numeric(front_population1$MCNEG)

    xaxis_label <- 'False Negatives'
    yaxis_label <- 'False Positives'
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
#' @param object An object of class "svm_fs_meta".
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
#' obj <- svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' res <- run.svm_fs_meta(obj)
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

#' Update df_solutions for an svm_fs_meta object
#' @param object An object of class "svm_fs_meta".
#' @return The updated svm_fs_meta object.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' obj <- svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
#'                   pop_size = 4, num_fea = 2, n_iter = 1)
#' obj$population <- generate_initial_population(obj$population, obj$data, obj$inputs, obj$output)
#' obj$population <- fnds(obj$population)
#' obj <- update_df_solutions.svm_fs_meta(obj)
update_df_solutions.svm_fs_meta <- function(object) {
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
