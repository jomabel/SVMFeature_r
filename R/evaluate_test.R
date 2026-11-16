#' Evaluate a Fixed Classifier on Test Data
#'
#' @description
#' Calculates held-out classification errors and metrics using an evaluated
#' solution's stored features, coefficients, and intermediate plane.
#'
#' @param solution A valid solution previously evaluated on training data.
#' @param data A nonempty data frame containing the selected predictors and target.
#' @param output The name of the target column, with numeric labels -1 and/or 1.
#'
#' @details
#' Preprocessing must be fitted on training data and applied unchanged to test.
#' This function does not normalize data, reconstruct planes, change training
#' metrics or objectives, or modify Pareto ranks. Feature columns are matched
#' by name in the solution's stored order; additional columns are ignored.
#' Scores are obtained through the solution's \code{predict()} method. Legacy
#' evaluated lists are given class \code{solution} locally for prediction.
#'
#' Class 1 is positive, including decision scores of zero. Confusion matrix
#' rows are actual classes and columns are predictions, ordered -1, 1.
#' Undefined precision and recall are NA/NaN, as in training; undefined F1 is
#' zero. Kappa is NA when undefined, and AUC is NA for a single true class.
#' Non-finite decision scores are rejected with an error.
#'
#' @return A separate list containing \code{n_samples}, \code{n_errors},
#'   \code{error_rate} (\code{n_errors / n_samples}, from 0 to 1),
#'   \code{fn}, \code{fp}, \code{accuracy}, \code{precision}, \code{recall},
#'   \code{f1}, \code{kappa}, \code{auc}, \code{confusion_matrix},
#'   \code{predictions}, and \code{decision_scores}. Predictions and scores
#'   are named vectors in test-row order, with names from the data's row names.
#'
#' @export
#' @examples
#' train <- data.frame(x = c(-2, -1, 1, 2), y = c(-1, -1, 1, 1))
#' sol <- solution(0, 1, "MC")
#' sol$features <- "x"
#' sol$vectors <- list(1L, 4L)
#' sol$plane_coord <- 1
#' sol <- evaluate_solution(sol, train, "y")
#' test <- data.frame(x = c(-3, -0.5, 0, 0.5, 3), y = c(-1, 1, -1, 1, 1))
#' metrics <- evaluate_test(sol, test, "y")
#' metrics$error_rate
#' metrics$confusion_matrix
evaluate_test <- function(solution, data, output) {
  finite_numeric <- function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x))
  }
  if (!is.list(solution) || !isTRUE(solution$successful_evaluation)) {
    stop("A valid evaluated solution is required before evaluate_test().", call. = FALSE)
  }
  prediction_object <- solution
  if (!inherits(prediction_object, "solution")) {
    class(prediction_object) <- c("solution", class(prediction_object))
  }
  scores <- stats::predict(prediction_object, newdata = data, type = "score")
  data <- as.data.frame(data)
  features <- unlist(solution$features, use.names = FALSE)
  if (!is.character(output) || !is.null(dim(output)) || length(output) != 1L || is.na(output) ||
      !output %in% names(data) || output %in% features) {
    stop("output must name a target column distinct from the selected features.", call. = FALSE)
  }
  if (!finite_numeric(data[[output]])) {
    stop("Test output must contain finite real numeric values.", call. = FALSE)
  }
  y_true <- data[[output]]
  if (!all(y_true %in% c(-1, 1))) {
    stop("Test output must contain only numeric classes -1 and 1.", call. = FALSE)
  }

  predictions <- .prediction_classes(scores)
  metrics <- calculate_metrics(list(), y_true, predictions, scores)
  if (length(unique(y_true)) < 2L) metrics$auc <- NA_real_
  fn <- unname(metrics$conf_matrix["1", "-1"])
  fp <- unname(metrics$conf_matrix["-1", "1"])
  names(predictions) <- names(scores) <- rownames(data)
  list(n_samples = nrow(data), n_errors = fn + fp,
       error_rate = (fn + fp) / nrow(data), fn = fn, fp = fp,
       accuracy = metrics$accuracy, precision = metrics$precision,
       recall = metrics$recall, f1 = metrics$f1, kappa = metrics$kappa,
       auc = metrics$auc, confusion_matrix = metrics$conf_matrix,
       predictions = predictions, decision_scores = scores)
}
