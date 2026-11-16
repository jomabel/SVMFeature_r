#' Predict Classes or Decision Scores from a Solution
#'
#' @description
#' Applies a solution's fixed classifier to one or more new observations without
#' requiring their true classes or refitting the classifier.
#'
#' @param object An evaluated object of class \code{solution}.
#' @param newdata A nonempty data frame with the selected predictors. A single
#'   observation must be supplied as a one-row data frame. Extra columns are ignored.
#' @param type Either \code{"class"} (default), returning -1 or 1, or
#'   \code{"score"}, returning signed intermediate-plane decision scores.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' Feature columns are matched by name in the stored coefficient order. Apply
#' the transformation fitted on training data before predicting; the method
#' does not normalize data or change the solution. Class 1 is positive and
#' scores of zero predict 1. Scores are not probabilities. Invalid classifiers,
#' missing predictors, non-finite values, and score overflow cause an error.
#'
#' @return A numeric vector in observation order, named using the row names of
#'   \code{newdata}.
#' @importFrom stats predict
#' @export
#' @examples
#' train <- data.frame(x = c(-2, -1, 1, 2), y = c(-1, -1, 1, 1))
#' sol <- solution(0, 1, "MC")
#' sol$features <- "x"
#' sol$vectors <- list(1L, 4L)
#' sol$plane_coord <- 1
#' sol <- evaluate_solution(sol, train, "y")
#' predict(sol, newdata = data.frame(x = 0.5))
#' predict(sol, newdata = data.frame(x = c(-0.5, 0, 0.5)), type = "score")
predict.solution <- function(object, newdata, type = c("class", "score"), ...) {
  type <- match.arg(type)
  finite_numeric <- function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x))
  }
  if (!is.list(object) || !isTRUE(object$successful_evaluation)) {
    stop("A valid evaluated solution is required before prediction.", call. = FALSE)
  }
  features <- unlist(object$features, use.names = FALSE)
  weights <- object$plane_coord
  if (!is.character(features) || !length(features) || anyNA(features) ||
      anyDuplicated(features) || !finite_numeric(weights) ||
      length(weights) != length(features) ||
      !is.finite(sum(weights^2)) || sum(weights^2) <= 0 ||
      !.has_distinct_margin_planes(object)) {
    stop("The solution must have valid selected features and finite, distinct planes.", call. = FALSE)
  }
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.", call. = FALSE)
  newdata <- as.data.frame(newdata)
  if (!nrow(newdata) || anyNA(names(newdata)) || any(!nzchar(names(newdata))) ||
      anyDuplicated(names(newdata))) {
    stop("newdata must be nonempty with unique, nonempty column names.", call. = FALSE)
  }
  if (!all(features %in% names(newdata))) {
    stop("newdata must contain every selected feature.", call. = FALSE)
  }
  if (!all(vapply(newdata[features], finite_numeric, logical(1)))) {
    stop("Selected predictors must contain finite real numeric values.", call. = FALSE)
  }
  distances <- numeric(nrow(newdata))
  for (i in seq_along(features)) distances <- distances + weights[i] * newdata[[features[i]]]
  orientation <- if (object$plane_term_b[1] > object$plane_term_b[2]) 1 else -1
  scores <- orientation * (distances + object$plane_term_b[3])
  if (any(!is.finite(scores))) {
    stop("Decision scores must be finite; check feature scaling.", call. = FALSE)
  }
  names(scores) <- rownames(newdata)
  if (type == "score") scores else .prediction_classes(scores)
}

.prediction_classes <- function(scores) ifelse(scores >= 0, 1, -1)
