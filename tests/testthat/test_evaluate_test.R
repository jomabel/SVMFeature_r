test_classifier <- function(direction = 1, shift = 0) {
  train <- data.frame(x = c(-2, -1, 1, 2) + shift, y = c(-1, -1, 1, 1))
  sol <- solution(0, 1, "MC")
  sol$features <- "x"
  sol$vectors <- list(1L, 4L)
  sol$plane_coord <- direction
  evaluate_solution(sol, train, "y")
}

test_that("held-out errors use the fixed plane, preserve rows and do not mutate", {
  sol <- test_classifier()
  original <- sol
  test <- data.frame(unused = letters[1:5], y = c(-1, 1, -1, 1, 1),
                     x = c(-3, -0.5, 0, 0.5, 3), row.names = c("e", "a", "c", "b", "d"))
  original_test <- test
  metrics <- evaluate_test(sol, test, "y")
  expect_equal(metrics$error_rate, 0.4)
  expect_equal(metrics$n_errors, 2)
  expect_equal(metrics$n_samples, 5)
  expect_equal(metrics$fn, 1)
  expect_equal(metrics$fp, 1)
  expect_equal(metrics$accuracy, 0.6)
  expect_equal(unname(metrics$predictions), c(-1, -1, 1, 1, 1))
  expect_equal(unname(metrics$decision_scores), test$x)
  expect_identical(metrics$predictions, predict(sol, test["x"], type = "class"))
  expect_identical(metrics$decision_scores, predict(sol, test["x"], type = "score"))
  expect_identical(names(metrics$predictions), rownames(test))
  expect_identical(names(metrics$decision_scores), rownames(test))
  expect_equal(as.vector(metrics$confusion_matrix), c(1, 1, 1, 2))
  expect_identical(sol, original)
  expect_identical(test, original_test)
})

test_that("held-out scores preserve orientation, nonzero intercepts and training parity", {
  train <- data.frame(x = c(2, 4, 8, 10), y = c(-1, -1, 1, 1))
  for (direction in c(1, -1)) {
    sol <- test_classifier(direction, shift = 6)
    sol <- evaluate_solution(sol, train, "y")
    metrics <- evaluate_test(sol, train, "y")
    expect_equal(unname(metrics$decision_scores), c(-4, -2, 2, 4))
    for (metric in c("accuracy", "precision", "recall", "f1", "kappa", "auc")) {
      expect_equal(metrics[[metric]], sol[[metric]])
    }
    expect_equal(metrics$confusion_matrix, sol$conf_matrix)
  }
})

test_that("one-class test sets return undefined metrics without failure", {
  metrics <- evaluate_test(test_classifier(), data.frame(x = c(-3, -1), y = c(-1, -1)), "y")
  expect_equal(metrics$error_rate, 0)
  expect_true(is.na(metrics$auc))
  expect_true(is.na(metrics$precision))
  expect_true(is.na(metrics$recall))
  expect_true(is.na(metrics$kappa))
  expect_equal(metrics$f1, 0)
})

test_that("invalid test data and invalid classifiers fail clearly", {
  sol <- test_classifier()
  test <- data.frame(x = c(-1, 1), y = c(-1, 1))
  bad <- list(test[FALSE, ], test["y"], transform(test, x = NA_real_),
              transform(test, x = Inf), transform(test, x = "1"),
              transform(test, x = 1i), transform(test, x = TRUE),
              transform(test, y = 0), transform(test, y = NA_real_),
              transform(test, y = TRUE), cbind(test, test["x"]))
  for (data in bad) expect_error(evaluate_test(sol, data, "y"))
  expect_error(evaluate_test(sol, as.matrix(test), "y"), "data frame")
  expect_error(evaluate_test(sol, test, "missing"), "output")
  expect_error(evaluate_test(sol, test, "x"), "output")
  expect_error(evaluate_test(solution(0, 1, "MC"), test, "y"), "evaluated")
  sol$plane_coord <- 0
  expect_error(evaluate_test(sol, test, "y"), "valid")
  sol <- test_classifier()
  sol$plane_coord <- 2
  expect_error(evaluate_test(sol, data.frame(x = .Machine$double.xmax, y = 1), "y"), "scores")
})

test_that("selected feature order is honored independently of test column order", {
  data <- data.frame(x = c(-2, -1, 1, 2), z = c(2, 1, -1, -2), y = c(-1, -1, 1, 1))
  sol <- test_classifier()
  sol$features <- c("z", "x")
  sol$plane_coord <- c(-1, 2)
  sol <- evaluate_solution(sol, data, "y")
  metrics <- evaluate_test(sol, data[c("y", "x", "z")], "y")
  expect_equal(unname(metrics$decision_scores), c(-6, -3, 3, 6))
  expect_equal(metrics$error_rate, 0)
})

test_that("the installed Colon example evaluates each fixed Pareto solution", {
  skip_if_not_installed("plsgenomics")
  example <- system.file("examples", "colon_train_test.R", package = "svmfsmeta")
  expect_true(nzchar(example))
  env <- new.env(parent = globalenv())
  capture.output(sys.source(example, envir = env))
  expect_equal(nrow(env$train), 41)
  expect_equal(nrow(env$test), 21)
  expect_length(intersect(rownames(env$train), rownames(env$test)), 0)
  expect_equal(as.integer(table(env$train$y)), c(26, 15))
  expect_equal(as.integer(table(env$test$y)), c(14, 7))
  expect_true(nrow(env$results) > 0)
  for (sol in env$front) {
    before <- sol
    metrics <- evaluate_test(sol, env$test, "y")
    expected <- table(factor(env$test$y, levels = c(-1, 1)),
                      factor(metrics$predictions, levels = c(-1, 1)))
    expect_equal(metrics$confusion_matrix, expected)
    expect_equal(metrics$error_rate, mean(metrics$predictions != env$test$y))
    expect_identical(metrics$predictions, predict(sol, env$test[sol$features]))
    expect_identical(metrics$decision_scores,
                     predict(sol, env$test[sol$features], type = "score"))
    expect_identical(sol, before)
    row <- env$results[env$results$SOL == sol$num, ]
    expect_equal(row$TEST_FN, metrics$fn)
    expect_equal(row$TEST_FP, metrics$fp)
    expect_equal(row$TRAIN_ERROR_RATE, sum(sol$objective[3:4]) / 41)
  }
  feature <- env$inputs[[1]]
  original_train <- env$data[-env$test_rows, feature]
  expect_equal(env$training_min[[feature]], min(original_train))
  expect_equal(env$test[[feature]],
               (env$data[env$test_rows, feature] - min(original_train)) /
                 (max(original_train) - min(original_train)))
})
