prediction_model <- function(direction = 1) {
  data <- data.frame(x = c(2, 4, 8, 10), y = c(-1, -1, 1, 1))
  sol <- solution(0, 1, "MC")
  sol$features <- "x"
  sol$plane_coord <- direction
  sol$vectors <- list(1L, 4L)
  evaluate_solution(sol, data, "y")
}

test_that("solutions dispatch predict for one or many unlabeled observations", {
  for (direction in c(1, -1)) {
    sol <- prediction_model(direction)
    before <- sol
    expect_s3_class(sol, "solution")
    data <- data.frame(unused = letters[1:3], x = c(2, 6, 10), row.names = c("c", "a", "b"))
    original <- data
    expect_equal(predict(sol, data), c(c = -1, a = 1, b = 1))
    expect_equal(predict(sol, data, type = "score"), c(c = -4, a = 0, b = 4))
    expect_equal(predict(sol, data[2, , drop = FALSE]), c(a = 1))
    expect_identical(sol, before)
    expect_identical(data, original)
    expect_error(predict(sol, data, type = "probability"), "arg")
  }
})

test_that("prediction validates predictors and evaluated planes", {
  sol <- prediction_model()
  data <- data.frame(x = 1)
  for (bad in list(data[FALSE, , drop = FALSE], data.frame(z = 1),
                   data.frame(x = NA_real_), data.frame(x = Inf),
                   data.frame(x = TRUE), data.frame(x = "1"), data.frame(x = 1i))) {
    expect_error(predict(sol, bad))
  }
  expect_error(predict(sol, 1), "data frame")
  expect_error(predict(solution(0, 1, "MC"), data), "evaluated")
  sol$plane_term_b <- c(1, 1, 1)
  expect_error(predict(sol, data), "distinct")
  sol <- prediction_model()
  sol$plane_coord <- 2
  expect_error(predict(sol, data.frame(x = .Machine$double.xmax)), "scores")
})

test_that("evaluation uses S3 prediction and supports old evaluated lists", {
  sol <- prediction_model()
  data <- data.frame(x = c(2, 6, 10), y = c(-1, 1, 1))
  # A method override proves evaluation calls prediction, rather than
  # independently recalculating the original stored planes.
  local_mocked_bindings(
    predict.solution = function(object, newdata, type = c("class", "score"), ...) {
      expect_identical(type, "score")
      stats::setNames(c(2, -1, -3), rownames(newdata))
    }, .package = "svmfsmeta"
  )
  result <- evaluate_test(sol, data, "y")
  expect_equal(unname(result$decision_scores), c(2, -1, -3))
  expect_equal(unname(result$predictions), c(1, -1, -1))
  expect_equal(result$error_rate, 1)
})

test_that("legacy evaluation uses a local class without changing saved lists", {
  legacy <- unclass(prediction_model())
  before <- legacy
  data <- data.frame(x = c(2, 6, 10), y = c(-1, 1, 1))
  result <- evaluate_test(legacy, data, "y")
  expect_equal(result$error_rate, 0)
  expect_identical(legacy, before)
})
