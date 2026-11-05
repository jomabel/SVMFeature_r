test_that("F1 is zero when precision and recall are both zero", {
  y_true <- rep(c(-1, 1), 20)
  y_pred <- -y_true
  result <- calculate_metrics(solution(1, 1, "MC"), y_true, y_pred, y_pred)

  expect_equal(result$precision, 0)
  expect_equal(result$recall, 0)
  expect_equal(result$f1, 0)
})

test_that("F1 is one for perfect predictions", {
  y_true <- rep(c(-1, 1), 20)
  result <- calculate_metrics(solution(1, 1, "MC"), y_true, y_true, y_true)

  expect_equal(result$precision, 1)
  expect_equal(result$recall, 1)
  expect_equal(result$f1, 1)
})

test_that("F1 agrees with nonzero confusion matrix counts", {
  # TP = 30, FN = 10, FP = 20, TN = 20.
  y_true <- c(rep(1, 40), rep(-1, 40))
  y_pred <- c(rep(1, 30), rep(-1, 10), rep(1, 20), rep(-1, 20))
  result <- calculate_metrics(solution(1, 1, "MC"), y_true, y_pred, y_pred)

  expect_equal(result$precision, 0.6)
  expect_equal(result$recall, 0.75)
  expect_equal(result$f1, 2 / 3)
})

test_that("F1 remains zero when there are no predicted positives", {
  y_true <- rep(c(-1, 1), 20)
  y_pred <- rep(-1, length(y_true))
  result <- calculate_metrics(solution(1, 1, "MC"), y_true, y_pred, y_pred)

  expect_true(is.na(result$precision))
  expect_equal(result$recall, 0)
  expect_equal(result$f1, 0)
})

test_that("F1 remains zero when there are no actual positives", {
  y_true <- rep(-1, 40)
  y_pred <- rep(c(-1, 1), 20)
  result <- calculate_metrics(solution(1, 1, "MC"), y_true, y_pred, y_pred)

  expect_equal(result$precision, 0)
  expect_true(is.na(result$recall))
  expect_equal(result$f1, 0)
})
