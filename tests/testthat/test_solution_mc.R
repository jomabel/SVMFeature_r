test_that("correct predictions inside the margins do not count as MC errors", {
  data <- data.frame(x = c(-2, -1, 1, 2), y = c(-1, -1, 1, 1))

  for (objective in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(num = 1, num_features = 1, obj = objective)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 4L)

      result <- evaluate_solution(sol, data, "y")

      # The central plane is x = 0. Both interior points are correct,
      # although each contributes one unit to the margin penalty.
      expect_true(result$successful_evaluation)
      expect_equal(result$objective, c(4, 2, 0, 0))
      expect_equal(as.integer(result$conf_matrix), c(2L, 0L, 0L, 2L))
      expect_equal(result$accuracy, 1)
    }
  }
})

test_that("MC stores false negatives before false positives for either orientation", {
  data <- data.frame(
    x = c(-2, -1, 1, -1, -0.5, 1, 2),
    y = c(-1, -1, -1, 1, 1, 1, 1)
  )

  for (objective in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(num = 1, num_features = 1, obj = objective)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 7L)

      result <- evaluate_solution(sol, data, "y")

      # With the central plane at zero, the positive points at -1 and
      # -0.5 are FN; the negative point at 1 is FP. The margin penalty
      # is 1 + 3 + 3 + 2.5 + 1, independently of those three errors.
      expect_true(result$successful_evaluation)
      expect_equal(result$objective, c(4, 10.5, 2, 1))
      expect_equal(as.integer(result$conf_matrix), c(2L, 2L, 1L, 2L))
      expect_equal(result$accuracy, 4 / 7)

      exported <- to_dict(result)
      expect_equal(exported$MCPOS, 2)
      expect_equal(exported$MCNEG, 1)
    }
  }
})

test_that("MC uses the positive prediction for points on the central plane", {
  data <- data.frame(x = c(-1, 0, 0, 1), y = c(-1, -1, 1, 1))

  for (objective in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(num = 1, num_features = 1, obj = objective)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 4L)

      result <- evaluate_solution(sol, data, "y")

      # Both zeros are predicted as +1, so only the negative zero is
      # misclassified. Each zero still contributes one unit of epsilon.
      expect_true(result$successful_evaluation)
      expect_equal(result$objective, c(2, 2, 0, 1))
      expect_equal(as.integer(result$conf_matrix), c(1L, 0L, 1L, 2L))
      expect_equal(result$accuracy, 0.75)
    }
  }
})
