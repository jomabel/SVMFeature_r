test_that("equivalent plane orientations give the same correct predictions", {
  data <- data.frame(x = c(-1, 1), y = c(-1, 1))

  for (objective in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(num = 1, num_features = 1, obj = objective)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 2L)

      result <- evaluate_solution(sol, data, "y")

      expect_true(result$successful_evaluation)
      expect_equal(result$objective, c(2, 0, 0, 0))
      expect_equal(as.integer(result$conf_matrix), c(1L, 0L, 0L, 1L))
      expect_equal(result$accuracy, 1)
      expect_equal(result$f1, 1)
      expect_equal(result$auc, 1)
    }
  }
})

test_that("oriented continuous scores preserve AUC for offset planes", {
  data <- data.frame(x = c(2, 8, 4, 10), y = c(-1, -1, 1, 1))

  for (objective in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(num = 1, num_features = 1, obj = objective)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 4L)

      result <- evaluate_solution(sol, data, "y")

      expect_equal(result$objective, c(8, 12, 1, 1))
      expect_equal(as.integer(result$conf_matrix), rep(1L, 4))
      expect_equal(result$accuracy, 0.5)
      # Continuous scores give 0.75; using predicted labels would give 0.5.
      expect_equal(result$auc, 0.75)
    }
  }
})

test_that("points on the central plane are assigned to the positive class", {
  data <- data.frame(
    x = rep(c(-1, 0, 0, 1), 10),
    y = rep(c(-1, -1, 1, 1), 10)
  )

  for (weight in c(-1, 1)) {
    sol <- solution(num = 1, num_features = 1, obj = "MC")
    sol$features <- "x"
    sol$plane_coord <- weight
    sol$vectors <- list(1L, 4L)

    result <- evaluate_solution(sol, data, "y")

    expect_equal(as.integer(result$conf_matrix), c(10L, 0L, 10L, 20L))
    expect_equal(result$accuracy, 0.75)
    expect_equal(result$auc, 0.875)
  }
})
