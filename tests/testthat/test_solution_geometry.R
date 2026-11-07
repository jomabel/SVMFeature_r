test_that("coincident margins are rejected for both orientations and objectives", {
  dat <- data.frame(
    x1 = c(0, -1, 1, 0),
    x2 = c(-1, -2, 2, 1),
    y = c(-1, -1, 1, 1)
  )

  for (mode in c("DIST-EPS", "MC")) {
    for (weight in c(-1, 1)) {
      sol <- solution(1, 2, mode)
      sol$features <- c("x1", "x2")
      sol$plane_coord <- c(weight, 0)
      sol$vectors <- list(1L, 4L)
      # Crossovers may carry an earlier evaluation into another evaluation.
      sol$successful_evaluation <- TRUE
      sol$objective <- c(4, 0, 0, 0)
      sol$accuracy <- 1
      sol$conf_matrix <- diag(2)

      expect_message(result <- evaluate_solution(sol, dat, "y"), "distinct")
      expect_false(result$successful_evaluation)
      expect_true(all(is.na(result$objective)))
      expect_false(is_valid_solution(result))
      expect_null(result$accuracy)
      expect_null(result$conf_matrix)

      planes <- construct_planes(sol, dat)
      expect_message(direct <- calculate_epsilon_objective(planes, dat, "y"),
                     "distinct")
      expect_false(direct$successful_evaluation)
      expect_true(all(is.na(direct$objective)))
      expect_null(direct$conf_matrix)
    }
  }
})

test_that("tiny distinct margins remain valid when the normal is rescaled", {
  dat <- data.frame(x = c(-1e-10, 1e-10), y = c(-1, 1))

  for (mode in c("DIST-EPS", "MC")) {
    for (weight in c(-1e6, -1, -1e-6, 1e-6, 1, 1e6)) {
      sol <- solution(1, 1, mode)
      sol$features <- "x"
      sol$plane_coord <- weight
      sol$vectors <- list(1L, 2L)

      result <- evaluate_solution(sol, dat, "y")

      expect_true(result$successful_evaluation)
      expect_equal(result$objective[1] / 2e-10, 1)
      expect_equal(result$objective[2:4], c(0, 0, 0))
      expect_equal(result$accuracy, 1)
    }
  }
})

test_that("margins separated only by floating point noise are rejected", {
  dat <- data.frame(x = c(1, 1 + .Machine$double.eps), y = c(-1, 1))
  sol <- solution(1, 1, "MC")
  sol$features <- "x"
  sol$plane_coord <- 1
  sol$vectors <- list(1L, 2L)

  expect_message(result <- evaluate_solution(sol, dat, "y"), "distinct")
  expect_false(result$successful_evaluation)
  expect_false(is_valid_solution(result))
})

test_that("a nonzero normal remains valid when its coordinates sum to zero", {
  dat <- data.frame(x1 = c(-1, 1), x2 = c(1, -1), y = c(-1, 1))
  sol <- solution(1, 2, "DIST-EPS")
  sol$features <- c("x1", "x2")
  sol$plane_coord <- c(1, -1)
  sol$vectors <- list(1L, 2L)

  result <- evaluate_solution(sol, dat, "y")

  expect_true(result$successful_evaluation)
  expect_equal(result$objective, c(2 * sqrt(2), 0, 0, 0))
})

test_that("failed reevaluation clears objectives and cached metrics", {
  dat <- data.frame(x = c(-1, 1), y = c(-1, 1))
  sol <- solution(1, 1, "MC")
  sol$features <- "x"
  sol$plane_coord <- 1
  sol$vectors <- list(1L, 2L)
  evaluated <- evaluate_solution(sol, dat, "y")
  evaluated$plane_coord <- 0

  expect_message(result <- evaluate_solution(evaluated, dat, "y"), "normal")
  expect_false(result$successful_evaluation)
  expect_false(is_valid_solution(result))
  expect_true(all(is.na(result$objective)))
  expect_null(result$accuracy)
  expect_null(result$conf_matrix)
})

test_that("overflowing margin intercepts cannot produce successful evaluations", {
  dat <- data.frame(x1 = c(-1e308, 1e308),
                    x2 = c(-1e308, 1e308), y = c(-1, 1))
  sol <- solution(1, 2, "MC")
  sol$features <- c("x1", "x2")
  sol$plane_coord <- c(1, 1)
  sol$vectors <- list(1L, 2L)

  expect_message(result <- evaluate_solution(sol, dat, "y"), "finite")
  expect_false(result$successful_evaluation)
  expect_false(is_valid_solution(result))
})
