population_edge_solution <- function(num, mode = "DIST-EPS", key = 1) {
  sol <- solution(num, num_features = 1, obj = mode)
  sol$features <- "x"
  sol$plane_coord <- key / 10
  sol$vectors <- list(1L, 4L)
  sol$objective <- c(10 - key, key, key, key)
  sol$successful_evaluation <- TRUE
  sol
}

population_edge_initial_mocks <- function(candidate, .env = parent.frame()) {
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  testthat::local_mocked_bindings(
    generate_random_solution = function(solution, data, inputs, output) {
      calls$count <- calls$count + 1L
      if (calls$count > 1000L) {
        stop("Test guard: initialization exceeded 1000 candidates.")
      }
      candidate(solution, calls$count)
    },
    evaluate_solution = function(solution, data, output) solution,
    .package = "svmfsmeta", .env = .env
  )
  calls
}

test_that("sorting accepts empty populations and clears obsolete fronts", {
  pop <- population(pop_size = 3, num_features = 1)
  pop$fronts <- list(list(solutions = 1L))

  result <- fnds(pop)

  expect_identical(result$solution_list, list())
  expect_identical(result$fronts, list())
})

test_that("a single candidate is assigned to the first front", {
  for (mode in c("DIST-EPS", "MC")) {
    pop <- population(pop_size = 3, num_features = 1, objective = mode)
    sol <- population_edge_solution(7, mode)
    sol$sol_dom_by <- 3
    sol$dominates_list <- list(2)
    sol$list_dominated_by <- list(3)
    sol$front <- 5
    pop$solution_list <- list(sol)

    result <- fnds(pop)

    expect_equal(result$solution_list[[1]]$front, 1)
    expect_equal(result$fronts[[1]]$solutions, 1)
    expect_equal(result$solution_list[[1]]$sol_dom_by, 0)
    expect_identical(result$solution_list[[1]]$dominates_list, list())
    expect_identical(result$solution_list[[1]]$list_dominated_by, list())
  }
})

test_that("reduction retains populations already at or below the target", {
  for (mode in c("DIST-EPS", "MC")) {
    for (count in 0:3) {
      pop <- population(pop_size = 3, num_features = 1, objective = mode,
                        clones = 1, p_mutation = 0.13, mut_coord = 1)
      pop$solution_list <- lapply(seq_len(count), function(i) {
        population_edge_solution(i + 10, mode, key = i)
      })
      pop <- fnds(pop)
      pop$df_solutions <- data.frame(stale = 1)

      result <- reduce_population(pop)

      expect_length(result$solution_list, count)
      expect_equal(vapply(result$solution_list, function(sol) sol$num, numeric(1)),
                   seq_len(count) - 1)
      expect_equal(lapply(result$solution_list, function(sol) sol$objective),
                   lapply(pop$solution_list, function(sol) sol$objective))
      expect_identical(result$fronts, list())
      expect_null(result$df_solutions)
      expect_equal(result$pop_size, 3)
      expect_equal(result$p_mutation, 0.13)
      expect_equal(result$clones, 1)
      expect_identical(result$objective, mode)
    }
  }
})

test_that("reduction reports exhausted fronts instead of indexing past them", {
  pop <- population(pop_size = 3, num_features = 1)
  pop$solution_list <- lapply(1:4, function(i) population_edge_solution(i, key = i))
  pop <- fnds(pop)
  pop$fronts <- pop$fronts[1]

  expect_error(reduce_population(pop), "fronts are incomplete.*fnds")

  pop$fronts <- list()
  expect_error(reduce_population(pop), "fronts are incomplete.*fnds")
})

test_that("coordinate mutation keeps either sign within the allowed interval", {
  dat <- data.frame(x1 = c(-2, -1, 1, 2), x2 = c(2, 1, -1, -2),
                    y = c(-1, -1, 1, 1))
  sol <- solution(1, num_features = 2, obj = "DIST-EPS")
  sol$features <- c("x1", "x2")
  sol$plane_coord <- c(-0.99, 0.99)
  sol$vectors <- list(1L, 4L)

  for (strategy in c(0, 1)) {
    pop <- population(pop_size = 3, num_features = 2, p_mut_ind = 0,
                      p_mut_fea = 0, p_mut_coord = 1, mut_coord = strategy)
    set.seed(42)
    coordinates <- replicate(50, {
      mutate_solution(sol, pop, dat, sol$features, "y")$plane_coord
    })
    expect_true(all(is.finite(coordinates)))
    expect_true(all(coordinates >= -1 & coordinates <= 1))
  }
})

test_that("valid offspring are preferred to candidates with invalid objectives", {
  set.seed(1)
  for (mode in c("DIST-EPS", "MC")) {
    valid <- population_edge_solution(1, mode)
    valid$successful_evaluation <- FALSE # Manually supplied objectives remain supported.
    invalid <- population_edge_solution(2, mode)
    invalid$objective <- rep(NA_real_, 4)

    first <- select_best_offspring(valid, invalid)
    second <- select_best_offspring(invalid, valid)

    expect_identical(first$winner, valid)
    expect_identical(first$loser, invalid)
    expect_identical(second$winner, valid)
    expect_identical(second$loser, invalid)
  }
})

test_that("one hundred consecutive invalid candidates stop initialization", {
  calls <- population_edge_initial_mocks(function(sol, count) {
    sol$successful_evaluation <- FALSE
    sol
  })
  pop <- population(pop_size = 3, num_features = 1, clones = 1)

  expect_error(generate_initial_population(pop, data.frame(), "x", "y"),
               "initial population after 100 consecutive rejections")
  expect_equal(calls$count, 100)
})

test_that("forbidden duplicates also stop initialization after one hundred rejections", {
  calls <- population_edge_initial_mocks(function(sol, count) {
    population_edge_solution(sol$num)
  })
  pop <- population(pop_size = 3, num_features = 1, clones = 0)

  expect_error(generate_initial_population(pop, data.frame(), "x", "y"),
               "initial population after 100 consecutive rejections")
  expect_equal(calls$count, 101) # The first candidate is accepted.
})

test_that("an accepted initial candidate resets the rejection counter", {
  calls <- population_edge_initial_mocks(function(sol, count) {
    sol <- population_edge_solution(sol$num)
    sol$successful_evaluation <- count %% 100L == 0L
    sol
  })
  pop <- population(pop_size = 2, num_features = 1, clones = 1)

  result <- generate_initial_population(pop, data.frame(), "x", "y")

  expect_length(result$solution_list, 2)
  expect_equal(calls$count, 200)
  expect_equal(vapply(result$solution_list, function(sol) sol$num, numeric(1)), 0:1)
})

test_that("initialization appends candidates without replacing a partial population", {
  calls <- population_edge_initial_mocks(function(sol, count) {
    population_edge_solution(sol$num, key = count + 1)
  })
  pop <- population(pop_size = 3, num_features = 1)
  existing <- population_edge_solution(0)
  pop$solution_list <- list(existing)

  result <- generate_initial_population(pop, data.frame(), "x", "y")

  expect_identical(result$solution_list[[1]], existing)
  expect_length(result$solution_list, 3)
  expect_equal(calls$count, 2)
  expect_equal(vapply(result$solution_list, function(sol) sol$num, numeric(1)), 0:2)
})

test_that("constant predictors stop initialization instead of retrying forever", {
  dat <- data.frame(x = rep(1, 4), y = c(-1, -1, 1, 1))
  pop <- population(pop_size = 3, num_features = 1, clones = 1)

  expect_error(suppressMessages(generate_initial_population(pop, dat, "x", "y")),
               "initial population after 100 consecutive rejections")
})
