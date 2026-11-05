generation_test_solution <- function(num, objective, key = 1,
                                     successful = TRUE) {
  sol <- solution(num, num_features = 1, obj = objective)
  sol$features <- "x"
  sol$plane_coord <- key / 10
  sol$vectors <- list(1L, 4L)
  sol$objective <- c(key, key + 1, key + 2, key + 3)
  sol$successful_evaluation <- successful
  sol
}

generation_test_population <- function(objective, clones = 0, probability = 0) {
  pop <- population(pop_size = 3, num_features = 1, clones = clones,
                    p_mutation = probability, objective = objective)
  pop$solution_list <- lapply(0:2, function(id) {
    generation_test_solution(id, objective, key = id + 1)
  })
  pop
}

generation_test_data <- function() {
  data.frame(x = c(-3, -2, -1, 1, 2, 3), y = rep(c(-1, 1), each = 3))
}

generation_test_mocks <- function(candidate, mutate = NULL, evaluate = NULL,
                                  .env = parent.frame()) {
  calls <- new.env(parent = emptyenv())
  calls$parents <- 0L
  calls$crossovers <- 0L
  calls$evaluations <- 0L
  calls$mutations <- 0L
  calls$random <- 0L

  testthat::local_mocked_bindings(
    tournament_select_parent = function(population) {
      calls$parents <- calls$parents + 1L
      population$solution_list[[1L + ((calls$parents - 1L) %% 2L)]]
    },
    crossover_solutions = function(population, parent, mother, num, data, output) {
      calls$crossovers <- calls$crossovers + 1L
      position <- 2L * calls$crossovers
      list(candidate(position - 1L, num), candidate(position, num + 1L))
    },
    evaluate_solution = function(solution, data, output) {
      calls$evaluations <- calls$evaluations + 1L
      # Keep regressions finite even if an earlier implementation retries forever.
      if (calls$evaluations > 1000L) {
        stop("Test guard: generation exceeded 1000 evaluations.")
      }
      if (is.null(evaluate)) solution else evaluate(solution)
    },
    mutate_solution = function(solution, population, data, inputs, output) {
      calls$mutations <- calls$mutations + 1L
      if (is.null(mutate)) stop("Unexpected forced mutation.")
      mutate(solution)
    },
    generate_random_solution = function(solution, data, inputs, output) {
      calls$random <- calls$random + 1L
      stop("Unexpected random immigrant.")
    },
    .package = "svmfsmeta", .env = .env
  )
  calls
}

test_that("allowed clones fill a generation without forcing mutation", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 1, probability = 0)
    calls <- generation_test_mocks(function(position, num) {
      generation_test_solution(num, objective, key = 1)
    })

    result <- new_population(pop, generation_test_data(), "x", "y")

    expect_length(result$solution_list, 6)
    expect_identical(result$solution_list[1:3], pop$solution_list)
    expect_true(all(vapply(result$solution_list[4:6], function(sol) {
      check_clones(pop, sol) == 1
    }, logical(1))))
    expect_equal(calls$evaluations, 3)
    expect_equal(calls$mutations, 0)
    expect_equal(calls$random, 0)
  }
})

test_that("one hundred consecutive forbidden clones stop generation", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 0, probability = 0)
    calls <- generation_test_mocks(function(position, num) {
      generation_test_solution(num, objective, key = 1)
    })

    expect_error(new_population(pop, generation_test_data(), "x", "y"), "100")
    expect_equal(calls$evaluations, 100)
    expect_equal(calls$crossovers, 50)
    expect_equal(calls$mutations, 0)
    expect_equal(calls$random, 0)
  }
})

test_that("one hundred consecutive invalid offspring stop generation", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 1, probability = 0)
    calls <- generation_test_mocks(function(position, num) {
      generation_test_solution(num, objective, key = 4, successful = FALSE)
    })

    expect_error(new_population(pop, generation_test_data(), "x", "y"), "100")
    expect_equal(calls$evaluations, 100)
    expect_equal(calls$crossovers, 50)
    expect_equal(calls$mutations, 0)
    expect_equal(calls$random, 0)
  }
})

test_that("rejected offspring are replaced by later crosses with consecutive ids", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 0, probability = 0)
    calls <- generation_test_mocks(function(position, num) {
      if (position == 3L) {
        return(generation_test_solution(num, objective, key = 5, successful = FALSE))
      }
      key <- if (position %% 2L == 0L) 3 + position / 2 else 1
      generation_test_solution(num, objective, key)
    })

    result <- new_population(pop, generation_test_data(), "x", "y")
    ids <- vapply(result$solution_list, function(sol) sol$num, numeric(1))

    expect_length(result$solution_list, 6)
    expect_equal(ids, 0:5)
    expect_equal(anyDuplicated(ids), 0)
    expect_true(all(vapply(result$solution_list, function(sol) {
      sol$successful_evaluation
    }, logical(1))))
    expect_equal(vapply(result$solution_list[4:6], function(sol) {
      check_clones(pop, sol)
    }, numeric(1)), rep(0, 3))
    expect_equal(calls$evaluations, 6)
    expect_equal(calls$crossovers, 3)
    expect_equal(calls$mutations, 0)
    expect_equal(calls$random, 0)
  }
})

test_that("accepting offspring resets the consecutive rejection count", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 0, probability = 0)
    calls <- generation_test_mocks(function(position, num) {
      # Three batches of 99 clones followed by one acceptable offspring.
      key <- if (position %% 100L == 0L) 3 + position / 100 else 1
      generation_test_solution(num, objective, key)
    })

    result <- new_population(pop, generation_test_data(), "x", "y")

    expect_length(result$solution_list, 6)
    expect_equal(vapply(result$solution_list, function(sol) sol$num, numeric(1)), 0:5)
    expect_equal(calls$evaluations, 300)
    expect_equal(calls$crossovers, 150)
    expect_equal(calls$mutations, 0)
    expect_equal(calls$random, 0)
  }
})

test_that("mutation probability one mutates every candidate before evaluation", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 1, probability = 1)
    calls <- generation_test_mocks(
      candidate = function(position, num) generation_test_solution(num, objective),
      mutate = function(sol) {
        sol$was_mutated <- TRUE
        sol$plane_coord <- sol$plane_coord + 0.01
        sol
      },
      evaluate = function(sol) {
        if (!isTRUE(sol$was_mutated)) stop("Mutation must precede evaluation.")
        sol
      }
    )

    result <- new_population(pop, generation_test_data(), "x", "y")

    expect_length(result$solution_list, 6)
    expect_true(all(vapply(result$solution_list[4:6], function(sol) {
      isTRUE(sol$was_mutated)
    }, logical(1))))
    expect_equal(calls$mutations, 3)
    expect_equal(calls$evaluations, 3)
    expect_equal(calls$random, 0)
  }
})

test_that("intermediate mutation probabilities follow their random draws", {
  for (objective in c("DIST-EPS", "MC")) {
    pop <- generation_test_population(objective, clones = 1, probability = 0.5)
    calls <- generation_test_mocks(
      candidate = function(position, num) generation_test_solution(num, objective),
      mutate = function(sol) {
        sol$was_mutated <- TRUE
        sol$plane_coord <- sol$plane_coord + 0.01
        sol
      }
    )
    set.seed(42)
    expected_mutations <- stats::runif(3) < 0.5
    set.seed(42)

    result <- new_population(pop, generation_test_data(), "x", "y")

    expect_identical(vapply(result$solution_list[4:6], function(sol) {
      isTRUE(sol$was_mutated)
    }, logical(1)), expected_mutations)
    expect_equal(calls$mutations, sum(expected_mutations))
    expect_equal(calls$evaluations, 3)
    expect_equal(calls$random, 0)
  }
})

test_that("real operators complete generations when all mutations are disabled", {
  x <- seq(-2, 2, length.out = 40)
  data <- data.frame(x1 = x, x2 = sin(2 * x) + x / 3,
                     y = rep(c(-1, 1), 20))
  expected_settings <- list(clones = 1, p_mutation = 0, p_mut_ind = 0,
                            p_mut_fea = 0, p_mut_coord = 0)
  original_evaluate <- get("evaluate_solution", envir = asNamespace("svmfsmeta"))
  evaluations <- 0L
  testthat::local_mocked_bindings(
    evaluate_solution = function(solution, data, output) {
      evaluations <<- evaluations + 1L
      if (evaluations > 500L) {
        stop("Test guard: real generation exceeded 500 evaluations.")
      }
      original_evaluate(solution, data, output)
    },
    .package = "svmfsmeta"
  )

  for (objective in c("DIST-EPS", "MC")) {
    evaluations <- 0L
    set.seed(572)
    args <- c(list(data = data, inputs = c("x1", "x2"), output = "y",
                   pop_size = 6, num_fea = 1, n_iter = 2,
                   objective = objective), expected_settings)
    result <- run.svm_fs_meta(do.call(svm_fs_meta, args))

    expect_length(result$population$solution_list, 6)
    expect_identical(result$population[names(expected_settings)], expected_settings)
    expect_true(all(vapply(result$population$solution_list, function(sol) {
      sol$successful_evaluation
    }, logical(1))))
    expect_lte(evaluations, 500)
  }
})
