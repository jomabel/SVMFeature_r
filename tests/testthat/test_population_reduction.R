reduction_test_settings <- function(mode, pop_size = 2) {
  list(
    pop_size = pop_size,
    num_features = 1,
    num_obj = 2,
    clones = 1,
    p_mutation = 0.17,
    p_mut_ind = 0.23,
    p_mut_fea = 0.31,
    p_mut_coord = 0.41,
    mut_coord = 1,
    objective = mode
  )
}

reduction_test_population <- function(mode, objectives) {
  pop <- do.call(population, reduction_test_settings(mode))
  pop$solution_list <- lapply(seq_along(objectives), function(i) {
    sol <- solution(num = i, num_features = 1, obj = mode)
    sol$objective <- objectives[[i]]
    sol
  })
  pop <- fnds(pop)
  pop$df_solutions <- data.frame(SOL = seq_along(objectives))
  pop
}

test_that("reduction preserves settings when selecting a complete front", {
  cases <- list(
    "DIST-EPS" = list(
      c(4, 1, 0, 0), c(5, 2, 0, 0),
      c(3, 3, 0, 0), c(2, 2, 0, 0)
    ),
    "MC" = list(
      c(0, 0, 1, 2), c(0, 0, 2, 1),
      c(0, 0, 3, 4), c(0, 0, 4, 3)
    )
  )

  for (mode in names(cases)) {
    pop <- reduction_test_population(mode, cases[[mode]])
    expected_settings <- reduction_test_settings(mode)

    reduced <- reduce_population(pop)

    expect_s3_class(reduced, "population")
    expect_identical(reduced[names(expected_settings)], expected_settings)
    expect_length(reduced$solution_list, 2)
    expect_equal(vapply(reduced$solution_list, function(sol) sol$num, numeric(1)),
                 c(0, 1))
    expect_identical(reduced$fronts, list())
    expect_null(reduced$df_solutions)
  }
})

test_that("reduction preserves settings when truncating a front by crowding", {
  cases <- list(
    "DIST-EPS" = list(
      c(0, 0, 1, 1), c(2, 3, 0, 3),
      c(5, 7, 3, 0), c(10, 10, 2, 2)
    ),
    "MC" = list(
      c(1, 1, 0, 10), c(0, 3, 2, 7),
      c(3, 0, 5, 3), c(2, 2, 10, 0)
    )
  )

  for (mode in names(cases)) {
    pop <- reduction_test_population(mode, cases[[mode]])
    expected_settings <- reduction_test_settings(mode)

    reduced <- reduce_population(pop)

    expect_identical(reduced[names(expected_settings)], expected_settings)
    expect_length(reduced$solution_list, 2)
    first_objective <- if (mode == "MC") 3L else 1L
    expect_equal(sort(vapply(reduced$solution_list, function(sol) {
      sol$objective[[first_objective]]
    }, numeric(1))), c(0, 10))
    expect_equal(vapply(reduced$solution_list, function(sol) sol$num, numeric(1)),
                 c(0, 1))
    expect_identical(reduced$fronts, list())
    expect_null(reduced$df_solutions)
  }
})

test_that("nondefault settings survive multiple optimization generations", {
  x <- seq(-2, 2, length.out = 40)
  data <- data.frame(
    x1 = x,
    x2 = sin(2 * x) + x / 3,
    y = rep(c(-1, 1), 20)
  )

  for (mode in c("DIST-EPS", "MC")) {
    set.seed(572)
    expected_settings <- reduction_test_settings(mode, pop_size = 6)
    args <- c(
      list(data = data, inputs = c("x1", "x2"), output = "y",
           num_fea = expected_settings$num_features, n_iter = 2, mode = "iters"),
      expected_settings[names(expected_settings) != "num_features"]
    )
    obj <- do.call(svm_fs_meta, args)

    result <- run.svm_fs_meta(obj)

    expect_identical(result$population[names(expected_settings)], expected_settings)
    expect_length(result$population$solution_list, 6)
    expect_true(all(vapply(result$population$solution_list, function(sol) {
      sol$successful_evaluation
    }, logical(1))))
  }
})
