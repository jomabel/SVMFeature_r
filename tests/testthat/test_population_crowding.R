crowding_test_population <- function(mode, objectives, pop_size = length(objectives),
                                     fronts = rep(1L, length(objectives))) {
  pop <- population(pop_size = pop_size, num_features = 1, objective = mode)
  pop$solution_list <- lapply(seq_along(objectives), function(i) {
    sol <- solution(num = i, num_features = 1, obj = mode)
    sol$objective <- objectives[[i]]
    sol$front <- fronts[[i]]
    sol$crowding_distance <- 42
    sol
  })
  pop
}

test_that("crowding uses the active objectives and normalized distances", {
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
    pop <- crowding_test_population(mode, cases[[mode]])
    ranked <- crowding_distance(pop, num_front = 1)
    ids <- vapply(ranked, function(sol) sol$num, integer(1))
    distances <- vapply(ranked, function(sol) sol$crowding_distance, numeric(1))

    expect_equal(distances[order(ids)], c(Inf, 1.2, 1.5, Inf))
    expect_setequal(ids[1:2], c(1L, 4L))

    # Repeating the calculation must reset, rather than accumulate, distances.
    pop$solution_list <- ranked
    repeated <- crowding_distance(pop, num_front = 1)
    repeated_ids <- vapply(repeated, function(sol) sol$num, integer(1))
    repeated_distances <- vapply(repeated, function(sol) sol$crowding_distance,
                                numeric(1))
    expect_equal(repeated_distances[order(repeated_ids)], c(Inf, 1.2, 1.5, Inf))
  }
})

test_that("crowding only includes solutions from the requested front", {
  for (mode in c("DIST-EPS", "MC")) {
    pop <- crowding_test_population(
      mode,
      list(c(0, 0, 0, 10), c(10, 10, 10, 0), c(100, 100, 100, 100)),
      fronts = c(1L, 1L, 2L)
    )

    ranked <- crowding_distance(pop, num_front = 1)

    expect_setequal(vapply(ranked, function(sol) sol$num, integer(1)), c(1L, 2L))
    expect_identical(crowding_distance(pop, num_front = 99), list())
  }
})

test_that("fronts with one or two solutions have infinite crowding distances", {
  for (mode in c("DIST-EPS", "MC")) {
    for (size in 1:2) {
      objectives <- list(c(0, 0, 0, 10), c(10, 10, 10, 0))[seq_len(size)]
      pop <- crowding_test_population(mode, objectives, pop_size = 4)

      ranked <- crowding_distance(pop, num_front = 1)

      expect_length(ranked, size)
      expect_true(all(vapply(ranked, function(sol) {
        is.infinite(sol$crowding_distance)
      }, logical(1))))
    }
  }
})

test_that("constant objectives do not produce undefined crowding distances", {
  for (mode in c("DIST-EPS", "MC")) {
    pop <- crowding_test_population(mode, rep(list(c(2, 2, 3, 3)), 4))
    ranked <- crowding_distance(pop, num_front = 1)
    ids <- vapply(ranked, function(sol) sol$num, integer(1))
    distances <- vapply(ranked, function(sol) sol$crowding_distance, numeric(1))

    expect_equal(distances[order(ids)], c(Inf, 0, 0, Inf))
  }
})

test_that("population reduction retains the extremes of the MC front", {
  pop <- crowding_test_population(
    "MC",
    list(c(1, 1, 0, 10), c(0, 3, 2, 7),
         c(3, 0, 5, 3), c(2, 2, 10, 0)),
    pop_size = 2
  )
  pop <- fnds(pop)

  reduced <- reduce_population(pop)

  expect_length(reduced$solution_list, 2)
  expect_setequal(vapply(reduced$solution_list, function(sol) {
    paste(sol$objective[3:4], collapse = "/")
  }, character(1)), c("0/10", "10/0"))
  expect_true(all(vapply(reduced$solution_list, function(sol) {
    is.infinite(sol$crowding_distance)
  }, logical(1))))
})
