dominance_solution <- function(num, mode, objectives) {
  sol <- solution(num = num, num_features = 1, obj = mode)
  sol$objective <- objectives
  sol
}

test_that("equal active objectives do not imply dominance", {
  cases <- list(
    "DIST-EPS" = list(c(2, 3, 90, 80), c(2, 3, 1, 2)),
    "MC" = list(c(90, 80, 2, 3), c(1, 2, 2, 3))
  )

  for (mode in names(cases)) {
    a <- dominance_solution(1, mode, cases[[mode]][[1]])
    b <- dominance_solution(2, mode, cases[[mode]][[2]])

    expect_false(dominate(a, a))
    expect_false(dominate(a, b))
    expect_false(dominate(b, a))
    expect_equal(dominate2(a, b), 0)
    expect_equal(dominate2(b, a), 0)
  }
})

test_that("dominance requires improvement without worsening another objective", {
  cases <- list(
    "DIST-EPS" = list(
      worse = c(2, 3, 0, 0),
      better = list(c(3, 3, 9, 9), c(2, 2, 9, 9), c(3, 2, 9, 9))
    ),
    "MC" = list(
      worse = c(9, 9, 2, 3),
      better = list(c(0, 0, 1, 3), c(0, 0, 2, 2), c(0, 0, 1, 2))
    )
  )

  for (mode in names(cases)) {
    worse <- dominance_solution(1, mode, cases[[mode]]$worse)
    for (objectives in cases[[mode]]$better) {
      better <- dominance_solution(2, mode, objectives)

      expect_true(dominate(better, worse))
      expect_false(dominate(worse, better))
      expect_equal(dominate2(better, worse), 1)
      expect_equal(dominate2(worse, better), 2)
    }
  }
})

test_that("solutions with an objective tradeoff do not dominate each other", {
  cases <- list(
    "DIST-EPS" = list(c(1, 1, 0, 0), c(2, 2, 0, 0)),
    "MC" = list(c(0, 0, 1, 2), c(0, 0, 2, 1))
  )

  for (mode in names(cases)) {
    a <- dominance_solution(1, mode, cases[[mode]][[1]])
    b <- dominance_solution(2, mode, cases[[mode]][[2]])

    expect_false(dominate(a, b))
    expect_false(dominate(b, a))
    expect_equal(dominate2(a, b), 0)
    expect_equal(dominate2(b, a), 0)

    pop <- population(pop_size = 2, num_features = 1, objective = mode)
    pop$solution_list <- list(a, b)
    ranked <- fnds(pop)
    expect_equal(vapply(ranked$solution_list, function(sol) sol$front, numeric(1)),
                 c(1, 1))
  }
})

test_that("equal objectives share fronts independently of population order", {
  cases <- list(
    "DIST-EPS" = list(
      c(4, 1, 90, 90), c(4, 1, 0, 0),
      c(3, 2, 90, 90), c(3, 2, 0, 0), c(5, 3, 1, 1)
    ),
    "MC" = list(
      c(90, 90, 1, 1), c(0, 0, 1, 1),
      c(90, 90, 2, 2), c(0, 0, 2, 2), c(1, 1, 0, 3)
    )
  )
  expected_fronts <- c(1, 1, 2, 2, 1)

  for (mode in names(cases)) {
    solutions <- lapply(seq_along(cases[[mode]]), function(i) {
      dominance_solution(i, mode, cases[[mode]][[i]])
    })
    for (indices in list(1:5, c(5, 3, 1, 4, 2))) {
      pop <- population(pop_size = 5, num_features = 1, objective = mode)
      pop$solution_list <- solutions[indices]
      ranked <- fnds(pop)

      expect_equal(
        vapply(ranked$solution_list, function(sol) sol$front, numeric(1)),
        expected_fronts[indices]
      )
    }
  }
})
