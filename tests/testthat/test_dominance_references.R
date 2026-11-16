reference_population <- function(mode, ids = c(40, 0, 9), target = length(ids)) {
  pop <- population(target, 1, objective = mode)
  pop$solution_list <- lapply(seq_along(ids), function(i) {
    sol <- solution(ids[[i]], 1, mode)
    sol$objective <- c(length(ids) - i + 1, i, i, i)
    sol
  })
  fnds(pop)
}

test_that("dominance references use SOL identifiers rather than positions", {
  for (mode in c("DIST-EPS", "MC")) {
    pop <- reference_population(mode)
    expect_equal(vapply(pop$solution_list, function(s) s$front, numeric(1)), 1:3)
    expect_equal(vapply(pop$solution_list, function(s) s$sol_dom_by, numeric(1)), 0:2)
    expect_equal(pop$fronts[[1]]$solutions, 1)
    expect_equal(pop$fronts[[2]]$solutions, 2)
    expect_equal(pop$fronts[[3]]$solutions, 3)
    records <- lapply(pop$solution_list, to_dict)
    expect_equal(vapply(records, function(s) s$SOL_DOM_BY, numeric(1)), 0:2)
    expect_equal(records[[1]]$SOL, 40)
    expect_equal(unlist(records[[1]]$DOMINATES_TO), c(0, 9))
    expect_identical(records[[1]]$DOMINATED_BY, list())
    expect_equal(unlist(records[[2]]$DOMINATES_TO), 9)
    expect_equal(unlist(records[[2]]$DOMINATED_BY), 40)
    expect_identical(records[[3]]$DOMINATES_TO, list())
    expect_equal(unlist(records[[3]]$DOMINATED_BY), c(40, 0))

    repeated <- fnds(pop)
    expect_identical(repeated, pop)
    reordered <- pop
    reordered$solution_list <- reordered$solution_list[c(3, 1, 2)]
    reordered <- fnds(reordered)
    expect_equal(vapply(reordered$solution_list, function(s) s$front, numeric(1)), c(3, 1, 2))
    expect_equal(vapply(reordered$solution_list, function(s) s$sol_dom_by, numeric(1)), c(2, 0, 1))
    expect_equal(unlist(reordered$solution_list[[1]]$list_dominated_by), c(40, 0))
    expect_equal(unlist(reordered$solution_list[[2]]$dominates_list), c(9, 0))
    expect_equal(unlist(reordered$solution_list[[3]]$dominates_list), 9)
  }
})

test_that("survivor renumbering clears obsolete references before rebuilding", {
  for (mode in c("DIST-EPS", "MC")) {
    for (count in 2:4) {
      pop <- reference_population(mode, c(40, 0, 9, 77)[seq_len(count)], target = 3)
      original <- pop
      reduced <- reduce_population(pop)
      expect_identical(pop, original)
      expect_equal(vapply(reduced$solution_list, function(s) s$num, numeric(1)),
                   seq_len(min(count, 3)) - 1)
      for (sol in reduced$solution_list) {
        expect_identical(sol$dominates_list, list())
        expect_identical(sol$list_dominated_by, list())
        expect_equal(sol$front, -1)
      }
      rebuilt <- fnds(reduced)
      expect_equal(unlist(rebuilt$solution_list[[1]]$dominates_list),
                   seq_len(length(rebuilt$solution_list) - 1))
      expect_equal(unlist(rebuilt$solution_list[[2]]$list_dominated_by), 0)
      expect_equal(vapply(rebuilt$solution_list, function(s) s$front, numeric(1)),
                   seq_along(rebuilt$solution_list))
      expect_equal(vapply(rebuilt$solution_list, function(s) s$sol_dom_by, numeric(1)),
                   seq_along(rebuilt$solution_list) - 1)
    }
  }
})

test_that("optimized summaries and CSVs reference current population IDs", {
  data <- data.frame(x = seq(-2, 2, length.out = 20),
                     z = sin(seq_len(20)), y = rep(c(-1, 1), each = 10))
  for (mode in c("DIST-EPS", "MC")) {
    set.seed(72)
    result <- run.svm_fs_meta(svm_fs_meta(
      data, c("x", "z"), "y", pop_size = 4, num_fea = 1,
      n_iter = 2, objective = mode
    ))
    solutions <- result$population$solution_list
    for (sol in solutions) {
      expected <- vapply(Filter(function(other) dominate(sol, other), solutions),
                         function(other) other$num, numeric(1))
      expect_setequal(as.numeric(unlist(sol$dominates_list)), expected)
      expected <- vapply(Filter(function(other) dominate(other, sol), solutions),
                         function(other) other$num, numeric(1))
      expect_setequal(as.numeric(unlist(sol$list_dominated_by)), expected)
      expect_equal(sol$sol_dom_by, length(expected))
    }
    output <- save_results(result, tempdir(), "reference-test")
    exported <- utils::read.csv(output, colClasses = "character", na.strings = character())
    for (i in seq_len(nrow(exported))) {
      sol <- Filter(function(s) s$num == as.numeric(exported$SOL[i]), solutions)[[1]]
      expect_identical(exported$DOMINATES_TO[i], paste(unlist(sol$dominates_list), collapse = ","))
      expect_identical(exported$DOMINATED_BY[i], paste(unlist(sol$list_dominated_by), collapse = ","))
      expect_equal(as.numeric(exported$SOL_DOM_BY[i]), length(sol$list_dominated_by))
    }
  }
})
