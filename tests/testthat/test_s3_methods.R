s3_test_solution <- function(direction = 1) {
  data <- data.frame(x = c(2, 4, 8, 10), y = c(-1, -1, 1, 1))
  sol <- solution(0, 1, "MC")
  sol$features <- "x"
  sol$vectors <- list(1L, 4L)
  sol$plane_coord <- direction
  evaluate_solution(sol, data, "y")
}

test_that("solution prediction dispatches and preserves ordinary list compatibility", {
  for (direction in c(1, -1)) {
    sol <- s3_test_solution(direction)
    data <- data.frame(x = c(2, 6, 10), row.names = c("a", "b", "c"))
    expect_identical(predict(sol, data), c(a = -1, b = 1, c = 1))
    expect_identical(predict(sol, data, type = "score"), c(a = -4, b = 0, c = 4))
    expect_identical(predict(sol, data), predict_solution(unclass(sol), data))
    before <- sol
    text <- capture.output(printed <- withVisible(print(sol)))
    expect_match(paste(text, collapse = "\n"), "Accuracy")
    expect_false(printed$visible)
    expect_identical(printed$value, sol)
    expect_identical(sol, before)
  }
  expect_output(print(solution(1, 2, "MC")), "Not evaluated")
  expect_error(predict(solution(1, 2, "MC"), data.frame(x = 1)), "evaluated")
})

test_that("summaries handle empty and populated optimizers without modifying them", {
  data <- data.frame(x = c(2, 4, 8, 10), y = c(-1, -1, 1, 1))
  obj <- svm_fs_meta(data, "x", "y", pop_size = 4, num_fea = 1,
                     n_iter = 1, clones = 1, objective = "MC")
  before <- obj
  empty <- summary(obj)
  expect_s3_class(empty, "summary.svm_fs_meta")
  expect_equal(empty$n_samples, 4)
  expect_equal(empty$population$size, 0)
  expect_equal(nrow(empty$population$solutions), 0)
  expect_output(print(obj), "Population: 0 solutions")
  expect_output(print(obj$population), "Population: 0 solutions")
  expect_output(print(empty), "SVM feature selection")
  expect_identical(obj, before)

  set.seed(42)
  capture.output(fit <- run(obj))
  expect_s3_class(fit, "svm_fs_meta")
  result <- summary(fit)
  expect_s3_class(result$population, "summary.population")
  expect_equal(result$population$size, 4)
  expect_equal(sum(result$population$front_sizes$n_solutions), 4)
  rows <- result$population$solutions
  expect_equal(rows$SOL, vapply(fit$population$solution_list, function(s) s$num, numeric(1)))
  expect_equal(rows$MCPOS, vapply(fit$population$solution_list, function(s) s$objective[3], numeric(1)))
  expect_true(all(rows$EVALUATED))
  expect_output(print(result, row.names = FALSE), "front")
  set.seed(42)
  capture.output(legacy <- svmfsmeta::run.svm_fs_meta(obj))
  expect_identical(fit, legacy)

  grDevices::pdf(tempfile(fileext = ".pdf"))
  tryCatch({
    expect_message(expect_null(plot(obj)), "No solutions")
    for (mode in c("MC", "DIST-EPS")) {
      fit$objective <- fit$population$objective <- mode
      first <- rows[rows$FRONT == 1, , drop = FALSE]
      expected_x <- if (mode == "MC") first$MCPOS else first$DIST
      expected_y <- if (mode == "MC") first$MCNEG else first$EPS
      for (graph in list(plot(fit), plot(fit$population))) {
        expect_s3_class(graph, "ggplot")
        coords <- ggplot2::ggplot_build(graph)$data[[1]]
        expect_equal(coords$x, expected_x)
        expect_equal(coords$y, expected_y)
      }
    }
  }, finally = grDevices::dev.off())
})

test_that("mutations reject stale predictions until the solution is reevaluated", {
  data <- data.frame(x = c(2, 4, 8, 10), y = c(-1, -1, 1, 1))
  sol <- s3_test_solution()
  pop <- population(3, 1, p_mut_ind = 0, p_mut_fea = 0,
                    p_mut_coord = 1, mut_coord = 1)
  set.seed(42)
  changed <- mutate_solution(sol, pop, data, "x", "y")
  expect_false(changed$successful_evaluation)
  expect_true(all(is.na(changed$objective)))
  expect_null(changed$accuracy)
  expect_null(changed$conf_matrix)
  expect_error(predict(changed, data.frame(x = 6)), "evaluated")
  expect_error(predict_solution(changed, data.frame(x = 6)), "evaluated")
  changed <- evaluate_solution(changed, data, "y")
  expect_equal(unname(predict(changed, data.frame(x = 6), type = "score")), 0)

  changed <- mutate_vectors(sol, data, "y")
  expect_false(changed$successful_evaluation)
  expect_error(predict(changed, data), "evaluated")
  pop$p_mut_coord <- 0
  expect_identical(mutate_solution(sol, pop, data, "x", "y"), sol)
  single <- data.frame(x = c(-1, 1), y = c(-1, 1))
  sol$vectors <- list(1L, 2L)
  expect_identical(mutate_vectors(sol, single, "y"), sol)
})
