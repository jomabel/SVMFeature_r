test_that("class-vector sampling handles singleton and absent classes", {
  dat <- data.frame(y = c(-1, -1, -1, 1))
  sol <- solution(1, 1, "MC")

  set.seed(42)
  positive_rows <- replicate(50, get_class_vector(sol, 1, dat, "y"))
  expect_identical(positive_rows, rep(4L, 50))
  negative_rows <- replicate(50, get_class_vector(sol, -1, dat, "y"))
  expect_true(all(dat$y[negative_rows] == -1))
  expect_identical(get_class_vector(sol, 2, dat, "y"), NA_integer_)
})

test_that("vector mutation terminates when neither class has an alternative", {
  dat <- data.frame(y = c(-1, 1))
  sol <- solution(1, 1, "MC")
  sol$vectors <- list(1L, 2L)
  original_sampler <- get_class_vector
  calls <- 0L

  # Keep this regression finite with the old rejection-sampling loop.
  local_mocked_bindings(
    get_class_vector = function(...) {
      calls <<- calls + 1L
      if (calls > 10L) stop("Test guard: singleton mutation kept retrying.")
      original_sampler(...)
    },
    .package = "svmfsmeta"
  )

  result <- mutate_vectors(sol, dat, "y")
  expect_identical(result$vectors, sol$vectors)
})

test_that("vector mutation changes rows while preserving their classes", {
  dat <- data.frame(y = c(-1, -1, -1, 1, 1, 1))
  sol <- solution(1, 1, "MC")
  sol$vectors <- list(1L, 6L)

  set.seed(42)
  for (i in seq_len(30)) {
    result <- mutate_vectors(sol, dat, "y")
    rows <- unlist(result$vectors)
    expect_true(all(rows != unlist(sol$vectors)))
    expect_equal(dat$y[rows], c(-1, 1))
  }

  # Both remaining candidates are singletons greater than one where possible.
  pair <- data.frame(y = c(-1, -1, 1, 1))
  sol$vectors <- list(1L, 3L)
  expect_identical(mutate_vectors(sol, pair, "y")$vectors, list(2L, 4L))
})

test_that("vector mutation selects valid rows for missing or mismatched vectors", {
  dat <- data.frame(y = c(-1, 1))
  sol <- solution(1, 1, "MC")
  expect_identical(mutate_vectors(sol, dat, "y")$vectors, list(1L, 2L))
  sol$vectors <- list(2L, 1L)
  expect_identical(mutate_vectors(sol, dat, "y")$vectors, list(1L, 2L))
})

test_that("solution validation rejects non-finite or malformed active objectives", {
  invalid_values <- list(
    Inf, -Inf, NA_real_, NaN, 1 + 1i, "1", TRUE, NULL,
    list(1, 2, 3, 4), matrix(1:4, nrow = 2)
  )

  for (mode in c("DIST-EPS", "MC")) {
    active_index <- if (mode == "DIST-EPS") 1L else 3L
    valid <- solution(1, 1, mode)
    valid$objective <- c(2, 3, 4, 5)
    expect_true(is_valid_solution(valid))

    for (value in invalid_values) {
      invalid <- valid
      if (is.logical(value)) {
        invalid$objective <- rep(value, 4L)
      } else if (length(value) == 1L && is.atomic(value) && is.null(dim(value))) {
        invalid$objective[active_index] <- value
      } else {
        invalid$objective <- value
      }
      expect_false(is_valid_solution(invalid))
      expect_false(dominate(invalid, valid))
    }

    short <- valid
    short$objective <- if (mode == "DIST-EPS") 1 else c(1, 2, 3)
    expect_false(is_valid_solution(short))
  }
})

test_that("solution validation checks objective names and only active entries", {
  sol <- solution(1, 1, "MC")
  sol$objective <- c(Inf, NA_real_, 2, 3)
  expect_true(is_valid_solution(sol))
  sol$obj_fn <- "DIST-EPS"
  sol$objective <- c(2, 3, Inf, NA_real_)
  expect_true(is_valid_solution(sol))

  for (mode in list("unknown", NA_character_, character(), c("MC", "DIST-EPS"),
                    matrix("MC", 1), 1, NULL)) {
    sol$obj_fn <- mode
    expect_false(is_valid_solution(sol))
  }

  expect_false(is_valid_solution(NULL))
  expect_false(is_valid_solution(1))
  expect_false(is_valid_solution(list()))
})
