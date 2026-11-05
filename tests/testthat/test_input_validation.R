input_validation_data <- function() {
  data.frame(
    x1 = c(1, 2, 3, 5, 6, 7),
    x2 = c(6L, 5L, 4L, 3L, 2L, 1L),
    y = c(-1, -1, -1, 1, 1, 1)
  )
}

input_validation_model <- function(...) {
  args <- list(data = input_validation_data(), inputs = c("x1", "x2"),
               output = "y", pop_size = 3, num_fea = 1, n_iter = 1)
  supplied <- list(...)
  args[names(supplied)] <- supplied
  do.call(svm_fs_meta, args)
}

test_that("valid data and boundary sizes construct a model", {
  model <- input_validation_model()
  expect_s3_class(model, "svm_fs_meta")
  expect_equal(model$pop_size, 3)
  expect_equal(model$n_iter, 1)
  expect_equal(model$num_fea, 1)
  expect_equal(input_validation_model(num_fea = 2)$num_fea, 2)

  two_per_class <- input_validation_data()[c(1, 2, 4, 5), ]
  expect_s3_class(input_validation_model(data = two_per_class), "svm_fs_meta")
  expect_s3_class(input_validation_model(max_time = 0.01), "svm_fs_meta")
})

test_that("data must be a nonempty data frame with unambiguous names", {
  d <- input_validation_data()
  duplicate_names <- empty_names <- missing_names <- d
  names(duplicate_names) <- c("x1", "x1", "y")
  names(empty_names)[1] <- ""
  names(missing_names)[1] <- NA_character_
  invalid_data <- list(NULL, as.matrix(d), unclass(d), d[FALSE, ], d[, FALSE])
  for (value in invalid_data) {
    expect_error(input_validation_model(data = value), "data")
  }
  expect_error(input_validation_model(data = duplicate_names, inputs = "x1"), "data")
  for (value in list(empty_names, missing_names)) {
    expect_error(input_validation_model(data = value, inputs = "x2"), "data")
  }
})

test_that("data frame subclasses are normalized without altering the caller", {
  d <- input_validation_data()
  class(d) <- c("test_frame", "data.frame")
  original <- d
  model <- input_validation_model(data = d)
  expect_identical(class(model$data), "data.frame")
  expect_identical(d, original)
  expect_equal(model$data$x1, d$x1)
})

test_that("inputs must select distinct predictors and exclude the response", {
  invalid_inputs <- list(
    NULL, character(), c("x1", "x1"), c("x1", "absent"),
    c("x1", "y"), "", NA_character_, c("x1", NA_character_),
    1:2, factor("x1"), list("x1"), matrix(c("x1", "x2"), nrow = 1)
  )
  for (value in invalid_inputs) {
    expect_error(input_validation_model(inputs = value), "inputs")
  }
})

test_that("selected predictors must be finite real numeric columns", {
  invalid_columns <- list(
    as.character(1:6), factor(1:6), rep(TRUE, 6),
    as.complex(1:6), c(NA_real_, 2:6), c(NaN, 2:6),
    c(Inf, 2:6), c(-Inf, 2:6),
    I(as.list(1:6)), I(matrix(1:12, nrow = 6))
  )
  for (value in invalid_columns) {
    d <- input_validation_data()
    d$x1 <- value
    expect_error(input_validation_model(data = d), "inputs")
  }
})

test_that("unused columns need not satisfy predictor requirements", {
  d <- input_validation_data()
  d$unused_text <- letters[1:6]
  d$unused_missing <- rep(NA_real_, 6)
  model <- input_validation_model(data = d)
  expect_s3_class(model, "svm_fs_meta")
  expect_identical(model$inputs, c("x1", "x2"))
  expect_identical(model$data$unused_text, d$unused_text)
  expect_identical(model$data$unused_missing, d$unused_missing)
})

test_that("output must name exactly one existing response column", {
  invalid_outputs <- list(
    NULL, character(), c("y", "x1"), "absent", "", NA_character_,
    1, TRUE, factor("y"), list("y"), matrix("y", nrow = 1)
  )
  for (value in invalid_outputs) {
    expect_error(input_validation_model(output = value), "output")
  }
})

test_that("response values must encode both classes with enough observations", {
  invalid_responses <- list(
    as.character(c(-1, -1, -1, 1, 1, 1)),
    factor(c(-1, -1, -1, 1, 1, 1)), rep(TRUE, 6),
    as.complex(c(-1, -1, -1, 1, 1, 1)),
    c(NA_real_, -1, -1, 1, 1, 1), c(NaN, -1, -1, 1, 1, 1),
    c(Inf, -1, -1, 1, 1, 1), c(-Inf, -1, -1, 1, 1, 1),
    c(0, 0, 0, 1, 1, 1), c(-1, -1, 0, 1, 1, 1),
    rep(1, 6), rep(-1, 6), c(-1, 1, 1, 1, 1, 1),
    c(-1, -1, -1, -1, -1, 1),
    I(as.list(c(-1, -1, -1, 1, 1, 1))),
    I(matrix(c(-1, -1, -1, 1, 1, 1), ncol = 1))
  )
  for (value in invalid_responses) {
    d <- input_validation_data()
    d$y <- value
    expect_error(input_validation_model(data = d), "output")
  }
})

input_validation_bad_scalars <- list(
  NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf,
  "1", TRUE, 1 + 0i, list(1), matrix(1, nrow = 1)
)

for (parameter in c("pop_size", "n_iter", "num_fea")) {
  test_that(paste(parameter, "must be a finite scalar integer in range"), {
    extra_values <- switch(parameter,
      pop_size = list(-1, 0, 1, 2, 3.5),
      n_iter = list(-1, 0, 1.5),
      num_fea = list(-1, 0, 1.5, 3)
    )
    for (value in c(input_validation_bad_scalars, extra_values)) {
      args <- setNames(list(value), parameter)
      expect_error(do.call(input_validation_model, args), parameter)
    }
  })
}

test_that("max_time must be a finite positive numeric scalar", {
  for (value in c(input_validation_bad_scalars, list(-1, 0))) {
    expect_error(input_validation_model(max_time = value), "max_time")
  }
})

for (parameter in c("p_mutation", "p_mut_ind", "p_mut_fea", "p_mut_coord")) {
  test_that(paste(parameter, "must be a numeric probability"), {
    for (value in c(input_validation_bad_scalars, list(-0.01, 1.01))) {
      args <- setNames(list(value), parameter)
      expect_error(do.call(input_validation_model, args), parameter)
    }
    for (value in c(0, 1)) {
      args <- setNames(list(value), parameter)
      model <- do.call(input_validation_model, args)
      expect_equal(model$population[[parameter]], value)
    }
  })
}

test_that("all mutation probabilities can be zero at construction", {
  model <- input_validation_model(p_mutation = 0, p_mut_ind = 0,
                                  p_mut_fea = 0, p_mut_coord = 0)
  expect_s3_class(model, "svm_fs_meta")
})

for (parameter in c("clones", "mut_coord")) {
  test_that(paste(parameter, "accepts only scalar binary flags"), {
    invalid_flags <- c(input_validation_bad_scalars[-9], list(-1, 0.5, 2, NA))
    for (value in invalid_flags) {
      args <- setNames(list(value), parameter)
      expect_error(do.call(input_validation_model, args), parameter)
    }
    for (value in list(0, 1, FALSE, TRUE)) {
      args <- setNames(list(value), parameter)
      model <- do.call(input_validation_model, args)
      expect_equal(as.numeric(model$population[[parameter]]), as.numeric(value))
    }
  })
}

for (parameter in c("mode", "objective")) {
  test_that(paste(parameter, "must be one supported character value"), {
    valid_option <- if (parameter == "mode") "iters" else "DIST-EPS"
    invalid_options <- list(NULL, character(), NA_character_, "", "unknown",
                            c("iters", "time"), TRUE, 1, factor("iters"), list("iters"),
                            matrix(valid_option, nrow = 1))
    for (value in invalid_options) {
      args <- setNames(list(value), parameter)
      expect_error(do.call(input_validation_model, args), parameter)
    }
  })
}

test_that("every documented mode and objective combination is accepted", {
  for (mode in c("iters", "time")) {
    for (objective in c("DIST-EPS", "MC")) {
      model <- input_validation_model(mode = mode, objective = objective)
      expect_identical(model$mode, mode)
      expect_identical(model$objective, objective)
    }
  }
})
