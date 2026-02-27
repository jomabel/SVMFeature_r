library(testthat)
library(SVMFeature)

test_that("SVMFeature initializes correctly", {
  data <- data.frame(x1 = 1:4, x2 = 11:14, y = 5:8)

  inputs <- c('x1', 'x2')
  output <- 'y'
  costs <- rep(1, length(inputs))
  pop_size <- 50
  num_fea <- 1
  n_iter <- 5
  max_time <- 10
  mode <- "iters"

  svm <- SVMFeature(data, inputs, output, costs, pop_size, num_fea, n_iter, max_time, mode)

  expect_s3_class(svm, "SVMFeature")
  expect_equal(svm$mode, "iters")
  expect_equal(svm$num_fea, num_fea)
})

test_that("SVMFeature fails with incorrect inputs", {

  data <- data.frame(x1 = 1:10, x2 = 11:20, y = 21:30)

  test_that("mode must be either 'iters' or 'time'", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 2),
                 pop_size = 10, num_fea = 2, mode = "wrong_mode"),
      regexp = 'Error: "mode" must be either "iters" or "time"',
      fixed = TRUE
    )
  })

  test_that("inputs must be a character vector of column names from data", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "z"), output = "y", costs = c(1, 2),
                 pop_size = 10, num_fea = 2),
      regexp = 'Error: "inputs" must be a character vector of existing column names in "data"',
      fixed = TRUE
    )
  })

  test_that("output must be a string and a column in data", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "z", costs = c(1, 2),
                 pop_size = 10, num_fea = 2),
      regexp = 'Error: "output" must be a single column name in "data"',
      fixed = TRUE
    )
  })

  test_that("pop_size must be a positive number", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 2),
                 pop_size = -10, num_fea = 2),
      regexp = 'Error: "pop_size" must be a positive number',
      fixed = TRUE
    )
  })

  test_that("n_iter must be a positive number", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 2),
                 pop_size = 10, num_fea = 2, n_iter = -10),
      regexp = 'Error: "n_iter" must be a positive number',
      fixed = TRUE
    )
  })

  test_that("num_fea must be a positive number", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 2),
                 pop_size = 10, num_fea = -2),
      regexp = 'Error: "num_fea" must be a positive number',
      fixed = TRUE
    )
  })

  test_that("max_time must be a positive number", {
    expect_error(
      SVMFeature(data, inputs = c("x1", "x2"), output = "y", costs = c(1, 2),
                 pop_size = 10, num_fea = 2, max_time = -300),
      regexp = 'Error: "max_time" must be a positive number',
      fixed = TRUE
    )
  })

})

test_that("SVMFeature creates a population in a max number of iterations", {

  set.seed(123)
  data <- data.frame(
    index = 1:100,  # Sequential index from 1 to 100
    y1 = sample(0:1, 100, replace = TRUE),  # Binary target variable (0 or 1)
    matrix(runif(100 * 278, min = -1, max = 1), nrow = 100, ncol = 278)  # Random values between -1 and 1
  )
  colnames(data) <- c("index", "y1", paste0("x", 0:277))

  svm <- SVMFeature(data, inputs = paste0('x', 0:277), output = 'y1', costs = rep(1, 277),
                    pop_size = 50, num_fea = 5, n_iter = 5, max_time = 10,
                    mode = 'iters')
  svm <- run.SVMFeature(svm)

  expect_s3_class(svm, "SVMFeature")
  expect_false(is.null(svm$population))
  expect_s3_class(svm$population, "Population")
})

