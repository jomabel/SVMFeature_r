library(testthat)
library(svmfsmeta)

test_that("svmfsmeta initializes correctly", {
  data <- data.frame(x1 = 1:4, x2 = 11:14, y = c(-1, -1, 1, 1))

  inputs <- c('x1', 'x2')
  output <- 'y'
  pop_size <- 50
  num_fea <- 1
  n_iter <- 5
  max_time <- 10
  mode <- "iters"

  svm <- svmfsmeta::svm_fs_meta(data, inputs, output, pop_size, num_fea, n_iter, max_time, mode)

  expect_s3_class(svm, "svm_fs_meta")
  expect_equal(svm$mode, "iters")
  expect_equal(svm$num_fea, num_fea)
})

test_that("svmfsmeta fails with incorrect inputs", {

  data <- data.frame(x1 = 1:10, x2 = 11:20, y = rep(c(-1, 1), 5))

  test_that("mode must be either 'iters' or 'time'", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
                 pop_size = 10, num_fea = 2, mode = "wrong_mode"),
      regexp = '"mode" must be either "iters" or "time"',
      fixed = TRUE
    )
  })

  test_that("inputs must be a character vector of column names from data", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "z"), output = "y",
                 pop_size = 10, num_fea = 2),
      regexp = '"inputs" must be a non-empty character vector of unique existing column names in "data"',
      fixed = TRUE
    )
  })

  test_that("output must be a string and a column in data", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "z",
                 pop_size = 10, num_fea = 2),
      regexp = '"output" must name a single column in "data"',
      fixed = TRUE
    )
  })

  test_that("pop_size must be an integer of at least three", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
                 pop_size = -10, num_fea = 2),
      regexp = '"pop_size" must be a finite integer greater than or equal to 3',
      fixed = TRUE
    )
  })

  test_that("n_iter must be a positive integer", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
                 pop_size = 10, num_fea = 2, n_iter = -10),
      regexp = '"n_iter" must be a finite integer greater than or equal to 1',
      fixed = TRUE
    )
  })

  test_that("num_fea must be a positive integer", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
                 pop_size = 10, num_fea = -2),
      regexp = '"num_fea" must be a finite integer greater than or equal to 1',
      fixed = TRUE
    )
  })

  test_that("max_time must be a positive number", {
    expect_error(
      svmfsmeta::svm_fs_meta(data, inputs = c("x1", "x2"), output = "y",
                 pop_size = 10, num_fea = 2, max_time = -300),
      regexp = '"max_time" must be a positive finite number',
      fixed = TRUE
    )
  })

})

test_that("svmfsmeta creates a population in a max number of iterations", {

  set.seed(123)
  data <- data.frame(
    index = 1:100,  # Sequential index from 1 to 100
    y1 = sample(c(-1, 1), 100, replace = TRUE),  # Binary target variable (-1 or 1)
    matrix(runif(100 * 278, min = -1, max = 1), nrow = 100, ncol = 278)  # Random values between -1 and 1
  )
  colnames(data) <- c("index", "y1", paste0("x", 0:277))

  svm <- svmfsmeta::svm_fs_meta(data, inputs = paste0('x', 0:277), output = 'y1',
                    pop_size = 50, num_fea = 5, n_iter = 5, max_time = 10,
                    mode = 'iters')
  svm <- run.svm_fs_meta(svm)

  expect_s3_class(svm, "svm_fs_meta")
  expect_false(is.null(svm$population))
  expect_s3_class(svm$population, "population")
})

