library(testthat)
library(SVMFeature)

test_that("SVMFeature initializes correctly", {
  data <- data.frame(x1 = 1:4, x2 = 11:14, y = 5:8)

  inputs <- c('x1', 'x2')
  output <- 'y'
  costs <- rep(1, length(inputs))
  pop_size <- 50
  num_fea <- 1
  n_iter <- 10
  max_time <- 300
  mode <- "time"

  svm <- SVMFeature(data, inputs, output, costs, pop_size, num_fea, n_iter, max_time, mode)

  expect_s3_class(svm, "SVMFeature")
  expect_equal(svm$mode, "time")
  expect_equal(svm$num_fea, num_fea)
})

test_that("SVMFeature creates a population based on time-limit", {

  set.seed(123)

  data <- data.frame(
    index = 1:6,  # Sequential index
    y = sample(c(1, -1), 6, replace = TRUE),  # Binary target variable (1 or -1)
    x1 = sample(1:15, 6, replace = TRUE),  # Random integer values for x1 (range: 1 to 15)
    x2 = sample(1:5, 6, replace = TRUE),   # Random integer values for x2 (range: 1 to 5)
    x3 = sample(1:10, 6, replace = TRUE)   # Random integer values for x3 (range: 1 to 10)
  )

  svm <- SVMFeature(data, inputs = c('x1', 'x2', 'x3'), output = "y",
                    costs = c(5, 27, 10), pop_size = 10, num_fea = 2,
                    mode = "time", n_iter = 5, max_time = 10)
  svm <- run.SVMFeature(svm)

  expect_s3_class(svm, "SVMFeature")
  expect_false(is.null(svm$population))
  expect_s3_class(svm$population, "Population")
})
