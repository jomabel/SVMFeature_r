kappa_test_result <- function(counts, objective) {
  # Counts are ordered TN, FP, FN, TP; expected kappas are independent oracles.
  y_true <- rep(c(-1, -1, 1, 1), counts)
  y_pred <- rep(c(-1, 1, -1, 1), counts)
  calculate_metrics(solution(1, 1, objective), y_true, y_pred, y_pred)
}

test_that("Cohen kappa is one for perfect agreement", {
  for (objective in c("DIST-EPS", "MC")) {
    for (counts in list(c(3, 0, 0, 3), c(7, 0, 0, 2))) {
      expect_equal(kappa_test_result(counts, objective)$kappa, 1)
    }
  }
})

test_that("Cohen kappa is zero for chance-level agreement", {
  for (objective in c("DIST-EPS", "MC")) {
    for (counts in list(c(10, 10, 10, 10), c(18, 6, 6, 2))) {
      expect_equal(kappa_test_result(counts, objective)$kappa, 0)
    }
  }
})

test_that("Cohen kappa matches known values for mixed class proportions", {
  for (objective in c("DIST-EPS", "MC")) {
    # TN=20, FP=20, FN=10, TP=30 gives kappa=1/4.
    expect_equal(kappa_test_result(c(20, 20, 10, 30), objective)$kappa, 1 / 4)
    # TN=70, FP=10, FN=10, TP=10 gives kappa=3/8.
    expect_equal(kappa_test_result(c(70, 10, 10, 10), objective)$kappa, 3 / 8)
  }
})

test_that("Cohen kappa preserves negative agreement values", {
  for (objective in c("DIST-EPS", "MC")) {
    expect_equal(kappa_test_result(c(0, 20, 20, 0), objective)$kappa, -1)
    expect_equal(kappa_test_result(c(4, 16, 16, 4), objective)$kappa, -3 / 5)
  }
})

test_that("constant predictions give zero kappa when both classes occur", {
  for (objective in c("DIST-EPS", "MC")) {
    for (counts in list(c(3, 0, 3, 0), c(0, 3, 0, 3))) {
      expect_equal(kappa_test_result(counts, objective)$kappa, 0)
    }
  }
})

test_that("small-sample point estimates do not emit interval warnings", {
  for (objective in c("DIST-EPS", "MC")) {
    # These n=4 and n=6 cases formerly emitted confidence-interval warnings.
    expect_warning(result_n4 <- kappa_test_result(c(2, 0, 1, 1), objective), NA)
    expect_equal(result_n4$kappa, 1 / 2)
    expect_warning(result_n6 <- kappa_test_result(c(3, 0, 1, 2), objective), NA)
    expect_equal(result_n6$kappa, 2 / 3)
  }
})

test_that("perfect single-class agreement has explicitly undefined kappa", {
  for (objective in c("DIST-EPS", "MC")) {
    for (counts in list(c(6, 0, 0, 0), c(0, 0, 0, 6))) {
      expect_warning(result <- kappa_test_result(counts, objective), NA)
      expect_identical(result$kappa, NA_real_)
    }
  }
})

test_that("empty inputs have explicitly undefined kappa without warnings", {
  for (objective in c("DIST-EPS", "MC")) {
    expect_warning(result <- kappa_test_result(c(0, 0, 0, 0), objective), NA)
    expect_identical(result$kappa, NA_real_)
  }
})
