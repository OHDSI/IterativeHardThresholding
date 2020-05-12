library("testthat")

#
# IHT regression
#

test_that("Numeric penalty", {
  prior <- createIhtPrior(K = 1, penalty = 10)
  expect_equal(IterativeHardThresholding:::getPenalty(NULL, prior), 10)
})

test_that("Unhandled penalty", {
  prior <- createIhtPrior(K = 1, penalty = "aic")
  expect_error(IterativeHardThresholding:::getPenalty(NULL, prior))
})
