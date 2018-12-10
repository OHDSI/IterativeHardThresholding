library("testthat")

#
# IHT regression
#

test_that("Numeric penalty", {
  prior <- createIhtPrior(penalty = 10)
  expect_equal(IterativeHardThresholding:::getPenalty(NULL, prior), 10)
})

test_that("Unhandled penalty", {
  prior <- createIhtPrior(penalty = "aic")
  expect_error(IterativeHardThresholding:::getPenalty(NULL, prior))
})
