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

library(ParallelLogger)

test_that("JSON serialization", {
  prior <- createIhtPrior(K = 2)
  fileName <- tempfile()
  ParallelLogger::saveSettingsToJson(prior, fileName)
  rm(prior)
  restoredPrior <- ParallelLogger::loadSettingsFromJson(fileName)
  unlink(fileName)

  set.seed(666)
  p <- 20
  n <- 1000

  beta1 <- c(0.5, 0, 0, -1, 1.2)
  beta2 <- seq(0, 0, length = p - length(beta1))
  beta <- c(beta1,beta2)

  x <- matrix(rnorm(p * n, mean = 0, sd = 1), ncol = p)

  exb <- exp(x %*% beta)
  prob <- exb / (1 + exb)
  y <- rbinom(n, 1, prob)

  cyclopsData <- createCyclopsData(y ~ x - 1,modelType = "lr")

  iht <- fitCyclopsModel(cyclopsData,
                         prior = restoredPrior,
                         control = createControl(noiseLevel = "silent"))



})
