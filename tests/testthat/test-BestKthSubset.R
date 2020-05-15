library("testthat")

#
# IHT regression
#

test_that("IHT simulated logistic regression - no intercept", {
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

    iht <- fitCyclopsModel(cyclopsData, prior = createIhtPrior(K = 2, "bic", fitBestSubset = TRUE),
                               control = createControl(noiseLevel = "silent"))

    expect_lt(sum(coef(iht) != 0.0), 3)

    # Determine MLE
    non_zero <- which(coef(iht) != 0.0)
    glm <- glm(y ~ x[,non_zero] - 1, family = binomial())
    expect_equal(as.vector(coef(iht)[which(coef(iht) != 0.0)]), as.vector(coef(glm)), tolerance = 1E-6)
})

test_that("IHT simulated logistic regression - with intercept", {
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

  cyclopsData <- createCyclopsData(y ~ x,modelType = "lr")

  expect_warning(
    iht <- fitCyclopsModel(cyclopsData, prior = createIhtPrior(K = 2, "bic", fitBestSubset = TRUE),
                           control = createControl(noiseLevel = "silent"))
  )

  expect_lt(sum(coef(iht) != 0.0), 4)

  # Determine MLE
  non_zero <- which(coef(iht) != 0.0)
  glm <- glm(y ~ x[,non_zero[-1] - 1], family = binomial())
  expect_equal(as.vector(coef(iht)[which(coef(iht) != 0.0)]), as.vector(coef(glm)), tolerance = 1E-6)
})

test_that("Fast IHT simulated logistic regression - with intercept", {
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

  cyclopsData <- createCyclopsData(y ~ x,modelType = "lr")

  expect_warning(
    ihtFast <- fitCyclopsModel(cyclopsData,
                               prior = createFastIhtPrior(K = 2, "bic", fitBestSubset = TRUE),
                               control = createControl(noiseLevel = "silent"))
  )

  expect_warning(
    ihtSlow <- fitCyclopsModel(cyclopsData, forceNewObject = TRUE,
                               prior = createIhtPrior(K = 2, "bic", fitBestSubset = TRUE),
                               control = createControl(noiseLevel = "silent"))
  )

  expect_lt(sum(coef(iht) != 0.0), 4)

  # Determine MLE
  non_zero <- which(coef(iht) != 0.0)
  glm <- glm(y ~ x[,non_zero[-1] - 1], family = binomial())
  expect_equal(as.vector(coef(iht)[which(coef(iht) != 0.0)]), as.vector(coef(glm)), tolerance = 1E-6)
})

# test_that("IHT simulated logistic regression - no convergence", {
#   set.seed(666)
#   p <- 20
#   n <- 1000
#
#   beta1 <- c(0.5, 0, 0, -1, 1.2)
#   beta2 <- seq(0, 0, length = p - length(beta1))
#   beta <- c(beta1,beta2)
#
#   x <- matrix(rnorm(p * n, mean = 0, sd = 1), ncol = p)
#
#   exb <- exp(x %*% beta)
#   prob <- exb / (1 + exb)
#   y <- rbinom(n, 1, prob)
#
#   cyclopsData <- createCyclopsData(y ~ x,modelType = "lr")
#
#   expect_error(
#     iht <- fitCyclopsModel(cyclopsData, prior = createIhtPrior("bic", exclude = "(Intercept)", fitBestSubset = TRUE,
#                                                                maxIterations = 1),
#                            control = createControl(noiseLevel = "silent")),
#     "Algorithm did not converge"
#   )
# })
