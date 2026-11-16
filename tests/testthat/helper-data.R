# Simulated designs shared by the tests
simulateBinary <- function(n = 1500, p = 40, seed = 3, intercept = TRUE) {
  set.seed(seed)
  x <- matrix(stats::rbinom(n * p, 1, stats::runif(p, 0.05, 0.5)), n, byrow = TRUE)
  beta <- c(1, -1, 0.8, -0.6, 0.5, rep(0, p - 5))
  y <- stats::rbinom(n, 1, stats::plogis(-1.5 + x %*% beta))
  list(x = x, y = y, data = if (intercept) Cyclops::createCyclopsData(y ~ x, modelType = "lr")
                            else Cyclops::createCyclopsData(y ~ x - 1, modelType = "lr"))
}

silent <- function() Cyclops::createControl(noiseLevel = "silent")

penalizedNonZero <- function(fit, unpenalized = 1) sum(coef(fit)[-unpenalized] != 0)

# Old (1.0.3) test design
legacyDesign <- function(intercept) {
  set.seed(666); p <- 20; n <- 1000; beta <- c(0.5, 0, 0, -1, 1.2, rep(0, p - 5))
  x <- matrix(stats::rnorm(p * n), ncol = p); y <- stats::rbinom(n, 1, stats::plogis(x %*% beta))
  if (intercept) Cyclops::createCyclopsData(y ~ x, modelType = "lr") else Cyclops::createCyclopsData(y ~ x - 1, modelType = "lr")
}
