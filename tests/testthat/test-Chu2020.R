test_that("chu2020 keeps K covariates and converges", {
  s <- simulateBinary()
  fit <- fitIht(s$data, createIhtPrior(K = 5, variance = Inf, algorithm = "chu2020"))
  expect_equal(penalizedNonZero(fit), 5)
  expect_true(fit$iht$converged)
})

test_that("chu2020 with a ridge term keeps K covariates", {
  s <- simulateBinary()
  fit <- fitIht(s$data, createIhtPrior(K = 5, variance = 1, algorithm = "chu2020"))
  expect_equal(penalizedNonZero(fit), 5)
})

test_that("chu2020 reproduces MendelIHT.jl", {
  # Reference from MendelIHT.jl 1.4.11 (fit_iht with its defaults) on the same simulated data
  ref <- readRDS(test_path("mendel-fixture.rds"))
  set.seed(42)
  n <- 500
  p <- 30
  x <- matrix(stats::rbinom(n * p, 1, rep(stats::runif(p, 0.1, 0.5), each = n)), n)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + x %*% c(1, -0.8, 0.6, rep(0, p - 3))))
  d <- Cyclops::createCyclopsData(y ~ x, modelType = "lr")
  fit <- fitIht(d, createIhtPrior(K = 5, variance = Inf, algorithm = "chu2020"))
  expect_equal(as.numeric(coef(fit)), ref$coef, tolerance = 1e-8)
  expect_equal(fit$iht$iterations, as.integer(ref$iter))
  expect_equal(fit$log_likelihood, ref$logl, tolerance = 1e-8)
})

# Why the default goes beyond the paper: same objective, better optimizer

test_that("swaps reach optima the paper's algorithm and threshold steps alone miss", {
  # Plain L0 (variance Inf), K = 10, 6 simulated designs; chu2020 run to convergence. On these,
  # the default is at least as good as chu2020 and as the default without swaps on every design,
  # and strictly better on most.
  better <- 0
  for (seed in 7:12) {
    fit <- function(...) {
      s <- simulateBinary(n = 1500, p = 60, seed = seed)
      fitIht(s$data, createIhtPrior(K = 10, variance = Inf, ...))
    }
    default <- fit()
    noSwaps <- fit(control = createIhtControl(splicing = FALSE))
    chu <- fit(algorithm = "chu2020", control = createIhtControl(maxIterations = 5000))
    expect_true(chu$iht$converged)
    expect_gte(default$log_likelihood, chu$log_likelihood - 1e-6)
    expect_gte(default$log_likelihood, noSwaps$log_likelihood - 1e-6)
    better <- better + (default$log_likelihood > chu$log_likelihood + 0.1)
  }
  expect_gte(better, 5)
})

test_that("exact refits converge where the paper's algorithm needs more than 200 iterations", {
  set.seed(77)
  n <- 4000
  p <- 300
  x <- matrix(stats::rbinom(n * p, 1, rep(stats::runif(p, 0.003, 0.05), each = n)), n)
  beta <- numeric(p)
  beta[sample(p, 10)] <- sample(c(-1, 1), 10, TRUE) * stats::runif(10, 0.2, 0.8)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + x %*% beta))
  d <- function() Cyclops::createCyclopsData(y ~ x, modelType = "lr")
  expect_warning(
    chu <- fitIht(d(), createIhtPrior(K = 60, variance = Inf, algorithm = "chu2020")),
    "did not converge"
  )
  default <- fitIht(d(), createIhtPrior(K = 60, variance = Inf))
  expect_false(chu$iht$converged) # stops at MendelIHT's default cap of 200 iterations
  expect_true(default$iht$converged)
  expect_gt(default$log_likelihood, chu$log_likelihood)
})
