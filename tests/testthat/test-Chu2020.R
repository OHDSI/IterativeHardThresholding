test_that("chu2020 keeps K covariates and converges", {
  s <- simulateBinary()
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 5, variance = Inf, algorithm = "chu2020"),
                                  control = silent())
  expect_equal(penalizedNonZero(fit), 5)
  expect_true(fit$iht$converged)
})

test_that("chu2020 is deterministic", {
  pr <- createIhtPrior(K = 8, variance = Inf, algorithm = "chu2020")
  a <- Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = pr, control = silent())
  b <- Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = pr, control = silent())
  expect_equal(coef(a), coef(b))
})

test_that("chu2020 reproduces MendelIHT.jl", {
  # Reference from MendelIHT.jl 1.4.11 (fit_iht with its defaults) on the same simulated data
  ref <- readRDS(test_path("mendel-fixture.rds"))
  set.seed(42); n <- 500; p <- 30
  x <- matrix(stats::rbinom(n * p, 1, rep(stats::runif(p, 0.1, 0.5), each = n)), n)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + x %*% c(1, -0.8, 0.6, rep(0, p - 3))))
  d <- Cyclops::createCyclopsData(y ~ x, modelType = "lr")
  fit <- Cyclops::fitCyclopsModel(d, prior = createIhtPrior(K = 5, variance = Inf, algorithm = "chu2020"), control = silent())
  expect_equal(as.numeric(coef(fit)), ref$coef, tolerance = 1e-8)
  expect_equal(fit$iht$iterations, as.integer(ref$iter))
  expect_equal(fit$log_likelihood, ref$logl, tolerance = 1e-8)
})
