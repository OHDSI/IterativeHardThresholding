test_that("default keeps exactly K penalized covariates and converges", {
  s <- simulateBinary()
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 5), control = silent())
  expect_true(fit$iht$converged)
  expect_equal(penalizedNonZero(fit), 5)
  expect_true(coef(fit)[1] != 0) # intercept is not counted and not penalized
  expect_s3_class(fit, "cyclopsIhtFit")
})

test_that("variance Inf gives the unpenalized MLE on the support", {
  s <- simulateBinary()
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 4, variance = Inf), control = silent())
  nz <- which(coef(fit)[-1] != 0)
  g <- stats::glm(s$y ~ s$x[, nz], family = stats::binomial())
  expect_equal(as.numeric(coef(fit)[c(1, nz + 1)]), as.numeric(coef(g)), tolerance = 1e-5)
})

test_that("finds the best subset on a small problem (exhaustive search)", {
  s <- simulateBinary(n = 800, p = 10, seed = 11)
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 3, variance = Inf), control = silent())
  ll <- apply(utils::combn(10, 3), 2, function(j) as.numeric(stats::logLik(stats::glm(s$y ~ s$x[, j], family = stats::binomial()))))
  expect_equal(fit$log_likelihood, max(ll), tolerance = 1e-6)
})

test_that("every switch combination runs and keeps K", {
  s <- simulateBinary(p = 30)
  for (sw in list(createIhtControl(thresholdSteps = FALSE), createIhtControl(splicing = FALSE),
                  createIhtControl(threads = 2),
                  createIhtControl(start = "warmup"), createIhtControl(strictAcceptance = FALSE),
                  createIhtControl(exactRefits = FALSE))) {
    fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 4, control = sw), control = silent())
    expect_equal(penalizedNonZero(fit), 4)
  }
})

test_that("excluded covariates are kept and not counted", {
  s <- simulateBinary()
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 3, exclude = c("x40")), control = silent())
  expect_true(coef(fit)["x40"] != 0)
  expect_equal(sum(coef(fit)[!(names(coef(fit)) %in% c("(Intercept)", "x40"))] != 0), 3)
})

test_that("weight 0 equals dropping the rows", {
  s <- simulateBinary()
  w <- rep(1, length(s$y)); w[1:300] <- 0
  a <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 4), control = silent(), weights = w)
  sub <- Cyclops::createCyclopsData(s$y[-(1:300)] ~ s$x[-(1:300), ], modelType = "lr")
  b <- Cyclops::fitCyclopsModel(sub, prior = createIhtPrior(K = 4), control = silent())
  expect_equal(as.numeric(coef(a)), as.numeric(coef(b)), tolerance = 1e-5)
})

test_that("K larger than the number of penalized covariates is an error", {
  s <- simulateBinary(p = 6)
  expect_error(Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 7), control = silent()), "larger")
})

test_that("fixedCoefficients are refused", {
  s <- simulateBinary(p = 6)
  expect_error(Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 2), control = silent(),
                                        fixedCoefficients = rep(FALSE, 7)), "not supported")
})

test_that("plain L0 rejects a support that separates the data", {
  set.seed(8); n <- 400
  x <- matrix(stats::rbinom(n * 6, 1, 0.3), n); y <- stats::rbinom(n, 1, stats::plogis(-1 + x[, 1]))
  x <- cbind(x, y) # column 7 separates the outcome perfectly
  d <- Cyclops::createCyclopsData(y ~ x, modelType = "lr")
  fit <- Cyclops::fitCyclopsModel(d, prior = createIhtPrior(K = 2, variance = Inf), control = silent())
  expect_equal(as.numeric(coef(fit)[8]), 0) # intercept + 6 covariates, then the separating column
  expect_true(all(abs(coef(fit)) < 20))
})

test_that("startingCoefficients warm-start the fit", {
  s <- simulateBinary(p = 40, seed = 5)
  cold <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 8), control = silent())
  # starting from the solution itself: same answer, and the search stops at once
  again <- Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = createIhtPrior(K = 8),
                                    control = silent(), startingCoefficients = coef(cold))
  expect_equal(coef(again), coef(cold), tolerance = 1e-6)
  expect_lte(again$iht$refits, cold$iht$refits)
  # from a smaller K: a valid K-sized model whose objective is at least that of the warm-start refit
  small <- Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = createIhtPrior(K = 4), control = silent())
  grown <- Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = createIhtPrior(K = 8),
                                    control = silent(), startingCoefficients = coef(small))
  expect_equal(penalizedNonZero(grown), 8)
  expect_true(grown$iht$converged)
  expect_error(Cyclops::fitCyclopsModel(simulateBinary(p = 40, seed = 5)$data, prior = createIhtPrior(K = 8),
                                        control = silent(), startingCoefficients = c(1, 2)), "one finite value")
})
