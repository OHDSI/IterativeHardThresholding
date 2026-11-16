test_that("Poisson with an offset", {
  set.seed(4); n <- 1500; p <- 20
  x <- matrix(stats::rbinom(n * p, 1, 0.3), n); time <- stats::runif(n, 0.5, 2)
  y <- stats::rpois(n, time * exp(-1 + x %*% c(0.8, -0.7, 0.5, rep(0, p - 3))))
  d <- Cyclops::createCyclopsData(y ~ x + offset(log(time)), modelType = "pr")
  fit <- Cyclops::fitCyclopsModel(d, prior = createIhtPrior(K = 3, variance = Inf), control = silent())
  nz <- which(coef(fit)[-1] != 0)
  expect_equal(length(nz), 3)
  g <- stats::glm(y ~ x[, nz] + offset(log(time)), family = stats::poisson())
  expect_equal(as.numeric(coef(fit)[c(1, nz + 1)]), as.numeric(coef(g)), tolerance = 1e-5)
})

test_that("Cox", {
  set.seed(5); n <- 1000; p <- 20
  x <- matrix(stats::rnorm(n * p), n); t <- stats::rexp(n, exp(x %*% c(0.8, -0.6, 0.4, rep(0, p - 3))))
  cens <- stats::runif(n, 0, 3); time <- pmin(t, cens); status <- as.numeric(t <= cens)
  d <- Cyclops::createCyclopsData(survival::Surv(time, status) ~ x, modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(d, prior = createIhtPrior(K = 3, variance = Inf), control = silent())
  nz <- which(coef(fit) != 0)
  expect_equal(length(nz), 3)
  g <- survival::coxph(survival::Surv(time, status) ~ x[, nz], ties = "breslow")
  expect_equal(as.numeric(coef(fit)[nz]), as.numeric(coef(g)), tolerance = 1e-4)
})
