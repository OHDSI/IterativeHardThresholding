test_that("tuning K and variance chooses from the grid and returns out-of-fold predictions", {
  s <- simulateBinary()
  fit <- Cyclops::fitCyclopsModel(s$data, control = silent(),
    prior = createIhtPrior(K = c(2, 5), variance = c(0.1, Inf), cv = createIhtCvControl(seed = 1, returnCvPredictions = TRUE)))
  tab <- fit$iht$cv
  expect_equal(nrow(tab), 4)
  expect_equal(sum(tab$chosen), 1)
  expect_equal(fit$iht$K, tab$K[tab$chosen])
  expect_equal(penalizedNonZero(fit), fit$iht$K)
  expect_equal(nrow(fit$cvPredictions), length(s$y))
  expect_true(all(fit$cvPredictions$prediction > 0 & fit$cvPredictions$prediction < 1))
})

test_that("given folds are used", {
  s <- simulateBinary()
  folds <- data.frame(rowId = seq_along(s$y), fold = rep(1:4, length.out = length(s$y)))
  fit <- Cyclops::fitCyclopsModel(s$data, control = silent(),
    prior = createIhtPrior(K = 3, cv = createIhtCvControl(cvFolds = folds, returnCvPredictions = TRUE)))
  p <- fit$cvPredictions
  expect_equal(p$fold[order(p$rowId)], folds$fold)
})

test_that("out-of-fold predictions do not use the held-out rows", {
  s <- simulateBinary()
  folds <- data.frame(rowId = seq_along(s$y), fold = rep(1:3, length.out = length(s$y)))
  fit <- Cyclops::fitCyclopsModel(s$data, control = silent(),
    prior = createIhtPrior(K = 3, cv = createIhtCvControl(cvFolds = folds, returnCvPredictions = TRUE)))
  tr <- folds$fold != 1
  sub <- Cyclops::createCyclopsData(s$y[tr] ~ s$x[tr, ], modelType = "lr")
  f1 <- Cyclops::fitCyclopsModel(sub, prior = createIhtPrior(K = 3), control = silent())
  eta <- cbind(1, s$x[!tr, ]) %*% coef(f1)
  p <- fit$cvPredictions[fit$cvPredictions$fold == 1, ]
  expect_equal(p$prediction[order(p$rowId)], as.numeric(stats::plogis(eta)), tolerance = 1e-5)
})

test_that("one-standard-error rule picks a K no larger than the best", {
  s <- simulateBinary()
  a <- Cyclops::fitCyclopsModel(s$data, control = silent(), prior = createIhtPrior(K = c(2, 4, 8), cv = createIhtCvControl(seed = 2)))
  b <- Cyclops::fitCyclopsModel(simulateBinary()$data, control = silent(), prior = createIhtPrior(K = c(2, 4, 8), cv = createIhtCvControl(seed = 2, oneStandardError = TRUE)))
  expect_lte(b$iht$K, a$iht$K)
})

test_that("fitCyclopsModel passes cvFolds and returnCvPredictions to the hook", {
  s <- simulateBinary()
  folds <- data.frame(rowId = seq_along(s$y), fold = rep(1:3, length.out = length(s$y)))
  fit <- Cyclops::fitCyclopsModel(s$data, prior = createIhtPrior(K = 3), control = silent(),
                                  cvFolds = folds, returnCvPredictions = TRUE)
  p <- fit$cvPredictions
  expect_equal(nrow(p), length(s$y))
  expect_equal(p$fold[order(p$rowId)], folds$fold)
  expect_equal(penalizedNonZero(fit), 3)
})

test_that("folds in parallel give the same results as one after another", {
  skip_on_os("windows")
  s <- simulateBinary()
  folds <- data.frame(rowId = seq_along(s$y), fold = rep(1:3, length.out = length(s$y)))
  pr <- createIhtPrior(K = c(2, 4), variance = c(1, Inf), cv = createIhtCvControl(cvFolds = folds, returnCvPredictions = TRUE))
  a <- Cyclops::fitCyclopsModel(s$data, prior = pr, control = Cyclops::createControl(noiseLevel = "silent", threads = 1))
  b <- Cyclops::fitCyclopsModel(simulateBinary()$data, prior = pr, control = Cyclops::createControl(noiseLevel = "silent", threads = 2)) # CRAN checks allow 2 cores
  expect_equal(a$iht$cv, b$iht$cv)
  expect_equal(a$cvPredictions, b$cvPredictions)
  expect_equal(coef(a), coef(b))
})
