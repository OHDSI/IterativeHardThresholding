test_that("arguments are checked", {
  expect_error(createIhtPrior(K = 0), "whole number")
  expect_error(createIhtPrior(K = 2.5), "whole number")
  expect_error(createIhtPrior(K = Inf), "whole number")
  expect_error(createIhtPrior(K = c(2, 5)), "one whole number")
  expect_error(createIhtPrior(K = 2, variance = -1), "positive")
  expect_error(createIhtPrior(K = 2, variance = c(0.1, 1)), "one positive number")
  expect_error(createIhtPrior(K = 2, variance = "bic"), "legacy")
  expect_error(createIhtPrior(K = 2, penalty = "bic"), "legacy")
  expect_error(createIhtPrior(K = 2, foo = 1), "Unknown")
  expect_error(createIhtPrior(K = 2, algorithm = "other"), "should be one of")
  expect_equal(createIhtPrior(K = 2, variance = "Inf")$variance, Inf)
  expect_equal(createIhtPrior(K = 2, variance = "inf")$variance, Inf)
})

test_that("control values are checked", {
  bad <- list(createIhtControl(maxSwap = 0), createIhtControl(maxHalve = -1),
              createIhtControl(maxIterations = 0), createIhtControl(tolerance = 0),
              createIhtControl(splicing = "yes"), createIhtControl(start = "other"))
  for (control in bad) expect_error(createIhtPrior(K = 2, control = control))
  expect_error(createIhtPrior(K = 2, control = createIhtControl(maxSwap = 0)), "maxSwap")
  expect_error(createIhtPrior(K = 2, control = createIhtControl(exactRefits = FALSE,
                                                                 splicing = TRUE)),
               "splicing needs")
})

test_that("presets and switches", {
  p <- createIhtPrior(K = 2)
  expect_true(p$control$exactRefits)
  expect_true(p$control$splicing)
  p <- createIhtPrior(K = 2, algorithm = "chu2020")
  expect_false(p$control$exactRefits)
  expect_false(p$control$splicing)
  expect_equal(p$control$maxHalve, 3)
  p <- createIhtPrior(K = 2, algorithm = "chu2020",
                      control = createIhtControl(maxIterations = 1000))
  expect_equal(p$control$maxIterations, 1000)
  expect_equal(createIhtPrior(K = 2)$control$start, "gradient")
  expect_equal(createIhtPrior(K = 2)$variance, 1)
})

test_that("JSON round trip", {
  prior <- createIhtPrior(K = 3, variance = Inf)
  fileName <- tempfile()
  on.exit(unlink(fileName))
  ParallelLogger::saveSettingsToJson(prior, fileName)
  # the hook does not carry createIhtPrior's frame into the settings
  expect_identical(environment(prior$fitHook), asNamespace("IterativeHardThresholding"))
  restored <- ParallelLogger::loadSettingsFromJson(fileName)
  s <- simulateBinary()
  a <- fitIht(s$data, prior)
  b <- fitIht(simulateBinary()$data, restored)
  expect_equal(coef(a), coef(b))
})
