test_that("legacy reproduces version 1.0.3", {
  # coefficients from IterativeHardThresholding 1.0.3
  ref <- readRDS(test_path("legacy-1.0.3.rds"))
  legacy <- function(data, ...) {
    # 1.0.3 logs every iteration
    fit <- NULL
    utils::capture.output(fit <- suppressWarnings(
      fitIht(data, createIhtPrior(..., algorithm = "legacy"))
    ))
    fit
  }
  fits <- list(
    a = legacy(legacyDesign(FALSE), K = 2, penalty = "bic", fitBestSubset = TRUE),
    b = legacy(legacyDesign(TRUE), K = 2, penalty = "bic", fitBestSubset = TRUE),
    c = legacy(legacyDesign(TRUE), K = 10, penalty = 1),
    d = legacy(legacyDesign(TRUE), K = 5, penalty = 0.5, delta = 0.5, initialRidgeVariance = 1)
  )
  for (k in names(ref)) {
    expect_equal(as.numeric(coef(fits[[k]])), as.numeric(ref[[k]]), tolerance = 1e-10)
  }
})

test_that("legacy warns that it is deprecated and checks its arguments", {
  expect_warning(createIhtPrior(K = 2, algorithm = "legacy"), "deprecated")
  expect_error(suppressWarnings(createIhtPrior(K = 2, algorithm = "legacy", foo = 1)), "Unknown")
  expect_error(suppressWarnings(createIhtPrior(K = 2, algorithm = "legacy", penalty = "aic")),
               "penalty")
  # the second argument of 1.0.3 (penalty) is not silently dropped
  expect_error(suppressWarnings(createIhtPrior(10, 1, algorithm = "legacy")), "by name")
})
