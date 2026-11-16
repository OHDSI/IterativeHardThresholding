test_that("legacy reproduces version 1.0.3", {
  ref <- readRDS(test_path("legacy-1.0.3.rds")) # coefficients from IterativeHardThresholding 1.0.3
  pr <- function(...) suppressWarnings(createIhtPrior(..., algorithm = "legacy"))
  fits <- list(
    a = suppressWarnings(Cyclops::fitCyclopsModel(legacyDesign(FALSE), prior = pr(K = 2, penalty = "bic", fitBestSubset = TRUE), control = silent())),
    b = suppressWarnings(Cyclops::fitCyclopsModel(legacyDesign(TRUE), prior = pr(K = 2, penalty = "bic", fitBestSubset = TRUE), control = silent())),
    c = suppressWarnings(Cyclops::fitCyclopsModel(legacyDesign(TRUE), prior = pr(K = 10, penalty = 1), control = silent())),
    d = suppressWarnings(Cyclops::fitCyclopsModel(legacyDesign(TRUE), prior = pr(K = 5, penalty = 0.5, delta = 0.5, initialRidgeVariance = 1), control = silent())))
  for (k in names(ref)) expect_equal(as.numeric(coef(fits[[k]])), as.numeric(ref[[k]]), tolerance = 1e-10)
})

test_that("legacy warns that it is deprecated and checks its arguments", {
  expect_warning(createIhtPrior(K = 2, algorithm = "legacy"), "deprecated")
  expect_error(suppressWarnings(createIhtPrior(K = 2, algorithm = "legacy", foo = 1)), "Unknown")
})
