IterativeHardThresholding 2.0.0
=============

Breaking changes: `createIhtPrior()` now runs iterative hard thresholding as described in Chu et al.
(2020). Results differ from 1.0.3.

- Versions up to 1.0.3 ran Broken Adaptive Ridge with a projection onto the K largest coefficients
  after each fit. Coefficients that reached zero were fixed there, so covariates could not re-enter;
  with the default `penalty = "bic"` the models often had fewer than K covariates.
- The new default maximizes the log-likelihood minus a small ridge term (`variance`, default 1;
  `Inf` is plain L0) with at most (normally exactly) `K` non-zero penalized covariates. It fits each
  support exactly (Newton), moves by gradient hard-thresholding steps and by swaps of weak and
  strong covariates, and accepts a move only if it improves the objective. Needs Cyclops >= 4.0.1.
- The second argument of `createIhtPrior()` is now `variance`; in 1.0.3 it was `penalty`.
- `algorithm = "chu2020"` runs the algorithm of Chu et al. (2020) as in MendelIHT.jl.
- `algorithm = "legacy"` runs 1.0.3 unchanged, with its arguments (`penalty`, `fitBestSubset`,
  ...), so earlier models can be reproduced. Deprecated. Without it, those arguments stop with an
  error.
- `K` no longer counts the intercept or `exclude`d covariates (in 1.0.3, `K = 25` kept at most 24
  covariates).
- `startingCoefficients` warm-start the search, e.g. from a solution at a smaller `K`. Choosing `K`
  and `variance` is left to the caller (PatientLevelPrediction tunes them by cross-validation).
- `createIhtControl()` has a switch for each part of the algorithm.
- `createFastIhtPrior()` is removed.
- The documentation no longer lists `"aic"` and upper-case penalties, which never worked.
  `algorithm = "legacy"` now refuses them when the prior is created instead of at fit time.

IterativeHardThresholding 1.0.3
=============

Changes: Fix documentation links to Cyclops


IterativeHardThresholding 1.0.2
=============

Changes: prepare for CRAN


IterativeHardThresholding 1.0.0
=============

Changes: initial HADES version
