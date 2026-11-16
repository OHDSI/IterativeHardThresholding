IterativeHardThresholding 2.0.0
=============

Breaking changes: `createIhtPrior()` now runs iterative hard thresholding as intended. Results differ from 1.0.3.

- Versions up to 1.0.3 ran Broken Adaptive Ridge with the K largest coefficients kept after each fit. Coefficients
  that reached zero were fixed there, so covariates could not re-enter; with the default `penalty = "bic"` only a
  handful of covariates survived whatever K was.
- The new default maximizes the log-likelihood minus a small ridge term (`variance`, default 1; `Inf` is plain L0)
  with exactly `K` non-zero penalized covariates. It fits each support exactly (Newton), moves by gradient
  hard-thresholding steps and by swaps of weak and strong covariates, and accepts a move only if it improves the
  objective. Covariates are ranked by |gradient| and |coefficient|, as in Chu et al. (2020). Needs Cyclops >= 4.0.1.
- `algorithm = "chu2020"` runs the algorithm of Chu et al. (2020) as in MendelIHT.jl.
- `algorithm = "legacy"` runs 1.0.3 unchanged, with its arguments (`penalty`, `fitBestSubset`, ...), so earlier
  models can be reproduced. Deprecated. Without it, those arguments stop with an error.
- `K` no longer counts the intercept or `exclude`d covariates (in 1.0.3, `K = 25` kept at most 24 covariates).
- `K` and `variance` can be tuned by cross-validation (`createIhtCvControl()`); each fold runs the whole selection,
  and out-of-fold predictions can be returned (`returnCvPredictions`).
- `createIhtControl()` has a switch for each part of the algorithm.
- `createFastIhtPrior()` is removed.
- Fixed: `"aic"` and upper-case penalties were documented but stopped with an error.

IterativeHardThresholding 1.0.3
=============

Changes: Fix documentation links to Cyclops


IterativeHardThresholding 1.0.2
=============

Changes: prepare for CRAN


IterativeHardThresholding 1.0.0
=============

Changes: initial HADES version
