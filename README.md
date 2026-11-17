IterativeHardThresholding
=======

[![Build Status](https://github.com/ohdsi/IterativeHardThresholding/workflows/R-CMD-check/badge.svg)](https://github.com/OHDSI/IterativeHardThresholding/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/OHDSI/IterativeHardThresholding/coverage.svg?branch=main)](https://app.codecov.io/github/OHDSI/IterativeHardThresholding?branch=main)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/IterativeHardThresholding)](https://CRAN.R-project.org/package=IterativeHardThresholding)
[![CRAN_Status_Badge](https://cranlogs.r-pkg.org/badges/IterativeHardThresholding)](https://cran.r-project.org/package=IterativeHardThresholding)

`IterativeHardThresholding` is part of the [HADES](https://ohdsi.github.io/Hades/).

Introduction
============

`IterativeHardThresholding` fits regression models with at most `K` non-zero coefficients (L0 regression), using
[`Cyclops`](https://github.com/OHDSI/Cyclops) for the fitting. It maximizes the log-likelihood minus a small ridge
term, with exactly `K` penalized covariates:

  maximize  log L(beta) - sum(beta_j^2) / (2 * variance)   subject to at most K non-zero beta_j

The intercept and any `exclude`d covariates are not penalized and do not count towards `K`.

Features
========

- **Exact fits on the support.** Each candidate support is fitted exactly (Newton, Cyclops >= 4.0.1), as in gradient
  hard thresholding pursuit (Yuan et al. 2018) and Newton hard-thresholding pursuit (Zhou et al. 2021). The search
  moves by gradient hard-thresholding steps and by swapping weak covariates in the support for strong ones outside
  ("splicing", Zhu et al. 2020; swaps as in Beck and Eldar 2013 and L0Learn), and accepts a move only if it improves
  the objective, so it always converges. The model is the one of Chu et al. (2020); only the optimizer is different.
- **`algorithm = "chu2020"`** runs iterative hard thresholding as in Chu et al. (2020), and gives the same results as
  their MendelIHT.jl.
- **Cross-validation** of `K` and the ridge `variance`, with the whole selection run inside each fold, and out-of-fold
  predictions (`fitCyclopsModel(..., cvFolds, returnCvPredictions = TRUE)`). Folds run in parallel on Linux and macOS
  when `createControl(threads)` allows it.
- Logistic, Poisson and Cox models (validated); other Cyclops model types run with a warning.
- **`algorithm = "legacy"`** reproduces version 1.0.3 (Broken Adaptive Ridge with a cap of `K`), so that earlier
  models can be refitted. Deprecated.

When to use it
==============

On real observational health data (8 prediction tasks), L0 models with `K = 25` predicted better than a lasso with the
same number of covariates (held-out AUC +0.007 on average, 7 of 8 tasks) and were better calibrated. With `K = 100`
the lasso was better. IHT suits small models. The default `variance = 1` is fast; tuning it by cross-validation
(`variance = "tune"`) improves calibration for large `K` at about 10-20 times the cost.

Examples
========

```r
library(Cyclops)
library(IterativeHardThresholding)

set.seed(1)
n <- 2000; p <- 50
x <- matrix(rbinom(n * p, 1, 0.3), n)
y <- rbinom(n, 1, plogis(-1 + x %*% c(1, -0.8, 0.6, rep(0, p - 3))))
cyclopsData <- createCyclopsData(y ~ x, modelType = "lr")

# a model with 5 covariates
fit <- fitCyclopsModel(cyclopsData, prior = createIhtPrior(K = 5))
coef(fit)[coef(fit) != 0]

# choose K and the ridge variance by 3-fold cross-validation, with out-of-fold predictions
fit <- fitCyclopsModel(cyclopsData,
                       prior = createIhtPrior(K = c(3, 5, 10), variance = "tune",
                                              cv = createIhtCvControl(seed = 1, returnCvPredictions = TRUE)))
fit$iht$K          # chosen K
fit$iht$cv         # cross-validation table
head(fit$cvPredictions)

# the algorithm of Chu et al. (2020)
fit <- fitCyclopsModel(cyclopsData, prior = createIhtPrior(K = 5, variance = Inf, algorithm = "chu2020"))
```

System Requirements
===================
Requires `R` (version 3.2.2 or higher) and `Cyclops` (version 4.0.1 or higher).

Getting Started
===============
1. On Windows, make sure [RTools](https://CRAN.R-project.org/bin/windows/Rtools/) is installed.
2. In R:

  ```r
  install.packages("remotes")
  remotes::install_github("OHDSI/IterativeHardThresholding")
  ```

References
==========
Chu BB, Keys KL, German CA, Zhou H, Zhou JJ, Sobel EM, Sinsheimer JS, Lange K (2020). Iterative hard thresholding in
genome-wide association studies: Generalized linear models, prior weights, and double sparsity. GigaScience 9(6),
giaa044. <https://doi.org/10.1093/gigascience/giaa044>

The default algorithm builds on:

- Bahmani S, Raj B, Boufounos PT (2013). Greedy sparsity-constrained optimization. Journal of Machine Learning
  Research 14(25), 807-841. (GraSP: exact fits on the support, for logistic regression.)
- Beck A, Eldar YC (2013). Sparsity constrained nonlinear optimization: optimality conditions and algorithms. SIAM
  Journal on Optimization 23(3), 1480-1509. <https://doi.org/10.1137/120869778> (Optimality conditions based on swaps.)
- Yuan XT, Li P, Zhang T (2018). Gradient hard thresholding pursuit. Journal of Machine Learning Research 18(166),
  1-43. (A gradient hard-thresholding step, then an exact fit on the new support.)
- Zhou S, Xiu N, Qi HD (2021). Global and quadratic convergence of Newton hard-thresholding pursuit. Journal of Machine
  Learning Research 22(12), 1-45. (Newton steps on the support.)
- Zhu J, Wen C, Zhu J, Zhang H, Wang X (2020). A polynomial algorithm for best-subset selection problem. Proceedings of
  the National Academy of Sciences 117(52), 33117-33123. <https://doi.org/10.1073/pnas.2014241117> (abess: splicing.)
- Hazimeh H, Mazumder R (2020). Fast best subset selection: coordinate descent and local combinatorial optimization
  algorithms. Operations Research 68(5), 1517-1537. <https://doi.org/10.1287/opre.2019.1919> (L0Learn: L0 with a
  ridge term, and swaps.)
- Dedieu A, Hazimeh H, Mazumder R (2021). Learning sparse classifiers: continuous and mixed integer optimization
  perspectives. Journal of Machine Learning Research 22(135), 1-47. (The same for logistic regression.)

Getting Involved
================
* Package manual: [IterativeHardThresholding manual](https://raw.githubusercontent.com/OHDSI/IterativeHardThresholding/master/extras/IterativeHardThresholding.pdf)
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="https://github.com/OHDSI/IterativeHardThresholding/issues">GitHub issue tracker</a> for all bugs/issues/enhancements

License
=======
`IterativeHardThresholding` is licensed under Apache License 2.0.

Development
===========
`IterativeHardThresholding` is being developed in R Studio.

Acknowledgments
================
- This project is supported in part through the National Institutes of Health grant R01 HG006139.
