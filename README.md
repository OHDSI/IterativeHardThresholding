IterativeHardThresholding
=======

[![Build Status](https://github.com/ohdsi/IterativeHardThresholding/workflows/R-CMD-check/badge.svg)](https://github.com/OHDSI/IterativeHardThresholding/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/OHDSI/IterativeHardThresholding/coverage.svg?branch=main)](https://codecov.io/github/OHDSI/IterativeHardThresholding?branch=main)

IterativeHardThresholding is part of the [HADES](https://ohdsi.github.io/Hades/).

Introduction
============

IterativeHardThresholding is an `R` package for performing L_0-based regressions using `Cyclops`

Features
========

Examples
========
 * Cox's Proportional Hazards Model
 ```r
library(Cyclops)
library(IterativeHardThresholding)
library(survival)

## data dimension
p <- 20    # number of covariates
n <- 300   # sample size

## tuning parameters
lambda <- log(n)  # BAR penalty (BIC)
xi     <- 0.1     # initial ridge penalty

## Cox model parameters
true.beta <- c(1, 0.1, 0, -1, 1, rep(0, p - 5))

## simulate data from an exponential model
x        <- matrix(rnorm(p * n, mean = 0, sd = 1), ncol = p)
ti       <- rweibull(n, shape = 1, scale = exp(-x%*%true.beta))
ui       <- runif(n, 0, 10) # Controls censoring
ci       <- rweibull(n, shape = 1, scale = ui * exp(-x%*%true.beta))
survtime <- pmin(ti, ci)
delta    <- ti == survtime; mean(delta)

cyclopsData <- createCyclopsData(Surv(survtime, delta) ~ x, modelType = "cox")
ihtPrior    <- createIhtPrior(K = 3, penalty = "bic")

cyclopsFit <- fitCyclopsModel(cyclopsData,
                              prior = ihtPrior)
coef(cyclopsFit)
 ```

* Generalized Linear Model
 ```r
library(Cyclops)
library(IterativeHardThresholding)

## data dimension
p <- 20    # number of covariates
n <- 300   # sample size

## tuning parameters
lambda <- log(n)  # BAR penalty (BIC)
xi     <- 0.1     # initial ridge penalty

## logistic model parameters
itcpt     <- 0.2 # intercept
true.beta <- c(1, 0.3, 0, -1, 1, rep(0, p - 5))

## simulate data from logistic model
x <- matrix(rnorm(p * n, mean = 0, sd = 1), ncol = p)
y <- rbinom(n, 1, 1 / (1 + exp(-itcpt - x%*%true.beta)))


# fit BAR model
cyclopsData <- createCyclopsData(y ~ x, modelType = "lr")
ihtPrior    <- createIhtPrior(K  = 3, penalty = "bic", exclude = c("(Intercept)"))

cyclopsFit <- fitCyclopsModel(cyclopsData,
                              prior = ihtPrior)
coef(cyclopsFit)
 ```
Technology
============

System Requirements
===================
Requires `R` (version 3.2.0 or higher). Installation on Windows requires [RTools]( https://CRAN.R-project.org/bin/windows/Rtools/) (`devtools >= 1.12` required for RTools34, otherwise RTools33 works fine).

Dependencies
============
 * `Cyclops`

Getting Started
===============
1. On Windows, make sure [RTools](https://CRAN.R-project.org/bin/windows/Rtools/) is installed.
2. In R, use the following commands to download and install IterativeHardThresholding:

  ```r
  install.packages("devtools")
  library(devtools)
  install_github("ohdsi/Cyclops")
  install_github("ohdsi/IterativeHardThresholding")
  ```

3. To perform a L_0-based Cyclops model fit with IHT, use the following commands in R:

  ```r
  library(IterativeHardThresholding)
  cyclopsData <- createCyclopsData(formula, modelType = "modelType") ## TODO: Update
  ihtPrior    <- createIhtPrior(K = 5, penalty = "bic")
  cyclopsFit  <- fitCyclopsModel(cyclopsData, prior = ihtPrior)
  coef(cyclopsFit) #Extract coefficients
  ```

Getting Involved
================
* Package manual: [IterativeHardThresholding manual](https://raw.githubusercontent.com/OHDSI/IterativeHardThresholding/master/extras/IterativeHardThresholding.pdf)
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="../../issues">GitHub issue tracker</a> for all bugs/issues/enhancements

License
=======
IterativeHardThresholding is licensed under Apache License 2.0.

Development
===========
IterativeHardThresholding is being developed in R Studio.

### Development status

[![Build Status](https://travis-ci.org/OHDSI/IterativeHardThresholding.svg?branch=master)](https://travis-ci.org/OHDSI/IterativeHardThresholding)
[![codecov.io](https://codecov.io/github/OHDSI/IterativeHardThresholding/coverage.svg?branch=master)](https://codecov.io/github/OHDSI/IterativeHardThresholding?branch=master)

Beta

Acknowledgments
================
- This project is supported in part through the National Institutes of Health grant R01 HG006139.
