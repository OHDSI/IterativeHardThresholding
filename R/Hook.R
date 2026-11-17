# @file Hook.R
#
# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of IterativeHardThresholding
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Called by Cyclops::fitCyclopsModel(cyclopsData, prior = createIhtPrior(...)).
# With returnCvPredictions = TRUE, Cyclops calls it once per fold itself.
ihtFitHook <- function(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                       startingCoefficients, fixedCoefficients) {
  if (!is.null(fixedCoefficients)) {
    stop("fixedCoefficients are not supported with an IHT prior; ",
         "use exclude for unpenalized covariates")
  }
  if (Cyclops::getNumberOfRows(cyclopsData) < 1) stop("No rows in the data")
  # Settings are read from the prior object (they may have been through a JSON round trip)
  K <- checkK(prior$K)
  variance <- checkVariance(prior$variance)
  ctl <- resolveControl(structure(as.list(prior$control), class = "ihtControl"), prior$algorithm)
  if (is.null(ctl$maxSwap) || identical(ctl$maxSwap, "NA")) ctl$maxSwap <- NA
  # prior$keepControl = TRUE: Cyclops passes the control as given (it sets coordinate descent
  # for other hooks). The refits use IHT's own threads.
  control$threads <- if (is.null(ctl$threads)) 1 else ctl$threads

  unpenalized <- unpenalizedIndices(cyclopsData, prior$exclude, prior$forceIntercept)
  nPenalized <- Cyclops::getNumberOfCovariates(cyclopsData) - length(unpenalized)
  if (K > nPenalized) {
    stop("K (", K, ") is larger than the number of penalized covariates (", nPenalized, ")")
  }

  pb <- newProblem(cyclopsData, unpenalized, variance, control, weights, prior$forceIntercept)
  res <- runEngine(pb, K, ctl, warm = warmStart(pb, startingCoefficients, ctl))
  fit <- finalFit(pb, res$b)
  class(fit) <- c(class(fit), "cyclopsIhtFit")
  fit$iht <- list(algorithm = prior$algorithm, K = K, variance = variance, objective = res$F,
                  converged = res$converged, iterations = res$iterations, steps = res$steps,
                  splices = res$splices, refits = pb$refits,
                  failedRefits = pb$failedRefits, support = pb$ids[res$T])
  fit$ihtConverged <- res$converged
  fit$ihtIterations <- res$iterations
  if (!res$converged) {
    warning("IHT did not converge in ", ctl$maxIterations, " iterations")
  }
  fit
}

# startingCoefficients (one per covariate, as fitCyclopsModel takes them) as a warm start, e.g. an
# earlier solution at a smaller K or another variance. Its non-zero penalized coefficients are the
# starting support (the K largest if there are more; the strongest outside covariates are added if
# fewer). Only the exact engine uses it; "chu2020" starts as the paper does.
warmStart <- function(pb, startingCoefficients, ctl) {
  if (is.null(startingCoefficients) || !ctl$exactRefits) return(NULL)
  b <- as.numeric(startingCoefficients)
  if (length(b) != pb$p || anyNA(b)) {
    stop("startingCoefficients must have one finite value per covariate (", pb$p, ")")
  }
  list(b = b, T = pb$pen[b[pb$pen] != 0])
}

runEngine <- function(pb, K, ctl, warm = NULL) {
  if (ctl$exactRefits) exactEngine(pb, K, ctl, warm = warm) else gradientEngine(pb, K, ctl)
}

# A Cyclops fit object at coefficients b (all fixed), so coef(), predict() and the interface
# agree with b
finalFit <- function(pb, b) {
  Cyclops::fitCyclopsModel(pb$data, prior = pb$prior, control = pb$control,
                           weights = pb$weights, startingCoefficients = b,
                           fixedCoefficients = rep(TRUE, pb$p), warnings = FALSE)
}

# Coefficient names in column order ("(Intercept)" for the intercept), from evaluateLikelihood
coefficientNames <- function(cyclopsData) {
  zero <- rep(0, Cyclops::getNumberOfCovariates(cyclopsData))
  names(Cyclops::evaluateLikelihood(cyclopsData, zero)$gradient)
}

interceptIndex <- function(cyclopsData) which(coefficientNames(cyclopsData) == "(Intercept)")

# exclude may hold covariate ids (numeric or integer64) or coefficient names
unpenalizedIndices <- function(cyclopsData, exclude, forceIntercept) {
  idx <- if (forceIntercept) integer(0) else interceptIndex(cyclopsData)
  if (length(exclude) > 0) {
    if (is.character(exclude)) {
      hit <- match(exclude, coefficientNames(cyclopsData))
    } else {
      hit <- match(as.numeric(exclude), as.numeric(Cyclops::getCovariateIds(cyclopsData)))
    }
    if (anyNA(hit)) {
      stop("exclude: covariate(s) not found: ", paste(exclude[is.na(hit)], collapse = ", "))
    }
    idx <- c(idx, hit)
  }
  sort(unique(idx))
}
