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

# The fitHook of both priors. Defined here, so that its environment is the package namespace:
# a closure made inside createIhtPrior() would carry that call's frame into saved settings (JSON).
# It looks the hook up by name, so settings saved earlier run the installed version.
# With returnCvPredictions = TRUE, Cyclops calls it once per fold itself.
fitHook <- function(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                    startingCoefficients, fixedCoefficients, ...) {
  name <- if (identical(prior$algorithm, "legacy")) "legacyHook" else "ihtFitHook"
  hook <- utils::getFromNamespace(name, "IterativeHardThresholding")
  hook(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
       startingCoefficients, fixedCoefficients)
}

ihtFitHook <- function(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                       startingCoefficients, fixedCoefficients) {
  if (!is.null(fixedCoefficients)) {
    stop("fixedCoefficients are not supported with an IHT prior; ",
         "use exclude for unpenalized covariates")
  }
  # (Cyclops checks for rows only after it hands the fit to the hook)
  if (Cyclops::getNumberOfRows(cyclopsData) < 1) stop("No rows in the data")
  # Settings are read from the prior object (they may have been through a JSON round trip)
  K <- checkK(prior$K)
  variance <- checkVariance(prior$variance)
  ctl <- resolveControl(structure(as.list(prior$control), class = "ihtControl"), prior$algorithm)
  # prior$keepControl = TRUE: Cyclops passes the control as given (it sets coordinate descent
  # for other hooks). The refits use IHT's own threads.
  control$threads <- ctl$threads

  names <- coefficientNames(cyclopsData)
  intercept <- which(names == "(Intercept)")
  unpenalized <- unpenalizedIndices(cyclopsData, names, prior$exclude, prior$forceIntercept)
  # an intercept that is forced into the penalty but also excluded is simply not penalized
  forceIntercept <- prior$forceIntercept && !any(intercept %in% unpenalized)
  nPenalized <- length(names) - length(unpenalized)
  if (K > nPenalized) {
    stop("K (", K, ") is larger than the number of penalized covariates (", nPenalized, ")")
  }

  pb <- newProblem(cyclopsData, unpenalized, variance, control, weights, forceIntercept,
                   intercept)
  warm <- warmStart(pb, startingCoefficients, ctl)
  res <- if (ctl$exactRefits) exactEngine(pb, K, ctl, warm = warm) else gradientEngine(pb, K, ctl)
  fit <- finalFit(pb, res$b)
  class(fit) <- c(class(fit), "cyclopsIhtFit")
  fit$iht <- list(algorithm = prior$algorithm, K = K, variance = variance,
                  objective = res$value, converged = res$converged,
                  iterations = res$iterations, steps = res$steps, splices = res$splices,
                  refits = pb$refits, failedRefits = pb$failedRefits,
                  support = pb$ids[res$support])
  # the names of version 1.0.3
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
  if (length(b) != pb$p || !all(is.finite(b))) {
    stop("startingCoefficients must have one finite value per covariate (", pb$p, ")")
  }
  list(b = b, support = pb$pen[b[pb$pen] != 0])
}

# A Cyclops fit object at coefficients b (all fixed), so coef(), predict() and the interface
# agree with b (rejected trials leave the data object at other coefficients)
finalFit <- function(pb, b) {
  Cyclops::fitCyclopsModel(pb$data, prior = pb$prior, control = pb$control,
                           weights = pb$weights, startingCoefficients = b,
                           fixedCoefficients = rep(TRUE, pb$p), warnings = FALSE)
}

# Coefficient names in column order ("(Intercept)" for the intercept). Cyclops has no exported
# accessor for them; evaluateLikelihood names its gradient.
coefficientNames <- function(cyclopsData) {
  zero <- rep(0, Cyclops::getNumberOfCovariates(cyclopsData))
  names(Cyclops::evaluateLikelihood(cyclopsData, zero)$gradient)
}

# exclude may hold covariate ids (numeric or integer64) or coefficient names
unpenalizedIndices <- function(cyclopsData, names, exclude, forceIntercept) {
  idx <- if (forceIntercept) integer(0) else which(names == "(Intercept)")
  if (length(exclude) > 0) {
    if (is.character(exclude)) {
      hit <- match(exclude, names)
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
