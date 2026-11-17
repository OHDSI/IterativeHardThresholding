# @file Legacy.R
#
# Copyright 2022 Observational Health Data Sciences and Informatics
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
#
# @author Marc A. Suchard

# Version 1.0.3 of this package (Broken Adaptive Ridge with a projection onto the K largest
# coefficients after each fit), kept unchanged apart from how the prior object is built, so that
# earlier models can be reproduced: createIhtPrior(K, algorithm = "legacy", penalty = "bic", ...).
# Deprecated. The 1.0.3 behaviour is kept on purpose, including that K counts the intercept and
# excluded covariates.

# The arguments of version 1.0.3 and their defaults
legacyDefaults <- list(penalty = "bic", fitBestSubset = FALSE, initialRidgeVariance = 1E-1,
                       tolerance = 1E-8, maxIterations = 1E4, threshold = 1E-6, delta = 0)

createLegacyPrior <- function(K, exclude, forceIntercept, args) {
  unknown <- setdiff(names(args), names(legacyDefaults))
  if (length(unknown) > 0) {
    stop("Unknown argument(s) for algorithm = \"legacy\": ", paste(unknown, collapse = ", "))
  }
  legacy <- utils::modifyList(legacyDefaults, args)
  if (length(K) != 1 || K < 1) {
    stop("Maximum # of covariates must be >= 1")
  }
  # 1.0.3 accepted only these, and stopped with "Unhandled IHT penalty type" at fit time otherwise
  if (!(is.numeric(legacy$penalty) || identical(legacy$penalty, "bic"))) {
    stop("penalty must be a number or \"bic\"")
  }
  warning("algorithm = \"legacy\" runs IterativeHardThresholding 1.0.3 (Broken Adaptive Ridge ",
          "with a projection onto the K largest coefficients) and is deprecated; it will be ",
          "removed in a future version.", call. = FALSE)
  structure(list(K = K,
                 penalty = legacy$penalty,
                 exclude = exclude,
                 forceIntercept = forceIntercept,
                 algorithm = "legacy",
                 legacy = legacy,
                 fitHook = fitHook),
            class = c("cyclopsPrior", "ihtLegacyPrior"))
}

legacyHook <- function(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                       startingCoefficients, fixedCoefficients) {
  a <- prior$legacy
  ihtHook(a$fitBestSubset, a$initialRidgeVariance, a$tolerance, a$maxIterations, a$threshold,
          a$delta, cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
          startingCoefficients, fixedCoefficients)
}

# Below are package-private functions

ihtHook <- function(fitBestSubset,
                    initialRidgeVariance,
                    tolerance,
                    maxIterations,
                    cutoff,
                    delta,
                    cyclopsData,
                    ihtPrior,
                    control,
                    weights,
                    forceNewObject,
                    returnEstimates,
                    startingCoefficients,
                    fixedCoefficients) {

  # Getting starting values
  startFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = createIhtStartingPrior(cyclopsData,
                                                                                   exclude = ihtPrior$exclude,
                                                                                   forceIntercept = ihtPrior$forceIntercept,
                                                                                   initialRidgeVariance = initialRidgeVariance),
                                       control, weights, forceNewObject, returnEstimates, startingCoefficients, fixedCoefficients)

  priorType <- createIhtPriorType(cyclopsData, ihtPrior$exclude, ihtPrior$forceIntercept)
  include <- setdiff(c(1:Cyclops::getNumberOfCovariates(cyclopsData)), priorType$excludeIndices)

  pre_coef <- coef(startFit)
  penalty <- getPenalty(cyclopsData, ihtPrior)
  K <- ihtPrior$K

  ParallelLogger::logTrace("Initial penalty: %f", penalty)

  continue <- TRUE
  count <- 0
  converged <- FALSE

  while (continue) {
    count <- count + 1

    working_coef <- ifelse(abs(pre_coef) <= cutoff, 0.0, pre_coef)
    fixed <- working_coef == 0.0
    variance <- abs(working_coef) ^ (2 - delta) / penalty

    if (!is.null(priorType$excludeIndices)) {
      working_coef[priorType$excludeIndices] <- pre_coef[priorType$excludeIndices]
      fixed[priorType$excludeIndices] <- FALSE
      variance[priorType$excludeIndices] <- 0
    }

    prior <- Cyclops::createPrior(priorType$types, variance = variance,
                                  forceIntercept = ihtPrior$forceIntercept)

    fit <- Cyclops::fitCyclopsModel(cyclopsData,
                                    prior = prior,
                                    control, weights, forceNewObject,
                                    startingCoefficients = working_coef,
                                    fixedCoefficients = fixed)

    coef <- coef(fit)

    # IHT projection
    if (!is.null(priorType$excludeIndices)) {
      tmp_coef <- coef[-priorType$excludeIndices]
      entry <- length(tmp_coef) - K - 1
      kThLargest <- -sort(-abs(coef), partial = K)[K]
    } else {
      entry <- length(coef) - K - 1
      kThLargest <- -sort(-abs(coef), partial = K)[K]
    }

    mask <- abs(coef) >= kThLargest

    if (!is.null(priorType$excludeIndices)) {
      mask[priorType$excludeIndices] <- TRUE
    }

    coef <- coef * mask

    end <- min(10, length(variance))
    ParallelLogger::logTrace("Itr: %d", count)
    ParallelLogger::logTrace("\tVar : ", variance[1:end], capture = TRUE)
    ParallelLogger::logTrace("\tCoef: ", coef[1:end], capture = TRUE)
    ParallelLogger::logTrace("")

    if (max(abs(coef - pre_coef)) < tolerance) {
      converged <- TRUE
    } else {
      pre_coef <- coef
    }

    if (converged || count >= maxIterations) {
      continue <- FALSE
    }
  }

  if (count >= maxIterations) {
    stop(paste0('Algorithm did not converge after ',
                maxIterations, ' iterations.',
                ' Estimates may not be stable.'))
  }

  if (fitBestSubset) {
    fit <- Cyclops::fitCyclopsModel(cyclopsData, prior = createPrior("none"),
                                    control, weights, forceNewObject, fixedCoefficients = fixed)
  }

  class(fit) <- c(class(fit), "cyclopsIhtFit")
  fit$ihtConverged <- converged
  fit$ihtIterations <- count
  fit$penalty <- penalty
  fit$ihtFinalPriorVariance <- variance

  return(fit)
}

createIhtStartingPrior <- function(cyclopsData,
                                   exclude,
                                   forceIntercept,
                                   initialRidgeVariance) {

  Cyclops::createPrior("normal", variance = initialRidgeVariance, exclude = exclude, forceIntercept = forceIntercept)
}

createIhtPriorType <- function(cyclopsData,
                               exclude,
                               forceIntercept) {

  exclude <- Cyclops:::.checkCovariates(cyclopsData, exclude)

  if (Cyclops:::.cyclopsGetHasIntercept(cyclopsData) && !forceIntercept) {
    interceptId <- Cyclops:::.cyclopsGetInterceptLabel(cyclopsData)
    warn <- FALSE
    if (is.null(exclude)) {
      exclude <- c(interceptId)
      warn <- TRUE
    } else {
      if (!interceptId %in% exclude) {
        exclude <- c(interceptId, exclude)
        warn <- TRUE
      }
    }
    if (warn) {
      warning("Excluding intercept from regularization")
    }
  }

  indices <- NULL
  if (!is.null(exclude)) {
    covariateIds <- Cyclops::getCovariateIds(cyclopsData)
    indices <- which(covariateIds %in% exclude)
  }

  types <- rep("normal", Cyclops::getNumberOfCovariates(cyclopsData))
  if (!is.null(exclude)) {
    types[indices] <- "none"
  }

  list(types = types,
       excludeCovariateIds = exclude,
       excludeIndices = indices)
}

getPenalty <- function(cyclopsData, ihtPrior) {

  if (is.numeric(ihtPrior$penalty)) {
    return(ihtPrior$penalty)
  }

  if (ihtPrior$penalty == "bic") {
    return(log(Cyclops::getNumberOfRows(cyclopsData)) / 2)
  } else {
    stop("Unhandled IHT penalty type")
  }
}
