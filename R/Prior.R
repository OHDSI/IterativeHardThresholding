# @file Prior.R
#
# Copyright 2022 Observational Health Data Sciences and Informatics
#
# This file is part of BrokenAdaptiveRidge
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

#' @title Create an IHT Cyclops prior object
#'
#' @description
#' \code{createIhtPrior} creates an IHT Cyclops prior object for use with \code{\link[Cyclops]{fitCyclopsModel}}.
#'
#' @param K              Maximum # of non-zero covariates
#' @param penalty        Specifies the IHT penalty; possible values are `BIC` or `AIC` or a numeric value
#' @param exclude        A vector of numbers or covariateId names to exclude from prior
#' @param forceIntercept Logical: Force intercept coefficient into regularization
#' @param fitBestSubset  Logical: Fit final subset with no regularization
#' @param initialRidgeVariance Numeric: variance used for algorithm initiation
#' @param tolerance Numeric: maximum abs change in coefficient estimates from successive iterations to achieve convergence
#' @param maxIterations Numeric: maximum iterations to achieve convergence
#' @param threshold     Numeric: absolute threshold at which to force coefficient to 0
#' @param delta         Numeric: change from 2 in ridge norm dimension
#'
#' @examples
#' prior <- createIhtPrior(K = 10)
#'
#' @return
#' An IHT Cyclops prior object of class inheriting from
#' \code{"cyclopsPrior"} for use with \code{fitCyclopsModel}.
#'
#' @import Cyclops
#' @importFrom stats coef
#'
#' @export
createIhtPrior <- function(K,
                           penalty = "bic",
                           exclude = c(),
                           forceIntercept = FALSE,
                           fitBestSubset = FALSE,
                           initialRidgeVariance = 1E-1,
                           tolerance = 1E-8,
                           maxIterations = 1E4,
                           threshold = 1E-6,
                           delta = 0) {

    # TODO Check that penalty (and other arguments) is valid

    if (K < 1) {
      stop("Maximum # of covariates must be >= 1")
    }

    fitHook <- function(...) {
      # closure to capture IHT parameters
      ihtHook(fitBestSubset, initialRidgeVariance, tolerance,
              maxIterations, threshold, delta, ...)
    }

    structure(list(K = K,
                   penalty = penalty,
                   exclude = exclude,
                   forceIntercept = forceIntercept,
                   fitHook = fitHook),
              class = "cyclopsPrior")
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

    # TODO: Save current log-likelihood
    # TODO: Set stepSizeMultiplier (we will hack in a function to do this for prototyping)

    fit <- Cyclops::fitCyclopsModel(cyclopsData,
                                    prior = prior,
                                    control, weights, forceNewObject,
                                    startingCoefficients = working_coef,
                                    fixedCoefficients = fixed)

    coef <- coef(fit)

    # TODO: Get new log-likelihood
    # TODO: Compare new and current ll
    # TODO: If moving in right direction
    # TODO:   Then stepSizeMultiplier <- 1
    # TODO:   Else stepSizeMultiplier <- stepSizeMultier / 2

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
    return(log(Cyclops::getNumberOfRows(cyclopsData)) / 2) # TODO Handle stratified models
  } else {
    stop("Unhandled IHT penalty type")
  }
}
