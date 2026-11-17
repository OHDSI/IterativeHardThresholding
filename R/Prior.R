# @file Prior.R
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

#' @title Create an IHT Cyclops prior object
#'
#' @description
#' \code{createIhtPrior} creates a prior for \code{\link[Cyclops]{fitCyclopsModel}} that keeps at most \code{K}
#' non-zero penalized coefficients. It maximizes the log-likelihood minus a small ridge term,
#' \deqn{\ell(\beta) - \sum_{j \in P} \beta_j^2 / (2 \sigma^2) \quad \textrm{subject to at most } K
#'   \textrm{ non-zero } \beta_j, j \in P,}
#' where \eqn{P} are the penalized covariates (all but the intercept and \code{exclude}) and \eqn{\sigma^2} is
#' \code{variance}. \code{variance = Inf} is plain L0 (best subset of size \code{K}).
#'
#' The default algorithm fits the support exactly (Newton, Cyclops >= 4.0), moves by gradient hard-thresholding steps
#' and by swaps of weak and strong covariates ("splicing"), and accepts a move only when it improves the objective.
#' \code{algorithm = "chu2020"} runs iterative hard thresholding as in Chu et al. (2020) and their MendelIHT.jl.
#' \code{startingCoefficients} given to \code{fitCyclopsModel} are used as a warm start (for example an earlier
#' solution at a smaller \code{K}): their non-zero penalized coefficients form the starting support.
#'
#' \code{algorithm = "legacy"} runs version 1.0.3 of this package (Broken Adaptive Ridge with a cap of \code{K};
#' deprecated) so that earlier models can be reproduced.
#'
#' @param K              Number of non-zero penalized covariates: one value, or several to choose from by
#'                       cross-validation. The intercept and \code{exclude} do not count.
#' @param variance       Ridge variance \eqn{\sigma^2}: one value (\code{Inf} for plain L0), several values to choose
#'                       from by cross-validation, or \code{"tune"} for the default grid
#'                       (0.01, 0.03, 0.1, 0.3, 1, 3 and \code{Inf}, extended when the best value is at its edge).
#' @param exclude        Covariate ids or names that are not penalized and do not count towards \code{K}.
#' @param forceIntercept Logical: penalize the intercept and let it count towards \code{K}.
#' @param algorithm      \code{"default"}, \code{"chu2020"} or \code{"legacy"}. Sets the defaults of \code{control}.
#' @param cv             Cross-validation settings, from \code{\link{createIhtCvControl}}.
#' @param control        Algorithm switches, from \code{\link{createIhtControl}}; \code{NULL} entries take the
#'                       defaults of \code{algorithm}.
#' @param ...            Only with \code{algorithm = "legacy"}: the arguments of version 1.0.3 (\code{penalty},
#'                       \code{fitBestSubset}, \code{initialRidgeVariance}, \code{tolerance}, \code{maxIterations},
#'                       \code{threshold}, \code{delta}).
#'
#' @references
#' Chu BB, Keys KL, German CA, Zhou H, Zhou JJ, Sobel EM, Sinsheimer JS, Lange K (2020). Iterative hard thresholding
#' in genome-wide association studies: Generalized linear models, prior weights, and double sparsity.
#' GigaScience 9(6), giaa044.
#'
#' @examples
#' prior <- createIhtPrior(K = 10)
#'
#' @return
#' An IHT Cyclops prior object of class inheriting from \code{"cyclopsPrior"} for use with
#' \code{fitCyclopsModel}.
#'
#' @import Cyclops
#' @importFrom stats coef predict
#'
#' @export
createIhtPrior <- function(K,
                           variance = 1,
                           exclude = c(),
                           forceIntercept = FALSE,
                           algorithm = "default",
                           cv = createIhtCvControl(),
                           control = createIhtControl(),
                           ...) {
  algorithm <- match.arg(algorithm, c("default", "chu2020", "legacy"))
  legacyArgs <- list(...)
  if (algorithm == "legacy") {
    return(createLegacyPrior(K, exclude, forceIntercept, legacyArgs))
  }
  barOnly <- c("penalty", "fitBestSubset", "initialRidgeVariance", "tolerance", "maxIterations", "threshold", "delta")
  if (length(legacyArgs) > 0) {
    if (is.null(names(legacyArgs)) || any(names(legacyArgs) == "")) {
      stop("createIhtPrior: unnamed extra arguments")
    }
    unknown <- setdiff(names(legacyArgs), barOnly)
    if (length(unknown) > 0) {
      stop("Unknown argument(s) to createIhtPrior: ", paste(unknown, collapse = ", "))
    }
    stop("Argument(s) ", paste(names(legacyArgs), collapse = ", "), " belong to the method of version 1.0.3 (Broken ",
         "Adaptive Ridge with a cap of K), which version 2.0.0 replaced. Use algorithm = \"legacy\" to reproduce ",
         "earlier models, or see NEWS for the new arguments.")
  }
  if (is.character(variance) && length(variance) == 1 && tolower(variance) %in% c("bic", "aic")) {
    stop("variance = \"", variance, "\": the BIC/AIC penalty of version 1.0.3 is gone. Use algorithm = \"legacy\" ",
         "with penalty = \"", variance, "\" to reproduce earlier models, or see NEWS.")
  }
  K <- checkK(K)
  variance <- checkVariance(variance)
  if (!is.logical(forceIntercept) || length(forceIntercept) != 1) {
    stop("forceIntercept must be TRUE or FALSE")
  }
  if (!inherits(cv, "ihtCvControl")) stop("cv must come from createIhtCvControl()")
  if (!inherits(control, "ihtControl")) stop("control must come from createIhtControl()")
  control <- resolveControl(control, algorithm)

  structure(list(K = K,
                 variance = variance,
                 exclude = exclude,
                 forceIntercept = forceIntercept,
                 algorithm = algorithm,
                 cv = cv,
                 control = control,
                 keepControl = TRUE,
                 # reads everything from the prior object, so the hook survives a JSON round trip
                 # names cvFolds and returnCvPredictions, so Cyclops passes them (other hooks get its per-fold fallback)
                 fitHook = function(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                                    startingCoefficients, fixedCoefficients, cvFolds = NULL,
                                    returnCvPredictions = FALSE, ...) {
                   utils::getFromNamespace("ihtFitHook", "IterativeHardThresholding")(
                     cyclopsData, prior, control, weights, forceNewObject, returnEstimates, startingCoefficients,
                     fixedCoefficients, cvFolds = cvFolds, returnCvPredictions = returnCvPredictions)
                 }),
            class = c("cyclopsPrior", "ihtPrior"))
}

#' @title Algorithm switches for IHT
#'
#' @description
#' Every switch defaults to \code{NULL}, which takes the value of the \code{algorithm} given to
#' \code{\link{createIhtPrior}}. The switches exist so that each part of the default algorithm can be turned off
#' and compared with the algorithm of Chu et al. (2020).
#'
#' @param exactRefits      Fit each support exactly (Newton). \code{FALSE} runs plain IHT iterations (gradient step,
#'                         top-\code{K} projection, backtracking) as in Chu et al. (2020); default \code{TRUE},
#'                         \code{"chu2020"}: \code{FALSE}.
#' @param thresholdSteps   With exact refits: try a gradient hard-thresholding step before splicing.
#' @param splicing         With exact refits: swap the weakest covariates in the support for the strongest outside.
#' @param start            With exact refits: \code{"gradient"} starts from the \code{K} covariates with the largest
#'                         \eqn{|g_j|} in the model with only the unpenalized covariates (as Chu et al. 2020 do); \code{"warmup"} runs plain
#'                         IHT iterations until the support is stable for 10 iterations.
#' @param strictAcceptance Accept a move only if it improves the objective by more than \code{tolerance}.
#'                         \code{FALSE} accepts the last backtracking trial, as Chu et al. (2020) do.
#' @param maxHalve         Maximum step halvings per iteration (default 5; \code{"chu2020"}: 3).
#' @param maxSwap          Largest swap size; default \code{ceiling(K / 4)}.
#' @param tolerance        With exact refits: smallest objective gain that counts as an improvement (default 1e-3).
#'                         Without: the convergence tolerance on the scaled change in coefficients (1e-4).
#' @param minIterations    Without exact refits: minimum iterations before testing convergence (5).
#' @param maxIterations    Maximum iterations (200).
#' @param threads          Threads for the refits and evaluations (default 1: the refits on a small support are faster
#'                         on one thread than on several).
#'
#' @return An object of class \code{"ihtControl"}.
#'
#' @export
createIhtControl <- function(exactRefits = NULL,
                             thresholdSteps = NULL,
                             splicing = NULL,
                             start = NULL,
                             strictAcceptance = NULL,
                             maxHalve = NULL,
                             maxSwap = NULL,
                             tolerance = NULL,
                             minIterations = NULL,
                             maxIterations = NULL,
                             threads = NULL) {
  ctl <- list(exactRefits = exactRefits, thresholdSteps = thresholdSteps, splicing = splicing,
              start = start, strictAcceptance = strictAcceptance,
              maxHalve = maxHalve, maxSwap = maxSwap, tolerance = tolerance, minIterations = minIterations,
              maxIterations = maxIterations, threads = threads)
  structure(ctl, class = "ihtControl")
}

#' @title Cross-validation settings for IHT
#'
#' @description
#' Cross-validation runs when \code{K} or \code{variance} has more than one candidate, or when out-of-fold
#' predictions are requested. Each fold runs the whole selection on its training rows, so the out-of-fold
#' predictions do not use the held-out rows.
#'
#' @param folds              Number of folds when \code{cvFolds} is not given.
#' @param cvFolds            Optional \code{data.frame(rowId, fold)} assigning rows to folds.
#' @param seed               Seed for the random fold assignment (\code{NULL}: not set).
#' @param oneStandardError   Choose the smallest \code{K} whose mean held-out log-likelihood is within one standard
#'                           error of the best, instead of the best.
#' @param returnCvPredictions Return the out-of-fold predictions at the chosen \code{K} and variance.
#'
#' @return An object of class \code{"ihtCvControl"}.
#'
#' @export
createIhtCvControl <- function(folds = 3,
                               cvFolds = NULL,
                               seed = NULL,
                               oneStandardError = FALSE,
                               returnCvPredictions = FALSE) {
  if (!is.null(cvFolds) && !all(c("rowId", "fold") %in% names(cvFolds))) {
    stop("cvFolds must be a data.frame with columns rowId and fold")
  }
  structure(list(folds = folds, cvFolds = cvFolds, seed = seed, oneStandardError = oneStandardError,
                 returnCvPredictions = returnCvPredictions),
            class = "ihtCvControl")
}

# ---- package-private

presets <- list(
  default = list(exactRefits = TRUE, thresholdSteps = TRUE, splicing = TRUE,
                 start = "gradient", strictAcceptance = TRUE, maxHalve = 5, maxSwap = NA,
                 tolerance = 1e-3, minIterations = 0, maxIterations = 200, threads = 1),
  chu2020 = list(exactRefits = FALSE, thresholdSteps = TRUE, splicing = FALSE,
                 start = "gradient", strictAcceptance = FALSE, maxHalve = 3, maxSwap = NA,
                 tolerance = 1e-4, minIterations = 5, maxIterations = 200, threads = 1)
)

resolveControl <- function(control, algorithm) {
  ctl <- presets[[algorithm]]
  for (k in names(control)) {
    if (!is.null(control[[k]])) ctl[[k]] <- control[[k]]
  }
  if (ctl$exactRefits) {
    ctl$start <- match.arg(ctl$start, c("gradient", "warmup"))
  } else if (isTRUE(ctl$splicing) && !identical(control$splicing, NULL)) {
    stop("splicing needs exactRefits = TRUE")
  } else {
    ctl$splicing <- FALSE
  }
  structure(ctl, class = "ihtControl")
}

checkK <- function(K) {
  if (!is.numeric(K) || length(K) < 1 || any(is.na(K)) || any(K < 1) || any(K != round(K))) {
    stop("K must be one or more whole numbers >= 1")
  }
  sort(unique(as.integer(K)))
}

checkVariance <- function(variance) {
  if (is.character(variance)) {
    if (length(variance) == 1 && variance == "tune") return("tune")
    v <- suppressWarnings(as.numeric(ifelse(tolower(variance) == "inf", "Inf", variance)))
    if (any(is.na(v))) stop("variance must be numeric (Inf allowed) or \"tune\"")
    variance <- v
  }
  if (!is.numeric(variance) || length(variance) < 1 || any(is.na(variance)) || any(variance <= 0)) {
    stop("variance must be positive (Inf allowed), or \"tune\"")
  }
  sort(unique(variance))
}
