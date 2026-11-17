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
#' \code{createIhtPrior} creates a prior for \code{\link[Cyclops]{fitCyclopsModel}} that keeps at
#' most \code{K} non-zero penalized coefficients. It maximizes the log-likelihood minus a small
#' ridge term,
#' \deqn{\ell(\beta) - \sum_{j \in P} \beta_j^2 / (2 \sigma^2) \quad \textrm{subject to at most } K
#'   \textrm{ non-zero } \beta_j, j \in P,}
#' where \eqn{P} are the penalized covariates (all but the intercept and \code{exclude}) and
#' \eqn{\sigma^2} is \code{variance}. \code{variance = Inf} is plain L0 (best subset of size
#' \code{K}).
#'
#' The model is the one of Chu et al. (2020); the default algorithm only optimizes it differently.
#' It fits each support exactly (Newton, Cyclops >= 4.0), as in gradient hard thresholding pursuit
#' (Yuan et al. 2018) and Newton hard-thresholding pursuit (Zhou et al. 2021). It moves by gradient
#' hard-thresholding steps and by swaps of weak and strong covariates ("splicing", Zhu et al. 2020;
#' see also Beck and Eldar 2013, Hazimeh and Mazumder 2020), and accepts a move only when it
#' improves the objective.
#' \code{algorithm = "chu2020"} runs iterative hard thresholding as in Chu et al. (2020) and their
#' MendelIHT.jl.
#'
#' \code{startingCoefficients} given to \code{fitCyclopsModel} are used as a warm start (for example
#' an earlier solution at a smaller \code{K}): their non-zero penalized coefficients form the
#' starting support. \code{K} and \code{variance} are chosen by the caller, for example by
#' cross-validation in PatientLevelPrediction. \code{fitCyclopsModel(..., returnCvPredictions =
#' TRUE)} returns out-of-fold predictions at the given \code{K} and \code{variance}.
#'
#' \code{algorithm = "legacy"} runs version 1.0.3 of this package (Broken Adaptive Ridge with a
#' projection onto the \code{K} largest coefficients after each fit; deprecated) so that earlier
#' models can be reproduced.
#'
#' @param K              Number of non-zero penalized covariates. The intercept and \code{exclude}
#'                       do not count.
#' @param variance       Ridge variance \eqn{\sigma^2} (\code{Inf} for plain L0).
#' @param exclude        Covariate ids or names that are not penalized and do not count towards
#'                       \code{K}.
#' @param forceIntercept Logical: penalize the intercept and let it count towards \code{K}.
#' @param algorithm      \code{"default"}, \code{"chu2020"} or \code{"legacy"}. Sets the defaults of
#'                       \code{control}.
#' @param control        Algorithm switches, from \code{\link{createIhtControl}}; \code{NULL}
#'                       entries take the defaults of \code{algorithm}.
#' @param ...            Only with \code{algorithm = "legacy"}: the arguments of version 1.0.3
#'                       (\code{penalty}, \code{fitBestSubset}, \code{initialRidgeVariance},
#'                       \code{tolerance}, \code{maxIterations}, \code{threshold}, \code{delta}).
#'
#' @references
#' Chu BB, Keys KL, German CA, Zhou H, Zhou JJ, Sobel EM, Sinsheimer JS, Lange K (2020). Iterative
#' hard thresholding in genome-wide association studies: Generalized linear models, prior weights,
#' and double sparsity. GigaScience 9(6), giaa044. \doi{10.1093/gigascience/giaa044}
#'
#' Bahmani S, Raj B, Boufounos PT (2013). Greedy sparsity-constrained optimization. Journal of
#' Machine Learning Research 14(25), 807-841.
#'
#' Beck A, Eldar YC (2013). Sparsity constrained nonlinear optimization: optimality conditions and
#' algorithms. SIAM Journal on Optimization 23(3), 1480-1509. \doi{10.1137/120869778}
#'
#' Yuan XT, Li P, Zhang T (2018). Gradient hard thresholding pursuit. Journal of Machine Learning
#' Research 18(166), 1-43.
#'
#' Zhou S, Xiu N, Qi HD (2021). Global and quadratic convergence of Newton hard-thresholding
#' pursuit. Journal of Machine Learning Research 22(12), 1-45.
#'
#' Zhu J, Wen C, Zhu J, Zhang H, Wang X (2020). A polynomial algorithm for best-subset selection
#' problem. Proceedings of the National Academy of Sciences 117(52), 33117-33123.
#' \doi{10.1073/pnas.2014241117}
#'
#' Hazimeh H, Mazumder R (2020). Fast best subset selection: coordinate descent and local
#' combinatorial optimization algorithms. Operations Research 68(5), 1517-1537.
#' \doi{10.1287/opre.2019.1919}
#'
#' Dedieu A, Hazimeh H, Mazumder R (2021). Learning sparse classifiers: continuous and mixed integer
#' optimization perspectives. Journal of Machine Learning Research 22(135), 1-47.
#'
#' @examples
#' prior <- createIhtPrior(K = 10)
#'
#' @return
#' An IHT Cyclops prior object of class inheriting from \code{"cyclopsPrior"} for use with
#' \code{fitCyclopsModel}.
#'
#' @import Cyclops
#' @importFrom stats coef
#'
#' @export
createIhtPrior <- function(K,
                           variance = 1,
                           exclude = c(),
                           forceIntercept = FALSE,
                           algorithm = "default",
                           control = createIhtControl(),
                           ...) {
  algorithm <- match.arg(algorithm, c("default", "chu2020", "legacy"))
  legacyArgs <- list(...)
  if (algorithm == "legacy") {
    return(createLegacyPrior(K, exclude, forceIntercept, legacyArgs))
  }
  barOnly <- c("penalty", "fitBestSubset", "initialRidgeVariance", "tolerance", "maxIterations",
               "threshold", "delta")
  if (length(legacyArgs) > 0) {
    if (is.null(names(legacyArgs)) || any(names(legacyArgs) == "")) {
      stop("createIhtPrior: unnamed extra arguments")
    }
    unknown <- setdiff(names(legacyArgs), barOnly)
    if (length(unknown) > 0) {
      stop("Unknown argument(s) to createIhtPrior: ", paste(unknown, collapse = ", "))
    }
    stop("Argument(s) ", paste(names(legacyArgs), collapse = ", "), " belong to the method of ",
         "version 1.0.3 (Broken Adaptive Ridge with a projection onto the K largest ",
         "coefficients), which version 2.0.0 replaced. Use algorithm = \"legacy\" to reproduce ",
         "earlier models, or see NEWS for the new arguments.")
  }
  if (is.character(variance) && length(variance) == 1 && tolower(variance) %in% c("bic", "aic")) {
    stop("variance = \"", variance, "\": the BIC/AIC penalty of version 1.0.3 is gone. Use ",
         "algorithm = \"legacy\" with penalty = \"", variance, "\" to reproduce earlier models, ",
         "or see NEWS.")
  }
  K <- checkK(K)
  variance <- checkVariance(variance)
  if (!is.logical(forceIntercept) || length(forceIntercept) != 1) {
    stop("forceIntercept must be TRUE or FALSE")
  }
  if (!inherits(control, "ihtControl")) stop("control must come from createIhtControl()")
  control <- resolveControl(control, algorithm)

  structure(list(K = K,
                 variance = variance,
                 exclude = exclude,
                 forceIntercept = forceIntercept,
                 algorithm = algorithm,
                 control = control,
                 keepControl = TRUE,
                 # reads everything from the prior object, so it survives a JSON round trip
                 fitHook = function(cyclopsData, prior, control, weights, forceNewObject,
                                    returnEstimates, startingCoefficients, fixedCoefficients,
                                    ...) {
                   hook <- utils::getFromNamespace("ihtFitHook", "IterativeHardThresholding")
                   hook(cyclopsData, prior, control, weights, forceNewObject, returnEstimates,
                        startingCoefficients, fixedCoefficients)
                 }),
            class = c("cyclopsPrior", "ihtPrior"))
}

#' @title Algorithm switches for IHT
#'
#' @description
#' Every switch defaults to \code{NULL}, which takes the value of the \code{algorithm} given to
#' \code{\link{createIhtPrior}}. The switches exist so that each part of the default algorithm can
#' be turned off and compared with the algorithm of Chu et al. (2020).
#'
#' @param exactRefits      Fit each support exactly (Newton). \code{FALSE} runs plain IHT
#'                         iterations (gradient step, top-\code{K} projection, backtracking) as in
#'                         Chu et al. (2020); default \code{TRUE}, \code{"chu2020"}: \code{FALSE}.
#' @param thresholdSteps   With exact refits: try a gradient hard-thresholding step before
#'                         splicing.
#' @param splicing         With exact refits: swap the weakest covariates in the support for the
#'                         strongest outside.
#' @param start            With exact refits: \code{"gradient"} starts from the \code{K} covariates
#'                         with the largest \eqn{|g_j|} in the model with only the unpenalized
#'                         covariates (as Chu et al. 2020 do); \code{"warmup"} runs plain IHT
#'                         iterations until the support is stable for 10 iterations.
#' @param strictAcceptance Accept a move only if it improves the objective by more than
#'                         \code{tolerance}. \code{FALSE} accepts the last backtracking trial, as
#'                         Chu et al. (2020) do.
#' @param maxHalve         Maximum step halvings per iteration (default 5; \code{"chu2020"}: 3).
#' @param maxSwap          Largest swap size; default \code{ceiling(K / 4)}.
#' @param tolerance        With exact refits: smallest objective gain that counts as an improvement
#'                         (default 1e-3). Without: the convergence tolerance on the scaled change
#'                         in coefficients (1e-4).
#' @param minIterations    Without exact refits: minimum iterations before testing convergence (5).
#' @param maxIterations    Maximum iterations (200).
#' @param threads          Threads for the refits and evaluations (default 1: the refits on a small
#'                         support are faster on one thread than on several).
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
              start = start, strictAcceptance = strictAcceptance, maxHalve = maxHalve,
              maxSwap = maxSwap, tolerance = tolerance, minIterations = minIterations,
              maxIterations = maxIterations, threads = threads)
  structure(ctl, class = "ihtControl")
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
  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K < 1 || K != round(K)) {
    stop("K must be one whole number >= 1")
  }
  as.integer(K)
}

# "Inf" as text is accepted, since a JSON round trip turns Inf into a string
checkVariance <- function(variance) {
  if (is.character(variance) && length(variance) == 1 && tolower(variance) == "inf") {
    variance <- Inf
  }
  if (!is.numeric(variance) || length(variance) != 1 || is.na(variance) || variance <= 0) {
    stop("variance must be one positive number (Inf allowed)")
  }
  variance
}
