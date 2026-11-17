# @file Engine.R
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

# The two fitting engines, for one K and one variance:
#   exactEngine:    exact fits on the support, gradient hard-thresholding steps and splicing
#                   (the default); covariates are ranked by |b| (in the support) and |gradient|
#                   (outside), as in Chu et al. (2020)
#   gradientEngine: plain IHT iterations as in Chu et al. (2020) and MendelIHT.jl
#                   (algorithm = "chu2020")
# Both maximize F(b) = loglik(b) - sum(b[pen]^2) / (2 * variance) with at most K non-zero b[pen].
# Cyclops does all the fitting and evaluation; the engines only choose the supports.
# A solution is list(b = coefficients, support = penalized indices with b != 0, value = F(b)).

# A problem: one data object, weights, penalized set and variance.
newProblem <- function(cyclopsData, unpenalized, variance, control, weights, forceIntercept,
                       intercept) {
  p <- Cyclops::getNumberOfCovariates(cyclopsData)
  ids <- Cyclops::getCovariateIds(cyclopsData)
  excludeIds <- as.numeric(ids[setdiff(unpenalized, intercept)])
  prior <- if (is.finite(variance)) {
    Cyclops::createPrior("normal", variance = variance,
                         exclude = if (length(excludeIds)) excludeIds else c(),
                         forceIntercept = forceIntercept)
  } else {
    Cyclops::createPrior("none")
  }
  pb <- new.env(parent = emptyenv())
  pb$data <- cyclopsData
  pb$p <- p
  pb$ids <- ids
  pb$pen <- setdiff(seq_len(p), unpenalized)
  pb$unpenalized <- unpenalized
  pb$lam <- if (is.finite(variance)) 1 / variance else 0
  pb$prior <- prior
  pb$control <- control
  pb$weights <- weights
  pb$threads <- control$threads
  pb$refits <- 0
  pb$failedRefits <- 0
  pb$lastFailure <- NA_character_
  pb
}

ridge <- function(pb, b) pb$lam * sum(b[pb$pen]^2) / 2

# The objective F, its gradient and the Fisher information, from Cyclops::evaluateLikelihood()
objective <- function(pb, b) {
  ll <- Cyclops::evaluateLikelihood(pb$data, b, weights = pb$weights, gradient = FALSE,
                                    threads = pb$threads)
  ll$logLikelihood - ridge(pb, b)
}

objectiveGradient <- function(pb, b) {
  g <- Cyclops::evaluateLikelihood(pb$data, b, weights = pb$weights, threads = pb$threads)
  g <- as.numeric(g$gradient)
  g[pb$pen] <- g[pb$pen] - pb$lam * b[pb$pen]
  g
}

# Fisher information of the log-likelihood on columns `cols` (indices); the caller adds the ridge
fisher <- function(pb, b, cols) {
  m <- Cyclops::evaluateLikelihood(pb$data, b, weights = pb$weights, gradient = FALSE,
                                   information = pb$ids[cols], threads = pb$threads)
  unname(m$information)
}

# The K largest entries of v among the penalized covariates (the hard-thresholding projection)
topK <- function(pb, v, K) pb$pen[order(-abs(v[pb$pen]))[seq_len(K)]]

# The k largest entries of v among penalized covariates outside the support
topOutside <- function(pb, v, support, k) {
  o <- setdiff(pb$pen, support)
  o[order(-v[o])[seq_len(min(k, length(o)))]]
}

failRefit <- function(pb, reason) {
  pb$failedRefits <- pb$failedRefits + 1
  pb$lastFailure <- reason
  NULL
}

# Exact maximizer of F with the penalized coefficients outside the support fixed at 0,
# warm-started from `start`. Returns NULL when Cyclops does not report SUCCESS (Cyclops reports a
# separation as ILLCONDITIONED, with fit$divergentCoefficients).
refit <- function(pb, support, start) {
  keep <- seq_len(pb$p) %in% c(pb$unpenalized, support)
  start[!keep] <- 0
  pb$refits <- pb$refits + 1
  fit <- Cyclops::fitCyclopsModel(pb$data, prior = pb$prior, control = pb$control,
                                  weights = pb$weights, startingCoefficients = start,
                                  fixedCoefficients = !keep, warnings = FALSE)
  if (fit$return_flag != "SUCCESS") return(failRefit(pb, fit$return_flag))
  b <- as.numeric(coef(fit))
  value <- fit$log_likelihood - ridge(pb, b)
  # Backstop for a separation that Cyclops does not flag: without a ridge term the fit can stop
  # at huge coefficients. Treat that as a failed refit.
  if (pb$lam == 0 && (!is.finite(value) || max(abs(b)) > separationBound)) {
    return(failRefit(pb, "SEPARATION"))
  }
  list(b = b, support = pb$pen[b[pb$pen] != 0], value = value)
}

# |coefficient| above which an unpenalized fit is taken to be diverging
# (log odds ratio of 20 = odds ratio 5e8)
separationBound <- 20

# The fit with only the unpenalized covariates
baseFit <- function(pb) {
  base <- refit(pb, integer(0), rep(0, pb$p))
  if (is.null(base)) {
    stop("The fit with only the unpenalized covariates failed: ", pb$lastFailure)
  }
  base
}

# Start support from the fit with only the unpenalized covariates
startSupport <- function(pb, base, K, ctl) {
  if (ctl$start == "warmup") {
    w <- gradientEngine(pb, K, ctl, from = base, stableStop = 10, maxIterations = 500)
    return(w$support)
  }
  topOutside(pb, abs(objectiveGradient(pb, base$b)), integer(0), K)
}

# A finite fit on the start support can fail without a ridge term (separation). Build it
# greedily instead: add covariates in the order of |gradient| and skip any that make it fail.
greedyStart <- function(pb, cur, candidates, K) {
  support <- integer(0)
  for (j in candidates) {
    new <- refit(pb, c(support, j), cur$b)
    if (!is.null(new)) {
      cur <- new
      support <- c(support, j)
    }
    if (length(support) == K) break
  }
  cur
}

# warm: an earlier solution list(b, support) (another K or variance); the start support keeps its
# covariates (the strongest of them if it has more than K) and adds the strongest outside ones.
exactEngine <- function(pb, K, ctl, warm = NULL) {
  base <- baseFit(pb)
  if (is.null(warm)) {
    start <- startSupport(pb, base, K, ctl)
    b0 <- base$b
  } else {
    kept <- warm$support[order(-abs(warm$b[warm$support]))]
    start <- kept[seq_len(min(K, length(kept)))]
    if (length(start) < K) {
      g <- objectiveGradient(pb, warm$b)
      start <- c(start, topOutside(pb, abs(g), start, K - length(start)))
    }
    b0 <- warm$b
  }
  cur <- refit(pb, start, b0)
  if (is.null(cur)) {
    candidates <- topOutside(pb, abs(objectiveGradient(pb, base$b)), integer(0), 10 * K)
    cur <- greedyStart(pb, base, unique(c(start, candidates)), K)
  }
  steps <- 0
  splices <- 0
  converged <- FALSE
  it <- 0
  better <- function(new) !is.null(new) && new$value > cur$value + ctl$tolerance
  for (it in seq_len(ctl$maxIterations)) {
    g <- objectiveGradient(pb, cur$b)
    strong <- topOutside(pb, abs(g), cur$support, K) # the strongest outside, by |g|
    moved <- FALSE
    if (ctl$thresholdSteps) {
      gO <- g[strong]
      curvature <- as.numeric(t(gO) %*% fisher(pb, cur$b, strong) %*% gO) + pb$lam * sum(gO^2)
      mu <- sum(gO^2) / curvature
      if (!is.finite(mu)) mu <- 0
      tried <- NULL # the last trial; a halving often gives the same support again
      for (h in 0:ctl$maxHalve) {
        q <- cur$b + mu * g
        trial <- sort(topK(pb, q, K))
        if (setequal(trial, cur$support)) break
        if (!identical(trial, tried$support)) {
          tried <- list(support = trial, fit = refit(pb, trial, q))
        }
        new <- tried$fit
        lastTrial <- !ctl$strictAcceptance && h == ctl$maxHalve && !is.null(new)
        if (better(new) || lastTrial) {
          cur <- new
          steps <- steps + 1
          moved <- TRUE
          break
        }
        mu <- mu / 2
      }
    }
    if (moved) next
    if (ctl$splicing && length(cur$support) > 0) {
      # swap the weakest insiders by |b| for the strongest outsiders by |g|
      weak <- cur$support[order(abs(cur$b[cur$support]))]
      if (length(cur$support) < K) {
        # a support short of K (greedy start): fill it up once
        new <- refit(pb, c(cur$support, strong[seq_len(K - length(cur$support))]), cur$b)
        if (better(new)) {
          cur <- new
          splices <- splices + 1
          moved <- TRUE
        }
      } else {
        s <- if (is.na(ctl$maxSwap)) ceiling(length(cur$support) / 4) else ctl$maxSwap
        s <- min(s, length(strong), length(weak))
        while (s >= 1) {
          trial <- c(setdiff(cur$support, weak[seq_len(s)]), strong[seq_len(s)])
          new <- refit(pb, trial, cur$b)
          if (better(new)) {
            cur <- new
            splices <- splices + 1
            moved <- TRUE
            break
          }
          if (s == 1) break
          s <- ceiling(s / 2)
        }
      }
    }
    if (!moved) {
      converged <- TRUE
      break
    }
  }
  list(b = cur$b, support = cur$support, value = cur$value, iterations = it,
       converged = converged, steps = steps, splices = splices)
}

# Plain IHT as in MendelIHT.jl 1.4 (fit_iht!):
# - start from the unpenalized MLE, with the support of the K largest |gradient|;
# - step size ||g_S||^2 / g_S' J_S g_S on the current support S and the unpenalized covariates;
# - full gradient step, then keep the K largest |b| (penalized);
# - backtrack (halve) while the objective is below the previous iterate's, at most maxHalve
#   times, then accept;
# - stop when the scaled change max|b - b0| / (max|b0| + 1) is below tolerance after at least
#   minIterations, or at maxIterations;
# - return the iterate with the best objective.
# With stableStop, stop instead when the support has not changed for that many iterations
# (warm-up use).
gradientEngine <- function(pb, K, ctl, from = NULL, stableStop = NULL,
                           maxIterations = ctl$maxIterations) {
  cur <- if (is.null(from)) baseFit(pb) else from
  pen <- pb$pen
  b <- cur$b
  g <- objectiveGradient(pb, b)
  S <- topK(pb, g, K)
  gStep <- g
  gStep[setdiff(pen, S)] <- 0 # MendelIHT projects the first score to the starting support
  project <- function(q) {
    q[setdiff(pen, topK(pb, q, K))] <- 0
    q
  }
  prevF <- -Inf
  best <- list(value = -Inf, b = b)
  converged <- FALSE
  stable <- 0
  it <- 0
  Fb <- objective(pb, b)
  for (it in seq_len(maxIterations)) {
    if (Fb > best$value) best <- list(value = Fb, b = b)
    # as fit_iht!: stop at the start of iteration max_iter
    if (is.null(stableStop) && it >= maxIterations) break
    b0 <- b
    v <- c(pb$unpenalized, S)
    gv <- gStep[v]
    den <- as.numeric(t(gv) %*% fisher(pb, b, v) %*% gv) + pb$lam * sum(gStep[S]^2)
    eta <- sum(gv^2) / den
    if (!is.finite(eta)) eta <- 1e-8
    bn <- project(b0 + eta * gStep)
    Fn <- objective(pb, bn)
    h <- 0
    while (prevF > Fn && h < ctl$maxHalve) {
      eta <- eta / 2
      bn <- project(b0 + eta * gStep)
      Fn <- objective(pb, bn)
      h <- h + 1
    }
    Sn <- pen[bn[pen] != 0]
    stable <- if (setequal(Sn, S)) stable + 1 else 0
    b <- bn
    S <- Sn
    Fb <- Fn
    prevF <- Fn
    gStep <- objectiveGradient(pb, b)
    if (!is.null(stableStop)) {
      if (stable >= stableStop) {
        converged <- TRUE
        break
      }
    } else {
      scaled <- max(abs(b - b0)) / (max(abs(b0)) + 1)
      if (it >= ctl$minIterations && scaled < ctl$tolerance) {
        converged <- TRUE
        break
      }
    }
  }
  if (Fb > best$value) best <- list(value = Fb, b = b)
  list(b = best$b, support = pen[best$b[pen] != 0], value = best$value, iterations = it,
       converged = converged, steps = it, splices = 0)
}
