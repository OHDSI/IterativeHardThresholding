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

# The two fitting engines, for one K and one variance on one set of weights:
#   exactEngine:    exact fits on the support + gradient hard-thresholding steps + splicing (the default); covariates
#                   are ranked by |b| (in the support) and |gradient| (outside), as in Chu et al. (2020)
#   gradientEngine: plain IHT iterations as in Chu et al. (2020) / MendelIHT.jl (algorithm = "chu2020")
# Both maximize F(b) = loglik(b) - sum(b[pen]^2) / (2 * variance) with at most K non-zero b[pen].

# A problem: one data object, weights, penalized set and variance. All Cyclops access goes through here.
newProblem <- function(cyclopsData, unpenalized, variance, control, weights, forceIntercept) {
  p <- Cyclops::getNumberOfCovariates(cyclopsData)
  ids <- Cyclops::getCovariateIds(cyclopsData)
  pen <- setdiff(seq_len(p), unpenalized)
  lam <- if (is.finite(variance)) 1 / variance else 0
  excludeIds <- as.numeric(ids[setdiff(unpenalized, interceptIndex(cyclopsData))])
  prior <- if (is.finite(variance)) {
    Cyclops::createPrior("normal", variance = variance, exclude = if (length(excludeIds)) excludeIds else c(),
                         forceIntercept = forceIntercept)
  } else {
    Cyclops::createPrior("none")
  }
  pb <- new.env(parent = emptyenv())
  pb$data <- cyclopsData; pb$p <- p; pb$ids <- ids; pb$pen <- pen; pb$unpenalized <- unpenalized
  pb$lam <- lam; pb$variance <- variance; pb$prior <- prior; pb$control <- control; pb$weights <- weights
  pb$refits <- 0; pb$failedRefits <- 0; pb$failures <- character(0)
  pb$threads <- if (is.null(control$threads) || control$threads < 1) 1 else control$threads
  pb
}

# All evaluations go through Cyclops::evaluateLikelihood() (Cyclops >= 4.0.1), with the problem's weights.
evaluate <- function(pb, b, gradient = FALSE, information = NULL, informationDiagonal = NULL, weights = pb$weights) {
  Cyclops::evaluateLikelihood(pb$data, b, weights = weights, gradient = gradient, information = information,
                              informationDiagonal = informationDiagonal, threads = pb$threads)
}

# Unpenalized log-likelihood
logLik_ <- function(pb, b) evaluate(pb, b)$logLikelihood

objective <- function(pb, b) logLik_(pb, b) - pb$lam * sum(b[pb$pen]^2) / 2

# Gradient of F
gradient <- function(pb, b) {
  g <- as.numeric(evaluate(pb, b, gradient = TRUE)$gradient)
  g[pb$pen] <- g[pb$pen] - pb$lam * b[pb$pen]
  g
}

# Fisher information of the log-likelihood on columns `cols` (indices); the caller adds the ridge term
fisher <- function(pb, b, cols) {
  m <- evaluate(pb, b, information = as.numeric(pb$ids[cols]))$information
  matrix(as.numeric(m), length(cols))
}

# Exact maximizer of F with the penalized coefficients outside T fixed at 0, warm-started from `start`.
# Returns NULL when Cyclops does not report SUCCESS (e.g. separation without a ridge term).
refit <- function(pb, T, start) {
  keep <- seq_len(pb$p) %in% c(pb$unpenalized, T)
  start[!keep] <- 0
  pb$refits <- pb$refits + 1
  fit <- suppressWarnings(Cyclops::fitCyclopsModel(pb$data, prior = pb$prior, control = pb$control,
                                                   weights = pb$weights, startingCoefficients = start,
                                                   fixedCoefficients = !keep, warnings = FALSE))
  if (fit$return_flag != "SUCCESS") {
    pb$failedRefits <- pb$failedRefits + 1
    pb$failures <- c(pb$failures, fit$return_flag)
    return(NULL)
  }
  b <- as.numeric(coef(fit))
  F <- objective(pb, b)
  # Without a ridge term, a support that separates the data has no finite maximizer, yet the fit can stop with
  # SUCCESS at huge coefficients. Treat that as a failed refit, as for a reported separation.
  if (pb$lam == 0 && (!is.finite(F) || max(abs(b)) > separationBound)) {
    pb$failedRefits <- pb$failedRefits + 1
    pb$failures <- c(pb$failures, "SEPARATION")
    return(NULL)
  }
  list(b = b, T = pb$pen[b[pb$pen] != 0], F = F)
}

# |coefficient| above which an unpenalized fit is taken to be diverging (log odds ratio of 20 = odds ratio 5e8)
separationBound <- 20

# The k largest entries of v among penalized covariates outside T
topOutside <- function(pb, v, T, k) {
  o <- setdiff(pb$pen, T)
  o[order(-v[o])[seq_len(min(k, length(o)))]]
}

# Start support from the fit with only the unpenalized covariates
startSupport <- function(pb, cur, K, ctl) {
  if (ctl$start == "warmup") {
    w <- gradientEngine(pb, K, ctl, from = cur, stableStop = 10, maxIterations = 500)
    return(w$T)
  }
  g <- gradient(pb, cur$b)
  topOutside(pb, abs(g), integer(0), K)
}

# A finite fit on the start support can fail without a ridge term (separation). Build it greedily instead:
# add covariates in the order of |gradient| and skip any that make the fit fail.
greedyStart <- function(pb, cur, order, K) {
  T <- integer(0)
  for (j in order) {
    new <- refit(pb, c(T, j), cur$b)
    if (!is.null(new)) { cur <- new; T <- c(T, j) }
    if (length(T) == K) break
  }
  cur
}

# warm: an earlier solution list(b, T) (another K or variance); the start support keeps its covariates (the strongest
# of them if it has more than K) and adds the strongest outside covariates.
exactEngine <- function(pb, K, ctl, warm = NULL) {
  cur <- refit(pb, integer(0), rep(0, pb$p))
  if (is.null(cur)) stop("The fit with only the unpenalized covariates failed: ", utils::tail(pb$failures, 1))
  if (is.null(warm)) {
    T0 <- startSupport(pb, cur, K, ctl); b0 <- cur$b
  } else {
    T0 <- warm$T[order(-abs(warm$b[warm$T]))][seq_len(min(K, length(warm$T)))]
    if (length(T0) < K) {
      g <- gradient(pb, warm$b)
      T0 <- c(T0, topOutside(pb, abs(g), T0, K - length(T0)))
    }
    b0 <- warm$b
  }
  first <- refit(pb, T0, b0)
  if (is.null(first)) {
    g <- gradient(pb, cur$b)
    candidates <- topOutside(pb, abs(g), integer(0), 10 * K)
    first <- greedyStart(pb, cur, unique(c(T0, candidates)), K)
  }
  cur <- first
  maxSwap <- if (is.na(ctl$maxSwap)) NULL else ctl$maxSwap
  steps <- 0; splices <- 0; converged <- FALSE; it <- 0
  better <- function(new) !is.null(new) && new$F > cur$F + ctl$tolerance
  for (it in seq_len(ctl$maxIterations)) {
    g <- gradient(pb, cur$b)
    moved <- FALSE
    if (ctl$thresholdSteps) {
      O <- topOutside(pb, abs(g), cur$T, K)
      gO <- g[O]
      mu <- sum(gO^2) / (as.numeric(t(gO) %*% fisher(pb, cur$b, O) %*% gO) + pb$lam * sum(gO^2))
      if (!is.finite(mu)) mu <- 0
      for (h in 0:ctl$maxHalve) {
        q <- cur$b + mu * g
        Tn <- pb$pen[order(-abs(q[pb$pen]))[seq_len(K)]]
        if (setequal(Tn, cur$T)) break
        new <- refit(pb, Tn, ifelse(seq_len(pb$p) %in% c(pb$unpenalized, Tn), q, 0))
        if (better(new) || (!ctl$strictAcceptance && h == ctl$maxHalve && !is.null(new))) {
          cur <- new; steps <- steps + 1; moved <- TRUE
          break
        }
        mu <- mu / 2
      }
    }
    if (moved) next
    if (ctl$splicing && length(cur$T) > 0) {
      # weakest insiders by |b|, strongest outsiders by |g| (as the thresholding step ranks them)
      weak <- cur$T[order(abs(cur$b[cur$T]))]
      strong <- topOutside(pb, abs(g), cur$T, K)
      s <- if (is.null(maxSwap)) ceiling(max(1, length(cur$T)) / 4) else maxSwap
      repeat {
        nOut <- min(s, length(weak))
        if (length(cur$T) < K) {
          Tn <- c(cur$T, strong[seq_len(K - length(cur$T))])
        } else {
          Tn <- c(setdiff(cur$T, weak[seq_len(nOut)]), strong[seq_len(s)])
        }
        new <- refit(pb, Tn, cur$b)
        if (better(new)) { cur <- new; splices <- splices + 1; moved <- TRUE; break }
        if (s == 1) break
        s <- ceiling(s / 2)
      }
    }
    if (!moved) { converged <- TRUE; break }
  }
  list(b = cur$b, T = cur$T, F = cur$F, iterations = it, converged = converged, steps = steps, splices = splices)
}

# Plain IHT as in MendelIHT.jl 1.4 (fit_iht!): start from the unpenalized MLE with the support of the K largest
# |gradient|; step size ||g_S||^2 / g_S' J_S g_S on the current support S and the unpenalized covariates; full
# gradient step, keep the K largest |b| (penalized); backtrack (halve) while the objective is below the previous
# iterate's, at most maxHalve times, then accept; stop when the scaled change max|b - b0| / (max|b0| + 1) is below
# tolerance after at least minIterations, or at maxIterations; return the iterate with the best objective.
# With stableStop, stop instead when the support has not changed for that many iterations (warm-up use).
gradientEngine <- function(pb, K, ctl, from = NULL, stableStop = NULL, maxIterations = ctl$maxIterations) {
  cur <- if (is.null(from)) refit(pb, integer(0), rep(0, pb$p)) else from
  if (is.null(cur)) stop("The fit with only the unpenalized covariates failed: ", utils::tail(pb$failures, 1))
  pen <- pb$pen; un <- pb$unpenalized
  b <- cur$b
  g <- gradient(pb, b)
  S <- pen[order(-abs(g[pen]))[seq_len(K)]]
  gStep <- g
  gStep[setdiff(pen, S)] <- 0 # MendelIHT projects the first score to the starting support
  project <- function(q) {
    keep <- pen[order(-abs(q[pen]))[seq_len(K)]]
    q[setdiff(pen, keep)] <- 0
    q
  }
  prevF <- -Inf; best <- list(F = -Inf, b = b); converged <- FALSE; stable <- 0; it <- 0
  Fb <- objective(pb, b)
  for (it in seq_len(maxIterations)) {
    if (Fb > best$F) best <- list(F = Fb, b = b)
    if (is.null(stableStop) && it >= maxIterations) break # as fit_iht!: stop at the start of iteration max_iter
    b0 <- b
    v <- c(un, S)
    gv <- gStep[v]
    den <- as.numeric(t(gv) %*% fisher(pb, b, v) %*% gv) + pb$lam * sum(gStep[intersect(v, pen)]^2)
    eta <- sum(gv^2) / den
    if (!is.finite(eta)) eta <- 1e-8
    bn <- project(b0 + eta * gStep); Fn <- objective(pb, bn)
    h <- 0
    while (prevF > Fn && h < ctl$maxHalve) {
      eta <- eta / 2
      bn <- project(b0 + eta * gStep); Fn <- objective(pb, bn)
      h <- h + 1
    }
    Sn <- pen[bn[pen] != 0]
    stable <- if (setequal(Sn, S)) stable + 1 else 0
    b <- bn; S <- Sn; Fb <- Fn; prevF <- Fn
    gStep <- gradient(pb, b)
    if (!is.null(stableStop)) {
      if (stable >= stableStop) { converged <- TRUE; break }
    } else {
      scaled <- max(abs(b - b0)) / (max(abs(b0)) + 1)
      if (it >= ctl$minIterations && scaled < ctl$tolerance) { converged <- TRUE; break }
    }
  }
  if (Fb > best$F) best <- list(F = Fb, b = b)
  list(b = best$b, T = pen[best$b[pen] != 0], F = best$F, iterations = it, converged = converged, steps = it,
       splices = 0)
}
