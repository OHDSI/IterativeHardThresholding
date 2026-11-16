# @file CrossValidation.R
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

# Cross-validation over K x variance. All folds use the one data object, with weight 0 on the held-out rows, and run
# the whole selection on their training rows. Within a fold: variances from small to large (Inf last), K from small
# to large, each fit warm-started from the previous one. Criterion: held-out log-likelihood per row, averaged over
# folds. variance = "tune": the default grid, extended past an edge winner until the criterion turns down.

defaultVarianceGrid <- c(0.01, 0.03, 0.1, 0.3, 1, 3, Inf)

crossValidate <- function(cyclopsData, K, variance, ctl, control, weights, unpenalized, forceIntercept, cvSettings,
                          foldThreads = 1) {
  n <- Cyclops::getNumberOfRows(cyclopsData)
  w <- if (is.null(weights)) rep(1, n) else weights
  rowIds <- rowIdsOf(cyclopsData, control, w)
  fold <- assignFolds(cyclopsData, rowIds, w, cvSettings)
  nFolds <- max(fold, na.rm = TRUE)
  grid <- if (identical(variance, "tune")) defaultVarianceGrid else variance
  results <- list()   # key "fold|variance|K" -> list(b, T, heldOutLogLik, nHeldOut)

  # One fold over the variances vs and all K: warm starts stay inside the fold, from its previous K, else from its
  # solution at the nearest variance already fitted (as Cyclops's own cross-validation warm-starts its folds)
  fitFold <- function(f, vs) {
    wf <- w * (!is.na(fold) & fold != f)
    ho <- w * (!is.na(fold) & fold == f)
    out <- list()
    get <- function(key) if (!is.null(out[[key]])) out[[key]] else results[[key]]
    for (v in sort(vs)) {
      pb <- newProblem(cyclopsData, unpenalized, v, control, wf, forceIntercept)
      prev <- NULL
      for (k in K) {
        key <- paste(f, v, k, sep = "|")
        if (!is.null(get(key))) { prev <- get(key); next }
        warm <- prev
        if (is.null(warm)) {
          keys <- strsplit(as.character(c(names(results), names(out))), "|", fixed = TRUE)
          tried <- as.numeric(vapply(Filter(function(z) z[1] == f && z[3] == k, keys), `[`, "", 2))
          if (length(tried)) {
            lv <- function(x) log(pmin(x, 1e6))
            warm <- get(paste(f, tried[which.min(abs(lv(tried) - lv(v)))], k, sep = "|"))
          }
        }
        res <- runEngine(pb, k, ctl, warm = if (ctl$exactRefits) warm else NULL)
        hll <- evaluate(pb, res$b, weights = ho)$logLikelihood
        out[[key]] <- list(b = res$b, T = res$T, heldOutLogLik = hll, nHeldOut = sum(ho))
        prev <- out[[key]]
      }
    }
    out
  }

  # Folds run in parallel processes (forked, so the data are shared) when the caller's control allows several threads,
  # at most one per fold, as Cyclops runs its cross-validation folds on copies of the model in parallel threads.
  # Windows has no fork, so the folds run one after another there.
  nParallel <- if (.Platform$OS.type == "unix") max(1L, min(nFolds, as.integer(foldThreads))) else 1L
  evaluateGrid <- function(vs) {
    outs <- if (nParallel > 1) {
      parallel::mclapply(seq_len(nFolds), fitFold, vs = vs, mc.cores = nParallel, mc.preschedule = FALSE)
    } else {
      lapply(seq_len(nFolds), fitFold, vs = vs)
    }
    for (o in outs) {
      if (inherits(o, "try-error") || !is.list(o)) stop("A cross-validation fold failed: ", paste(o, collapse = " "))
      results[names(o)] <<- o
    }
  }

  summarize <- function() {
    keys <- strsplit(names(results), "|", fixed = TRUE)
    d <- data.frame(fold = as.integer(sapply(keys, `[`, 1)), variance = as.numeric(sapply(keys, `[`, 2)),
                    K = as.integer(sapply(keys, `[`, 3)),
                    ll = sapply(results, function(r) r$heldOutLogLik / max(1, r$nHeldOut)))
    agg <- do.call(rbind, lapply(split(d, list(d$K, d$variance), drop = TRUE), function(z) {
      data.frame(K = z$K[1], variance = z$variance[1], heldOutLogLik = mean(z$ll),
                 se = if (nrow(z) > 1) stats::sd(z$ll) / sqrt(nrow(z)) else NA_real_, folds = nrow(z))
    }))
    agg[order(agg$variance, agg$K), ]
  }

  evaluateGrid(grid)
  tab <- summarize()
  if (identical(variance, "tune")) {
    for (side in c("low", "high")) {
      for (step in 1:4) {
        finite <- tab[is.finite(tab$variance), ]
        best <- finite[which.max(finite$heldOutLogLik), ]
        if (side == "low" && best$variance == min(finite$variance)) {
          evaluateGrid(best$variance / 3)
        } else if (side == "high" && best$variance == max(finite$variance) && best$variance < 1000 &&
                   best$heldOutLogLik > max(tab$heldOutLogLik[!is.finite(tab$variance)], -Inf)) {
          evaluateGrid(best$variance * 3)
        } else {
          break
        }
        tab <- summarize()
      }
    }
  }

  best <- tab[which.max(tab$heldOutLogLik), ]
  chosen <- best
  if (isTRUE(cvSettings$oneStandardError) && is.finite(best$se)) {
    ok <- tab[tab$heldOutLogLik >= best$heldOutLogLik - best$se, ]
    ok <- ok[ok$K == min(ok$K), ]
    chosen <- ok[which.max(ok$heldOutLogLik), ]
  }
  tab$chosen <- tab$K == chosen$K & tab$variance == chosen$variance

  predictions <- NULL
  if (isTRUE(cvSettings$returnCvPredictions)) {
    predictions <- do.call(rbind, lapply(seq_len(nFolds), function(f) {
      r <- results[[paste(f, chosen$variance, chosen$K, sep = "|")]]
      ho <- which(!is.na(fold) & fold == f & w > 0)
      pb <- newProblem(cyclopsData, unpenalized, chosen$variance, control, NULL, forceIntercept)
      fit <- finalFit(pb, r$b, weights = NULL)
      data.frame(rowId = rowIds[ho], fold = f, prediction = as.numeric(predict(fit))[ho])
    }))
    noFold <- which(is.na(fold) | w == 0)
    if (length(noFold) > 0) {
      predictions <- rbind(predictions, data.frame(rowId = rowIds[noFold], fold = NA_integer_, prediction = NA_real_))
    }
  }
  list(K = chosen$K, variance = chosen$variance, table = tab, predictions = predictions, folds = fold)
}

# Row ids in the data's row order (from a quick fit with all coefficients fixed at 0)
rowIdsOf <- function(cyclopsData, control, w) {
  fit <- suppressWarnings(Cyclops::fitCyclopsModel(cyclopsData, prior = Cyclops::createPrior("none"), control = control,
                                                   startingCoefficients = rep(0, Cyclops::getNumberOfCovariates(cyclopsData)),
                                                   fixedCoefficients = rep(TRUE, Cyclops::getNumberOfCovariates(cyclopsData)),
                                                   warnings = FALSE))
  ids <- names(predict(fit))
  num <- suppressWarnings(as.numeric(ids))
  if (all(!is.na(num))) num else ids
}

# Folds from cvFolds (data.frame(rowId, fold)) or at random.
# Rows with weight 0 get no fold.
assignFolds <- function(cyclopsData, rowIds, w, cvSettings) {
  n <- length(rowIds)
  if (!is.null(cvSettings$cvFolds)) {
    fold <- cvSettings$cvFolds$fold[match(rowIds, cvSettings$cvFolds$rowId)]
    if (any(is.na(fold) & w > 0)) stop("cvFolds does not assign every row with positive weight to a fold")
    fold <- as.integer(factor(fold))
  } else {
    if (!is.null(cvSettings$seed)) set.seed(cvSettings$seed)
    k <- cvSettings$folds
    fold <- sample(rep(seq_len(k), length.out = n)) # like Cyclops's own fold assignment (not stratified)
  }
  fold[w == 0] <- NA_integer_
  fold
}
