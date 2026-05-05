# Wrapper functions for internal summary methods ------------------------------

#' R-hat convergence diagnostic
#'
#' Computes the potential scale reduction factor (\eqn{\hat{R}}) for each
#' parameter from an MCMC run. Supports both the legacy [coda::gelman.diag()]
#' implementation (Gelman & Rubin, 1992; Brooks & Gelman, 1998) and the
#' rank-normalized, folded implementation from Vehtari et al. (2021).
#'
#' @param mcmc_list A `coda::mcmc.list` object containing MCMC draws.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::gelman.diag()] with split chains. `"new"` uses a vendored
#'    implementation of rank-normalized, folded Rhat from Vehtari et al. (2021).
#' @param omit_mpsrf Boolean, only relevant if `version = "old"`. If `TRUE`
#'    (default) the multivariate PSRF is not returned.
#'
#' @details
#' The `"old"` method produces the traditional Gelman-Rubin diagnostic but
#' uses split chains for improved detection of non-stationarity.
#'
#' The `"new"` method implements rank-normalized, folded, split-chain Rhat
#' (Vehtari et al., 2021), which better detects convergence failure in the
#' presence of heavy tails or non-stationarity. Code is adapted from the
#' \pkg{posterior} package under BSD 3-Clause License (compatible with GPL-3).
#'
#' @return A named numeric vector of \eqn{\hat{R}} values. If
#'   `version = "old"` and `omit_mpsrf = FALSE`, the last element is `"mpsrf"`.
#'
#' @references
#' Gelman, A., & Rubin, D.B. (1992). Inference from iterative simulation using
#' multiple sequences. *Statistical Science*, *7*, 457-511.
#'
#' Brooks, S.P., & Gelman, A. (1998). General methods for monitoring convergence
#' of iterative simulations. *Journal of Computational and Graphical Statistics*,
#' *7*, 434-455.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667-718.
#'
#' @keywords internal
r_hat <- function(mcmc_list, version = c("old", "new"), omit_mpsrf = TRUE) {
  version <- match.arg(version)
  if (version == "old") {
    result <- gelman_diag_robust(mcmc_list, omit_mpsrf = omit_mpsrf)
  } else {
    mcmc_mats <- prep_mcmc_diagnostics(mcmc_list)
    result <- vapply(
      X = mcmc_mats,
      FUN = function(x) {
        out <- try(rhat_new(x), silent = TRUE)
        if (is(out, "try-error")) return(Inf)
        return(out)
      },
      FUN.VALUE = numeric(1)
    )
  }
  return(result)
}

#' Effective sample size (ESS)
#'
#' Computes the effective sample size for each parameter from an MCMC run.
#' Supports both [coda::effectiveSize()] and the improved version from
#' Vehtari et al. (2021).
#'
#' @param mcmc_list A `coda::mcmc.list` object containing MCMC draws.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::effectiveSize()]. `"new"` uses a vendored implementation of
#'    rank-normalized bulk ESS from Vehtari et al. (2021).
#'
#' @return A named numeric vector of effective sample sizes.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667-718.
#'
#' @keywords internal
n_eff <- function(mcmc_list, version = c("old", "new")) {
  version <- match.arg(version)
  if (version == "old") {
    result <- coda::effectiveSize(mcmc_list)
  } else {
    mcmc_mats <- prep_mcmc_diagnostics(mcmc_list)
    result <- vapply(
      X = mcmc_mats,
      FUN = function(x) {
        out <- try(ess_basic(x), silent = TRUE)
        if (is(out, "try-error")) return(0)
        return(out)
      },
      FUN.VALUE = numeric(1)
    )
  }
  return(result)
}

#' @noRd
prep_mcmc_diagnostics <- function(mcmc_list) {
  stopifnot(is(mcmc_list, "mcmc.list"))
  n_chains <- length(mcmc_list)
  n_iter <- unique(vapply(mcmc_list, nrow, integer(1)))
  if (length(n_iter) > 1L) {
    stop("Chains have unequal numbers of iterations; please trim or pad first.")
  }
  n_iter <- n_iter[1]
  param_names <- colnames(mcmc_list[[1]])
  result <- setNames(vector("list", length(param_names)), param_names)
  for (param in param_names) {
    mat <- matrix(NA_real_, nrow = n_iter, ncol = n_chains)
    for (chain in seq_len(n_chains)) {
      mat[ , chain] <- mcmc_list[[chain]][ , param]
    }
    result[[param]] <- mat
  }
  return(result)
}

# Old Rhat convergence diagnostic ---------------------------------------------

gelman_diag_robust <- function(
    mcl,
    autoburnin = FALSE,
    transform = TRUE,
    omit_mpsrf = TRUE
) {
  mcl <- split_mcl(mcl)
  gd <- try(
    coda::gelman.diag(
      mcl,
      autoburnin = autoburnin,
      transform = transform,
      multivariate = !omit_mpsrf
    ),
    silent = TRUE
  )
  if (is(gd, "try-error")) {
    if (omit_mpsrf) return(list(psrf = matrix(Inf)))
    else            return(list(psrf = matrix(Inf), mpsrf = Inf))
  }
  gd_out <- gd[[1]][ , 1]  # drop CI column
  if (!omit_mpsrf) {
    gd_out <- c(gd_out, gd[["mpsrf"]])
    names(gd_out)[length(gd_out)] <- "mpsrf"
  }
  return(gd_out)
}

split_mcl <- function(mcl) {
  if (!is.list(mcl)) mcl <- list(mcl)
  mcl2 <- mcl
  half <- floor(unlist(lapply(mcl, nrow)) / 2)
  for (i in seq_along(half)) {
    mcl[[i]]  <- coda::as.mcmc(mcl[[i]][1:half[i], ])
    mcl2[[i]] <- coda::as.mcmc(mcl2[[i]][(half[i] + 1):(2 * half[i]), ])
  }
  coda::as.mcmc.list(c(mcl, mcl2))
}

# -----------------------------------------------------------------------------
# All of the following code was adapted from the `posterior` package
# (posterior/R/convergence.R), implementing methods from Vehtari et al. (2021).
# https://doi.org/10.1214/20-BA1221
#
# The `posterior` package is released under BSD 3-Clause License, which is
# compatible with the GPL-3 license of this package. The original BSD license
# text is reproduced below.
#
# Copyright (c) 2021, posterior package authors;
# Stan Developers and their Assignees; Trustees of Columbia University
#
# BSD 3-Clause License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------

# Rhat convergence diagnostic -------------------------------------------------

#' @noRd
rhat_new <- function(x) {
  rhat_bulk <- rhat_basic(rank_normalise(split_chains(x)))
  rhat_tail <- rhat_basic(rank_normalise(split_chains(fold_draws(x))))
  max(rhat_bulk, rhat_tail)
}

#' @noRd
rhat_basic <- function(x) {
  if (has_bad_draws(x)) return(NA_real_)
  n_iter      <- NROW(x)
  chain_mean  <- colMeans(x)
  chain_var   <- apply(x, 2, stats::var)
  var_between <- n_iter * stats::var(chain_mean)
  var_within  <- mean(chain_var)
  sqrt((var_between / var_within + n_iter - 1) / n_iter)
}


# Effective Sample Size (ESS) diagnostics -------------------------------------

#' @noRd
ess_basic <- function(x) {
  ess(split_chains(x))
}

#' @noRd
ess_bulk <- function(x) {
  ess(rank_normalise(split_chains(x)))
}

#' @noRd
ess_tail <- function(x) {
  min(ess_quantile(x, probs = 0.05), ess_quantile(x, probs = 0.95))
}

#' @noRd
ess_mean <- function(x) {
  ess(split_chains(x))
}

#' @noRd
ess_median <- function(x) {
  ess_quantile(x, probs = 0.5, names = FALSE)
}

#' @noRd
ess_sd <- function(x) {
  ess(split_chains((x - mean(x))^2))
}

#' @noRd
ess_quantile <- function(x, probs, names = TRUE) {
  probs  <- check_quantile_probs(probs)
  result <- unlist(lapply(probs, ess_quantile_engine, x = x))
  if (names) names(result) <- paste0("ess_q", probs * 100)
  result
}

#' @noRd
ess_quantile_engine <- function(x, prob) {
  if (has_bad_draws(x)) return(NA_real_)
  if (prob == 1) prob <- (length(x) - 0.5) / length(x)
  ess(split_chains(x <= stats::quantile(x, prob)))
}

#' @noRd
ess <- function(x) {
  n_chain   <- NCOL(x)
  n_iter    <- NROW(x)
  if (n_iter < 3L || has_bad_draws(x)) return(NA_real_)
  n_samples <- n_chain * n_iter

  acov      <- apply(x, 2, autocovariance)
  acov_means <- rowMeans(acov)
  mean_var  <- acov_means[1] * n_iter / (n_iter - 1)
  var_plus  <- mean_var * (n_iter - 1) / n_iter
  if (n_chain > 1) var_plus <- var_plus + stats::var(colMeans(x))

  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_iter)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd  <- 1 - (mean_var - acov_means[t + 2]) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while ((t < NROW(acov) - 5) && !is.nan(rho_hat_even + rho_hat_odd) &&
         (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even <- 1 - (mean_var - acov_means[t + 1]) / var_plus
    rho_hat_odd  <- 1 - (mean_var - acov_means[t + 2]) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t
  if (rho_hat_even > 0) rho_hat_t[max_t + 1] <- rho_hat_even

  # Geyer's initial monotone sequence
  t <- 0
  while (t <= (max_t - 4)) {
    t <- t + 2
    if ((rho_hat_t[t + 1] + rho_hat_t[t + 2]) > (rho_hat_t[t - 1] + rho_hat_t[t])) {
      rho_hat_t[t + 1] <- (rho_hat_t[t - 1] + rho_hat_t[t]) / 2
      rho_hat_t[t + 2] <- rho_hat_t[t + 1]
    }
  }

  # Improved truncated estimate
  tau_hat   <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t + 1]
  tau_bound <- 1 / log10(n_samples)
  if (tau_hat < tau_bound) {
    warning("The ESS has been capped to avoid unstable estimates.")
    tau_hat <- tau_bound
  }
  n_samples / tau_hat
}


# Helper functions ------------------------------------------------------------

#' @noRd
rank_normalise <- function(x) {
  r <- rank(as.array(x), ties.method = "average")
  c <- 3 / 8
  p <- (r - c) / (length(r) - 2 * c + 1)
  z <- stats::qnorm(p)
  z[is.na(x)] <- NA_real_
  if (!is.null(dim(x))) z <- array(z, dim = dim(x), dimnames = dimnames(x))
  z
}

#' @noRd
split_chains <- function(x) {
  niter <- NROW(x)
  if (niter == 1L) return(x)
  half <- niter / 2
  cbind(x[1:floor(half), ], x[ceiling(half + 1):niter, ])
}

#' @noRd
fold_draws <- function(x) {
  abs(x - median(x))
}

#' @noRd
has_bad_draws <- function(x, tol = .Machine$double.eps) {
  anyNA(x) || any(is.infinite(x)) || abs(max(x) - min(x)) < tol
}

#' @noRd
autocovariance <- function(x) {
  N    <- length(x)
  var_x <- stats::var(x)
  if (var_x == 0) return(rep(0, N))
  M        <- stats::nextn(N)
  M_double <- 2 * M
  x_c      <- c(x - mean(x), rep.int(0, M_double - N))
  ac <- Re(stats::fft(abs(stats::fft(x_c))^2, inverse = TRUE)[1:N])
  ac / ac[1] * var_x * (N - 1) / N
}

#' @noRd
check_quantile_probs <- function(probs) {
  probs <- as.numeric(probs)
  if (any(probs < 0 | probs > 1)) stop("'probs' must contain values between 0 and 1.")
  probs
}
