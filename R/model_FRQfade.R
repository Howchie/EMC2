# FRQ with a fading evidence-registration opportunity.

.frqfade_check_cols <- function(pars) {
  need <- c("alpha", "beta", "lambda", "kappa", "t0")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("FRQfade requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

dFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- dfrqfade(dt[ok], pars[ok, "alpha"], pars[ok, "beta"],
      pars[ok, "lambda"], pars[ok, "kappa"], delta, cv_u)
  }
  out
}

pFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- pfrqfade(dt[ok], pars[ok, "alpha"], pars[ok, "beta"],
      pars[ok, "lambda"], pars[ok, "kappa"], delta, cv_u)
  }
  out
}

sFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  out <- rep(1, length(dt))
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- pfrqfade(dt[ok], pars[ok, "alpha"], pars[ok, "beta"],
      pars[ok, "lambda"], pars[ok, "kappa"], delta, cv_u,
      lower_tail = FALSE)
  }
  out
}

# Reference generator. The compiled simulator is the default route.
rFRQfade <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  .frqfade_check_cols(pars)
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop("FRQfade: pars rows must be a multiple of accumulators")
  if (length(ok) != nrow(pars))
    stop("FRQfade: ok must have one value per parameter row")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]
  alpha <- pars[, "alpha"]
  beta <- pars[, "beta"]
  lambda <- pars[, "lambda"]
  kappa <- pars[, "kappa"]
  delta <- if ("delta" %in% colnames(pars)) pars[, "delta"] else rep(0, nrow(pars))
  cv_u <- if ("cv_u" %in% colnames(pars)) pars[, "cv_u"] else rep(0, nrow(pars))
  valid <- ok & is.finite(t0) & is.finite(alpha) & alpha >= 1 &
    is.finite(beta) & beta >= 1 & is.finite(lambda) & lambda > 0 &
    is.finite(kappa) & kappa >= 0 & is.finite(delta) & delta >= 0 & delta <= 700 &
    is.finite(cv_u) & cv_u >= 0 & is.finite(cv_u^2)
  valid[is.na(valid)] <- FALSE
  idx <- which(valid)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    d <- delta[idx]
    c_u <- cv_u[idx]
    rank <- stats::rbeta(nrow(p), p[, "alpha"], p[, "beta"])
    varied <- d > 0
    if (any(varied)) {
      rank[varied] <- stats::qbeta(
        frq_h_inv_r(stats::runif(sum(varied)), d[varied]),
        p[varied, "alpha"], p[varied, "beta"])
    }
    c2 <- c_u^2
    y <- -log1p(-rank)
    A <- y / p[, "lambda"]
    frail <- c2 > 0
    if (any(frail)) {
      z <- c2[frail] * y[frail]
      A_frail <- y[frail] / p[frail, "lambda"]
      overflow <- z > log(.Machine$double.xmax)
      A_frail[overflow] <- Inf
      nz <- !overflow & z > 0
      A_frail[nz] <- A_frail[nz] * expm1(z[nz]) / z[nz]
      A[frail] <- A_frail
    }
    kappa <- p[, "kappa"]
    # Unlimited opportunity is proper: an exact rank of one (a zero-probability
    # RNG endpoint) must not become an omission, matching the compiled path.
    edge <- kappa == 0 & rank >= 1
    if (any(edge)) {
      r1 <- 1 - .Machine$double.eps / 2
      A[edge] <- -log1p(-r1) / p[edge, "lambda"]
      e_frail <- edge & c2 > 0
      if (any(e_frail)) {
        z1 <- c2[e_frail] * -log1p(-r1)
        A[e_frail] <- A[e_frail] * expm1(z1) / z1
      }
    }
    hit <- A
    fading <- kappa > 0
    if (any(fading)) {
      kA <- kappa[fading] * A[fading]
      live <- is.finite(kA) & kA < 1
      z <- rep(Inf, length(kA))
      z[live] <- -log1p(-kA[live]) / kappa[fading][live]
      hit[fading] <- z
    }
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  t0_safe <- t0
  t0_safe[!is.finite(t0_safe)] <- 0
  dt <- dt + matrix(t0_safe, nrow = nr)
  bad_col <- colSums(!is.infinite(dt)) == 0L
  R <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(R, seq_len(ncol(dt)))
  rt <- dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  ok_col <- colSums(matrix(ok, nrow = nr)) > 0L
  out <- data.frame(R = rep(NA_character_, n_trials), rt = rep(NA_real_, n_trials))
  out$R[ok_col] <- levels(lR)[R][ok_col]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt[ok_col] <- rt[ok_col]
  .apply_timed_guess_winner(out, levels(lR))
}

#' Finite reservoir quorum with a fading registration opportunity
#'
#' `FRQfade()` is a sibling of [FRQ()] in which response omissions arise from
#' the same registration process that determines finite response times.
#'
#' @details
#' Each accumulator has a quorum percentile with Beta shapes `alpha` and
#' `beta`. Evidence registers at mean rate `lambda`, while its total future
#' opportunity fades exponentially at rate `kappa`. For decision time
#' `u = rt - t0`, the integrated opportunity clock is
#' `A_kappa(u) = (1 - exp(-kappa * u)) / kappa`, with the limit `u` at
#' `kappa = 0`. Thus positive `kappa` limits total integrated opportunity to
#' `1/kappa`; it does not impose a finite cutoff in physical time. A quorum
#' above the available evidence percentile produces an omission.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *alpha* | 1 + exp | \[1, Inf\] | log(0) | | Quorum shape |
#' | *beta* | 1 + exp | \[1, Inf\] | log(0) | | Spare-capacity shape |
#' | *lambda* | log | \[0, Inf\] | log(log(2) / 0.5) | | Mean latent registration rate at time zero |
#' | *kappa* | log | \[0, Inf\] | log(0.5) | | Exponential decay rate of registration opportunity; exact zero is supported |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time |
#' | *delta* | log | \[0, Inf\] | log(0) | | Half-width of uniform log-odds threshold variability |
#' | *cv_u* | log | \[0, Inf\] | log(0) | | CV of independent Gamma-distributed unit rates |
#'
#' The fitted coordinates are `alpha`, `beta`, `lambda`, `kappa`, `t0`,
#' `delta`, and `cv_u`; `h` is derived and is not fitted. At `kappa = 0`,
#' opportunity is unlimited and `h = 1`, including with threshold variation
#' or unit-rate frailty. At positive `kappa`, `h` is the transformed Beta CDF
#' evaluated at the limiting evidence percentile `q_inf`. Both `delta` and
#' `cv_u` can change `h` under finite opportunity. `mapped_pars()` reports
#' `h`, `q_inf`, reservoir summaries `N = alpha + beta - 1` and
#' `d = alpha/N`, threshold SD `sQ = delta/sqrt(3)`, and unit-rate Gamma shape
#' `a_u = 1/cv_u^2` (`Inf` when `cv_u = 0`).
#'
#' With `cv_u = 0`, the model is the proper FRQ curve run on the internal clock
#' `A_kappa(u)`, which stops at `1/kappa`. A response occurs when the
#' corresponding unlimited-opportunity FRQ time is below this horizon, and
#' its physical decision time is `-log1p(-kappa * A)/kappa`. For fixed shapes,
#' `lambda/kappa` controls completion probability and `kappa` sets the time
#' scale, coupling omission probability to the conditional RT distribution.
#' This coupling remains when `alpha = beta = 1`.
#'
#' Optional nuisance parameters are `pContaminant`, for omissions outside the
#' evidence process, and `pGuess`, for uniform response-time outliers. The
#' `delta` and `cv_u` model extensions are included in every positional kernel
#' contract and default to zero.
#'
#' @return A model list defining the FRQfade race model.
#' @examples
#' design_FRQfade <- design(data = forstmann, model = FRQfade,
#'   matchfun = function(d) d$S == d$lR,
#'   formula = list(alpha ~ 1, beta ~ 1, lambda ~ lM + E,
#'                  kappa ~ lM, t0 ~ 1))
#' @export
FRQfade <- function() {
  p_types <- c("alpha" = log(0), "beta" = log(0),
    "lambda" = log(log(2) / 0.5), "kappa" = log(0.5),
    "t0" = log(0), "delta" = log(0), "cv_u" = log(0))
  transform <- c(alpha = "exp", beta = "exp", lambda = "exp",
    kappa = "exp", t0 = "exp", delta = "exp", cv_u = "exp")
  transform_lower <- c(alpha = 1, beta = 1)
  minmax <- cbind(alpha = c(1, Inf), beta = c(1, Inf),
    lambda = c(1e-4, Inf), kappa = c(1e-4, Inf), t0 = c(0.05, Inf),
    delta = c(1e-4, 6), cv_u = c(1e-4, Inf))
  exception <- c(alpha = 1, beta = 1, kappa = 0, t0 = 0,
                 delta = 0, cv_u = 0)
  nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- nuis$p_types
  transform <- nuis$transform
  minmax <- nuis$minmax
  exception <- nuis$exception

  list(
    type = "RACE",
    c_name = "FRQ_FADE",
    p_types = p_types,
    p_types_canonical = c("alpha", "beta", "lambda", "kappa", "t0"),
    transform = list(func = transform, lower = transform_lower),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      delta <- if ("delta" %in% colnames(pars)) pars[, "delta"] else 0
      cv_u <- if ("cv_u" %in% colnames(pars)) pars[, "cv_u"] else 0
      summary <- frq_fade_summary(pars[, "alpha"], pars[, "beta"],
        pars[, "lambda"], pars[, "kappa"], delta, cv_u)
      N <- pars[, "alpha"] + pars[, "beta"] - 1
      add <- cbind(q_inf = as.numeric(summary[, "q_inf"]),
        h = as.numeric(summary[, "h"]), N = as.numeric(N),
        d = as.numeric(pars[, "alpha"] / N),
        sQ = as.numeric(delta) / sqrt(3),
        a_u = ifelse(cv_u == 0, Inf, 1 / cv_u^2))
      rownames(add) <- NULL
      cbind(pars, add)
    },
    rfun = function(data, pars) {
      ok <- attr(pars, "ok")
      if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
      .rfun_FRQfade(data$lR, pars, ok = ok)
    },
    dfun = function(rt, pars) dFRQfade(rt, pars),
    pfun = function(rt, pars) pFRQfade(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("FRQfade: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
