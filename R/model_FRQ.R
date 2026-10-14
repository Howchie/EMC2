# ============================================================================
# FRQ: finite-reservoir quorum process, fitted in counter coordinates
#
# Generatively, N units are independently available with probability p and,
# when available, register after delays with mean rate lambda; the
# accumulator responds at the K-th registration.  With alpha = K and
# beta = N - K + 1,
#   F(x) = I_{q(x)}(alpha, beta),  q(x) = p G(x),
#   G(x) = 1 - (1 + clock^2 lambda x)^(-1/clock^2)   (1 - exp(-lambda x) at 0).
# The same law is a Poisson counter with threshold K, run on the clock
# phi = q / (1 - q), whose trial input rate is Gamma distributed with mean
# v = beta p lambda and SD sv = v/sqrt(beta).  The fitted coordinates are
# (K, v, sv, h, t0, delta, clock); beta, p and lambda are derived in the
# kernel (src/model_FRQ.h), and sv = 0 is the exact no-variability limit.
# ============================================================================

.frq_check_cols <- function(pars) {
  need <- c("K", "v", "sv", "h", "t0")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("FRQ requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

.frq_opt <- function(pars, ok, nm, default) {
  if (nm %in% colnames(pars)) pars[ok, nm] else default
}

dFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dfrq(t = dt[ok], K = pars[ok, "K"], v = pars[ok, "v"],
                    sv = pars[ok, "sv"], h = pars[ok, "h"],
                    delta = .frq_opt(pars, ok, "delta", 0),
                    clock = .frq_opt(pars, ok, "clock", 1))
  }
  out
}

pFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately kept: the CDF there is h (the complement of the
  # never-finish mass), not one.
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pfrq(t = dt[ok], K = pars[ok, "K"], v = pars[ok, "v"],
                    sv = pars[ok, "sv"], h = pars[ok, "h"],
                    delta = .frq_opt(pars, ok, "delta", 0),
                    clock = .frq_opt(pars, ok, "clock", 1))
  }
  out
}

# Survivor 1 - F, saturating at the defective mass 1 - h.  Computed through the
# reflection I_{1-q}(beta, alpha) (or the upper incomplete gamma at sv = 0)
# rather than by subtraction, so the deep tail stays accurate.
sFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  out <- rep(1, length(dt))
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- pfrq(t = dt[ok], K = pars[ok, "K"], v = pars[ok, "v"],
                    sv = pars[ok, "sv"], h = pars[ok, "h"],
                    delta = .frq_opt(pars, ok, "delta", 0),
                    clock = .frq_opt(pars, ok, "clock", 1),
                    lower_tail = FALSE)
  }
  out
}

# Inverse registration clock: G^{-1} for the unit CDF G with rate lambda and
# shape clock, given the unit survivor sr = 1 - G.
.frq_clock_inv <- function(sr, lambda, clock) {
  c2 <- clock * clock
  out <- -log(sr) / lambda
  pw <- c2 > 0
  out[pw] <- expm1(-c2[pw] * log(sr[pw])) / (c2[pw] * lambda[pw])
  out
}

# Pure-R reference simulator; .rfun_FRQ uses the C++ kernel unless
# options(emc2.cpp_rfun = FALSE).  Exact: the latent quorum percentile
# U ~ Beta(K, beta) satisfies U <= q(x) <=> T <= x, so T = G^{-1}(U / p), and
# U > p is an omission.  It is drawn as its complement W = 1 - U ~ Beta(beta, K)
# by inversion, because for beta < 1 much of U lies within rounding of one:
# V < 1 - h is the omission, and otherwise 1 - G(T) = (W - (1 - p)) / p.
# At sv = 0 the latent is the threshold epoch
# Y ~ Gamma(K): T = Y / v when h = 1, and T = G^{-1}(Y / E) when Y <= E.
# With delta > 0 the latent is drawn on the CDF scale and pulled back through
# H.  An all-Inf trial column is the package's omission convention.
rFRQ <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  .frq_check_cols(pars)
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop("FRQ: pars rows must be a multiple of accumulators")
  if (length(ok) != nrow(pars))
    stop("FRQ: ok must have one value per parameter row")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]
  idx <- which(ok & is.finite(t0))
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    dl <- if ("delta" %in% colnames(p)) p[, "delta"] else rep(0, nrow(p))
    cu <- if ("clock" %in% colnames(p)) p[, "clock"] else rep(1, nrow(p))
    g <- frq_rate(p[, "K"], p[, "v"], p[, "sv"], p[, "h"], dl, cu)
    valid <- !is.na(g[, "p"])
    pois <- valid & !is.finite(g[, "beta"])
    beta_r <- valid & !pois
    vary <- is.finite(dl) & dl > 0
    hit <- rep(Inf, nrow(p))
    if (any(beta_r)) {
      i <- which(beta_r)
      V <- runif(length(i))
      live <- !(p[i, "h"] < 1 & V < 1 - p[i, "h"])
      zz <- V
      iv <- vary[i]
      if (any(iv)) zz[iv] <- frq_h_inv_r(V[iv], dl[i][iv])
      W <- qbeta(zz, g[i, "beta"], p[i, "K"])
      om_p <- ifelse(p[i, "h"] == 1, 0,
                     g[i, "beta"] * g[i, "p"] / g[i, "Etot"])
      sr <- pmin(pmax((W - om_p) / g[i, "p"], .Machine$double.xmin), 1)
      hit[i[live]] <- .frq_clock_inv(sr[live], g[i, "lambda"][live], cu[i][live])
    }
    if (any(pois)) {
      i <- which(pois)
      Y <- rgamma(length(i), p[i, "K"])
      iv <- i[vary[i]]
      if (length(iv))
        Y[vary[i]] <- qgamma(frq_h_inv_r(runif(length(iv)), dl[iv]), p[iv, "K"])
      full <- p[i, "h"] == 1
      hit[i[full]] <- Y[full] / p[i[full], "v"]
      j <- which(!full)
      live <- Y[j] <= g[i[j], "Etot"]
      hit[i[j][live]] <- .frq_clock_inv(1 - Y[j][live] / g[i[j], "Etot"][live],
                                        g[i[j], "lambda"][live],
                                        cu[i[j]][live])
    }
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  # Invalid/non-finite t0 rows stay at +Inf rather than poisoning the race
  # column with NA when the additive shift is applied.
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
  # A trial remains simulable when at least one accumulator is valid; using
  # only the first row here made omission depend on response-level ordering.
  ok_col <- colSums(matrix(ok, nrow = nr)) > 0L
  out <- data.frame(R = rep(NA_character_, n_trials), rt = rep(NA_real_, n_trials))
  out$R[ok_col] <- levels(lR)[R][ok_col]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt[ok_col] <- rt[ok_col]
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Finite Reservoir Quorum Process (FRQ)
#'
#' Model file to estimate the Finite Reservoir Quorum (FRQ) process in EMC2.
#'
#' @details
#'
#' Each accumulator's decision time `u = rt - t0` has distribution function
#'
#' `F(u) = I_q(u)(K, beta)`, with `q(u) = p * G(u)`,
#'
#' where `I` is the regularized incomplete beta function and `G` the unit
#' registration law
#' `G(u) = 1 - (1 + clock^2 * lambda * u)^(-1/clock^2)`, or
#' `1 - exp(-lambda * u)` at `clock = 0`. The kernel derives
#' `beta = (v/sv)^2`, `p` from `h` by a beta quantile, and
#' `lambda = v / (beta * p)`. At `sv = 0` it uses the exact limit
#' `F(u) = P_K(v * u)` when `h = 1` and `F(u) = P_K(E * G(u))` otherwise, where
#' `P_K` is the regularized lower incomplete gamma function,
#' `E = P_K^{-1}(h)` and `lambda = v / E`. The kernel also uses this limit
#' whenever `sv / v < 1e-4`. An accumulator never finishes with probability
#' `1 - h`; a race omission has probability `prod(1 - h)` over accumulators.
#' The derivations are in `Math/CounterAndQuorumModels.tex`.
#'
#' | **Parameter** | **Transform** | **Bounds** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *K* | 1 + exp | \[1, Inf\] | log(2) (K = 3) | *alpha* = *K* | Event threshold |
#' | *v* | log | \[1e-4, Inf\] | log(6) | *lambda* = *v*/(*beta* *p*) | Mean initial input rate |
#' | *sv* | log | \[1e-4, Inf\] or 0 | log(0) | *beta* = (*v*/*sv*)^2 | Between-trial SD of the input rate |
#' | *h* | probit | \[1e-6, 1\] | qnorm(1) | *p* = qbeta(*h*, *K*, *beta*) | Probability that the accumulator ever responds |
#' | *t0* | log | \[0.05, Inf\] or 0 | log(0) | | Non-decision time |
#' | *delta* | log | \[1e-4, 6\] or 0 | log(0) | | Half-width of the log-odds threshold shift |
#' | *clock* | log | \[1e-4, Inf\] or 0 | log(1) | | Shape of the registration law |
#'
#' The values `K = 1`, `sv = 0`, `h = 1`, `t0 = 0`, `delta = 0` and
#' `clock = 0` are exact bound exceptions and can be set as constants.
#' `K` is real-valued.
#'
#' With `delta > 0` the model is `F(u) = H(I_q(u)(K, beta))`, the closed-form
#' marginal over a shift `eps ~ Uniform(-delta, delta)` of the log-odds of the
#' quorum percentile, and `p` is derived from `H^{-1}(h)`. `eps` is drawn
#' independently for each accumulator on each trial, including when `delta`
#' is shared across accumulators.
#'
#' `mapped_pars()` adds the derived columns `beta`, `p`, `lambda`, `Etot`
#' (expected total evidence, `beta * p / (1 - p)`, or `E` at `sv = 0`; `Inf`
#' at `h = 1`), `cv` = `sv / v` and `sQ` = `delta / sqrt(3)`.
#'
#' Optional nuisance parameters are `pContaminant`, an omission probability
#' outside the accumulators, and `pGuess`, a uniform-outlier probability.
#'
#' @return A model list defining the FRQ race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' # The design acts on the input rate and its SD; the threshold is shared
#' # across accumulators.
#' design_FRQ <- design(data = forstmann, model = FRQ, matchfun = matchfun,
#'                      formula = list(K ~ 1, v ~ lM + E, sv ~ lM, h ~ 1,
#'                                     t0 ~ 1),
#'                      contrasts = list(v = list(lM = ADmat)),
#'                      constants = c(h = qnorm(1))) # non-defective race
#' # For all parameters that are not defined in the formula, default values are
#' # assumed (see Table above).
#' @export
FRQ <- function() {
  # K on an OFFSET log scale, K = 1 + exp(x): the sampler coordinate is
  # unbounded while K keeps its K >= 1 support, below which the density is
  # singular at t0 and the likelihood unbounded as t0 -> min(rt).
  p_types <- c("K" = log(2), "v" = log(6),
               # sv = 0 (log(0)) is the exact no-rate-variability counter.
               "sv" = log(0),
               "h" = qnorm(1), "t0" = log(0),
               "delta" = log(0),
               # A constant-rate evidence clock.
               "clock" = log(1))
  transform <- c(K = "exp", v = "exp", sv = "exp", h = "pnorm", t0 = "exp",
                 delta = "exp", clock = "exp")
  transform_lower <- c(K = 1)
  # sv has no upper bound: sv > v (beta < 1) is a valid, very variable rate,
  # and the kernel solves 1 - p from its own quantile.  Its lower guard keeps
  # a sampled sv off the exact sv = 0 branch, which is reachable as a constant.
  # The free h upper bound stops short of one so that a sampler can keep a
  # finite defective mass; h = 1 is reachable through the exception.
  # delta is capped at 6: cosh(delta/2)^2 is already 101 there, and H^{-1}
  # shrinks 1 - h by delta/sinh(delta), so a bigger cap pushes p's complement
  # under double precision at h = 1 - 1e-9.
  minmax <- cbind(K = c(1, Inf), v = c(1e-4, Inf), sv = c(1e-4, Inf),
                  h = c(1e-6, 1), t0 = c(0.05, Inf), delta = c(1e-4, 6),
                  clock = c(1e-4, Inf))
  exception <- c(t0 = 0, delta = 0, clock = 0, h = 1, K = 1, sv = 0)

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  # pContaminant is redundant with 1 - h unless an omission mechanism outside
  # the decision process is wanted.
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = "FRQ",
    p_types = p_types,
    p_types_canonical = c("K", "v", "sv", "h", "t0"),
    transform = list(func = transform, lower = transform_lower),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: Ttransform does not run on the compiled likelihood
      # path (see src/utils.h); the kernel derives the same quantities itself.
      dl <- if ("delta" %in% colnames(pars)) pars[, "delta"] else rep(0, nrow(pars))
      cu <- if ("clock" %in% colnames(pars)) pars[, "clock"] else rep(1, nrow(pars))
      g <- frq_rate(pars[, "K"], pars[, "v"], pars[, "sv"], pars[, "h"], dl, cu)
      # Built as a bare matrix first: cbind()ing named scalars onto a one-row
      # `pars` makes R invent row names from the argument names, which would
      # then travel with the parameter matrix.
      add <- cbind(beta = as.numeric(g[, "beta"]), p = as.numeric(g[, "p"]),
                   lambda = as.numeric(g[, "lambda"]), Etot = as.numeric(g[, "Etot"]),
                   cv = as.numeric(pars[, "sv"] / pars[, "v"]),
                   sQ = as.numeric(dl) / sqrt(3))
      rownames(add) <- NULL
      cbind(pars, add)
    },
    rfun = function(data, pars) {
      .rfun_FRQ(data$lR, pars, ok = attr(pars, "ok"))
    },
    dfun = function(rt, pars) dFRQ(rt, pars),
    pfun = function(rt, pars) pFRQ(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("FRQ: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
