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
# (K, v, sv, h, t0, plo, clock); beta, p and lambda are derived in the
# kernel (src/model_FRQ.h), and sv = 0 is the exact no-variability limit.
# plo is the probability that an accumulator commits at the low threshold
# K_lo = 1 instead of K on a trial.  The accel chart samples the rate CV and
# runs the clock = 1 law on the operational time (exp(accel u) - 1) / accel.
# ============================================================================

.frq_check_cols <- function(pars) {
  need <- c("K", "v", "h", "t0")
  missing <- setdiff(need, colnames(pars))
  if (!any(c("sv", "cv") %in% colnames(pars))) missing <- c(missing, "sv")
  if (length(missing))
    stop("FRQ requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

# Rate SD per row: the sv column, or cv * v under the cv chart.
.frq_sv <- function(pars) {
  if ("sv" %in% colnames(pars)) pars[, "sv"] else pars[, "cv"] * pars[, "v"]
}

# pars with an sv column, for the entry points that take the rate SD.
.frq_with_sv <- function(pars) {
  if ("sv" %in% colnames(pars)) return(pars)
  out <- cbind(pars, sv = .frq_sv(pars))
  attr(out, "ok") <- attr(pars, "ok")
  out
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
                    sv = .frq_sv(pars)[ok], h = pars[ok, "h"],
                    plo = .frq_opt(pars, ok, "plo", 0),
                    clock = .frq_opt(pars, ok, "clock", 1),
                    accel = .frq_opt(pars, ok, "accel", 0))
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
                    sv = .frq_sv(pars)[ok], h = pars[ok, "h"],
                    plo = .frq_opt(pars, ok, "plo", 0),
                    clock = .frq_opt(pars, ok, "clock", 1),
                    accel = .frq_opt(pars, ok, "accel", 0))
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
                    sv = .frq_sv(pars)[ok], h = pars[ok, "h"],
                    plo = .frq_opt(pars, ok, "plo", 0),
                    clock = .frq_opt(pars, ok, "clock", 1),
                    accel = .frq_opt(pars, ok, "accel", 0),
                    lower_tail = FALSE)
  }
  out
}

# Inverse registration clock: G^{-1} for the unit CDF G with rate lambda and
# shape clock, given the log unit survivor log(1 - G).
.frq_clock_inv <- function(log_sr, lambda, clock) {
  c2 <- clock * clock
  out <- -log_sr / lambda
  pw <- c2 > 0
  out[pw] <- expm1(-c2[pw] * log_sr[pw]) / (c2[pw] * lambda[pw])
  out
}

.frq_log1mexp <- function(y) ifelse(y > -log(2), log(-expm1(y)), log1p(-exp(y)))

# log Q_{a,b}(u); as frq_log_qbeta in src/model_FRQ.h, the small-x side uses
# the leading term of I_x(a, b), since Q underflows long before u does.
.frq_log_qbeta <- function(u, a, b) {
  ly0 <- (log(u) + log(a) + lbeta(a, b)) / a
  small <- ly0 < log(1e-17 * a * (a + 1) / (abs(1 - b) + 1))
  out <- ly0
  big <- !small
  if (any(big)) {
    x <- qbeta(u[big], a[big], b[big])
    out[big] <- ifelse(x > 0, log(x), ly0[big])
  }
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
# With plo > 0 each row first draws its threshold, K_lo = 1 with probability
# plo and K otherwise, and the same construction runs at that threshold with
# the shared p.  An all-Inf trial column is the package's omission convention.
rFRQ <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  .frq_check_cols(pars)
  ok <- .frq_check_race(lR, pars, ok, "FRQ")
  dt <- rep(Inf, nrow(pars))
  idx <- which(ok & is.finite(pars[, "t0"]))
  if (length(idx)) dt[idx] <- .frq_draw_dt(.frq_with_sv(pars[idx, , drop = FALSE]))
  .frq_race(lR, dt, pars[, "t0"], ok)
}

.frq_check_race <- function(lR, pars, ok, model) {
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop(model, ": pars rows must be a multiple of accumulators")
  if (length(ok) != nrow(pars))
    stop(model, ": ok must have one value per parameter row")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  ok
}

# FRQ decision times (t0 excluded) for parameter rows p; +Inf is an omission
# or an invalid row.
.frq_draw_dt <- function(p) {
  pl <- if ("plo" %in% colnames(p)) p[, "plo"] else rep(0, nrow(p))
  cu <- if ("clock" %in% colnames(p)) p[, "clock"] else rep(1, nrow(p))
  g <- frq_rate(p[, "K"], p[, "v"], p[, "sv"], p[, "h"], pl, cu)
  valid <- !is.na(g[, "p"])
  pois <- valid & !is.finite(g[, "beta"])
  beta_r <- valid & !pois
  # Threshold per row: K_lo = 1 in the low state.  The draw is made only where
  # plo > 0, so the base model consumes the same uniforms as before.
  kk <- p[, "K"]
  mix <- valid & is.finite(pl) & pl > 0
  if (any(mix)) kk[mix][runif(sum(mix)) < pl[mix]] <- 1
  hit <- rep(Inf, nrow(p))
  if (any(beta_r)) {
    i <- which(beta_r)
    V <- runif(length(i))
    # W = 1 - U ~ Beta(beta, k) and 1 - p are held as logs: both underflow
    # for beta < 1.  W <= 1 - p is the omission.
    lw <- .frq_log_qbeta(V, g[i, "beta"], kk[i])
    log_om_p <- g[i, "log1mp"]
    live <- p[i, "h"] == 1 | lw > log_om_p
    floor <- log(.Machine$double.xmin)
    log_sr <- rep(floor, length(i))
    above <- lw > log_om_p
    log_sr[above] <- lw[above] + .frq_log1mexp(log_om_p[above] - lw[above]) -
      log(g[i, "p"][above])
    log_sr[!is.finite(log_sr)] <- floor
    log_sr <- pmin(log_sr, 0)
    hit[i[live]] <- .frq_clock_inv(log_sr[live], g[i, "lambda"][live],
                                   cu[i][live])
  }
  if (any(pois)) {
    i <- which(pois)
    Y <- rgamma(length(i), kk[i])
    full <- p[i, "h"] == 1
    hit[i[full]] <- Y[full] / p[i[full], "v"]
    j <- which(!full)
    live <- Y[j] <= g[i[j], "Etot"]
    hit[i[j][live]] <- .frq_clock_inv(log1p(-Y[j][live] / g[i[j], "Etot"][live]),
                                      g[i[j], "lambda"][live],
                                      cu[i[j]][live])
  }
  if ("accel" %in% colnames(p)) {
    # Operational time back to physical time, u = log1p(a x) / a.
    a <- p[, "accel"]
    hit[!(is.finite(a) & a >= 0)] <- Inf
    w <- is.finite(hit) & is.finite(a) & a > 0
    hit[w] <- log1p(a[w] * hit[w]) / a[w]
  }
  hit[!is.finite(hit) | hit < 0] <- Inf
  hit
}

# Resolve a race from per-row decision times dt (accumulators fastest within
# each trial column); an all-Inf trial column is the package's omission
# convention.
.frq_race <- function(lR, dt, t0, ok) {
  nr <- length(levels(lR))
  n_trials <- length(dt) / nr
  # Invalid/non-finite t0 rows stay at +Inf rather than poisoning the race
  # column with NA when the additive shift is applied.
  t0_safe <- t0
  t0_safe[!is.finite(t0_safe)] <- 0
  dt <- matrix(dt + t0_safe, nrow = nr)

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
#' | *K* | log | \[2, Inf\] | log(3) | *alpha* = *K* | Event threshold |
#' | *v* | log | \[1e-4, Inf\] | log(6) | *lambda* = *v*/(*beta* *p*) | Mean initial input rate |
#' | *sv* | log | \[1e-4, Inf\] or 0 | log(0) | *beta* = (*v*/*sv*)^2 | Between-trial SD of the input rate |
#' | *h* | probit | \[1e-6, 1\] | qnorm(1) | *p* = qbeta(*h*, *K*, *beta*) | Probability that the accumulator ever responds |
#' | *t0* | log | \[0.05, Inf\] or 0 | log(0) | | Non-decision time |
#' | *plo* | probit | \[1e-6, 1\] or 0 | qnorm(0) | | Probability of the low-threshold state |
#' | *clock* | log | \[1e-4, Inf\] or 0 | log(1) | | Shape of the registration law |
#'
#' The values `sv = 0`, `h = 1`, `t0 = 0`, `plo = 0` and
#' `clock = 0` are exact bound exceptions and can be set as constants.
#' `K` is real-valued. A sampled `K` is `2 + exp(x)`; its link follows the
#' lower bound, so a `bound` argument to [design()] moves both. A constant `K`
#' uses the plain log link.
#'
#' `parameterization = "cv"` samples *cv* = *sv*/*v*, the between-trial CV of
#' the input rate, in place of *sv*:
#'
#' | **Parameter** | **Transform** | **Bounds** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *cv* | log | \[1e-4, Inf\] or 0 | log(0) | *beta* = 1/*cv*^2, *sv* = *cv* *v* | Between-trial CV of the input rate |
#'
#' `cv = 0` is the same exact limit as `sv = 0`. At equal values of *sv* and
#' *cv* *v* the two parameterizations give the same likelihood.
#'
#' `parameterization = "accel"` samples *cv* and replaces *clock* by the clock
#' acceleration *accel*. The registration law is the `clock = 1` law run on the
#' operational time `A(u) = (exp(accel * u) - 1) / accel` (`A(u) = u` at
#' `accel = 0`):
#'
#' `F(u) = I_q(u)(K, beta)`, with `q(u) = p * x(u) / (1 + x(u))` and
#' `x(u) = lambda * A(u)`,
#'
#' so that the Gamma-mixed counter's intensity is `W * (v / beta) *
#' exp(accel * u)` given its trial rate multiplier `W ~ Gamma(beta, 1)`. At
#' `cv = 0` and `h = 1`, `F(u) = P_K(v * A(u))`. `h`, `p` and `plo` enter as in
#' the other parameterizations, since `A(Inf) = Inf`.
#'
#' | **Parameter** | **Transform** | **Bounds** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *cv* | log | \[1e-4, Inf\] or 0 | log(0) | *beta* = 1/*cv*^2, *sv* = *cv* *v* | Between-trial CV of the input rate |
#' | *accel* | log | \[1e-4, Inf\] or 0 | log(0) | | Exponential acceleration rate of the registration clock |
#'
#' `accel = 0` gives the same likelihood as `parameterization = "cv"` with
#' `clock = 1`.
#'
#' With `plo > 0` each accumulator independently commits at the low threshold
#' `K_lo = 1` with probability `plo` on each trial:
#' `F(u) = (1 - plo) * I_q(u)(K, beta) + plo * I_q(u)(1, beta)`, with
#' `I_q(1, beta) = 1 - (1 - q)^beta`. Both thresholds share `beta`, `p` and the
#' clock. `p` is solved from
#' `h = (1 - plo) * I_p(K, beta) + plo * I_p(1, beta)` by a bracketed root
#' search when `h < 1`.
#'
#' `mapped_pars()` adds the derived columns `beta`, `p`, `lambda`, `Etot`
#' (expected total evidence, `beta * p / (1 - p)`, or `E` at `sv = 0`; `Inf`
#' at `h = 1`) and `cv` = `sv / v` (`sv` = `cv * v` under
#' `parameterization = "cv"` or `"accel"`).
#'
#' Optional nuisance parameters are `pContaminant`, an omission probability
#' outside the accumulators, and `pGuess`, a uniform-outlier probability.
#'
#' @param parameterization Rate-variability and clock coordinates: `"sv"`
#'   (the default) samples the rate SD *sv*; `"cv"` samples the rate CV
#'   *cv* = *sv*/*v*; `"accel"` samples *cv* and the clock acceleration
#'   *accel* in place of *clock*.
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
#' # The cv parameterization.
#' design_FRQ_cv <- design(data = forstmann,
#'                         model = function() FRQ(parameterization = "cv"),
#'                         matchfun = matchfun,
#'                         formula = list(K ~ 1, v ~ lM + E, cv ~ lM, h ~ 1,
#'                                        t0 ~ 1),
#'                         contrasts = list(v = list(lM = ADmat)),
#'                         constants = c(h = qnorm(1)))
#' # The accel parameterization.
#' design_FRQ_accel <- design(data = forstmann,
#'                            model = function() FRQ(parameterization = "accel"),
#'                            matchfun = matchfun,
#'                            formula = list(K ~ 1, v ~ lM + E, cv ~ 1,
#'                                           accel ~ 1, h ~ 1, t0 ~ 1),
#'                            contrasts = list(v = list(lM = ADmat)),
#'                            constants = c(h = qnorm(1)))
#' @export
FRQ <- function(parameterization = c("sv", "cv", "accel")) {
  parameterization <- match.arg(parameterization)
  accel <- parameterization == "accel"
  # The rate-variability coordinate in kernel column 3: the rate SD sv, or the
  # rate CV cv = sv / v (the cv and accel charts).
  rv <- if (accel) "cv" else parameterization
  # Kernel column 7: the clock shape, or the clock acceleration (accel chart).
  ck <- if (accel) "accel" else "clock"
  # K has no explicit transform side: the bound links it, so a sampled K is
  # K = 2 + exp(x) and follows a user bound, while a constant K keeps the
  # plain log link.
  p_types <- c("K" = log(3), "v" = log(6),
               # sv = 0 (log(0)) is the exact no-rate-variability counter.
               "rv" = log(0),
               "h" = qnorm(1), "t0" = log(0),
               # No low-threshold state.
               "plo" = qnorm(0),
               # A constant-rate evidence clock (no acceleration).
               "ck" = if (accel) log(0) else log(1))
  transform <- c(K = "exp", v = "exp", rv = "exp", h = "pnorm", t0 = "exp",
                 plo = "pnorm", ck = "exp")
  # sv (or cv) has no upper bound: sv > v (beta < 1) is a valid, very variable
  # rate, and the kernel solves 1 - p from its own quantile.  Its lower guard
  # keeps a sampled value off the exact zero branch, which is reachable as a
  # constant.
  # The free h upper bound stops short of one so that a sampler can keep a
  # finite defective mass; h = 1 is reachable through the exception.
  # plo's lower guard keeps a sampled value off the exact plo = 0 base model,
  # which is reachable as a constant; likewise clock's and accel's.
  minmax <- cbind(K = c(2, Inf), v = c(1e-4, Inf), rv = c(1e-4, Inf),
                  h = c(1e-6, 1), t0 = c(0.05, Inf), plo = c(1e-6, 1),
                  ck = c(1e-4, Inf))
  exception <- c(t0 = 0, plo = 0, ck = 0, h = 1, rv = 0)
  .rn <- function(x) { x[x == "rv"] <- rv; x[x == "ck"] <- ck; x }
  names(p_types) <- .rn(names(p_types))
  names(transform) <- .rn(names(transform))
  colnames(minmax) <- .rn(colnames(minmax))
  names(exception) <- .rn(names(exception))

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  # pContaminant is redundant with 1 - h unless an omission mechanism outside
  # the decision process is wanted.
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = if (accel) "FRQ_ACCEL" else if (rv == "cv") "FRQ_CV" else "FRQ",
    p_types = p_types,
    p_types_canonical = c("K", "v", rv, "h", "t0"),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: Ttransform does not run on the compiled likelihood
      # path (see src/utils.h); the kernel derives the same quantities itself.
      pl <- if ("plo" %in% colnames(pars)) pars[, "plo"] else rep(0, nrow(pars))
      cu <- if ("clock" %in% colnames(pars)) pars[, "clock"] else rep(1, nrow(pars))
      sv <- .frq_sv(pars)
      g <- frq_rate(pars[, "K"], pars[, "v"], sv, pars[, "h"], pl, cu)
      # Built as a bare matrix first: cbind()ing named scalars onto a one-row
      # `pars` makes R invent row names from the argument names, which would
      # then travel with the parameter matrix.
      # The other rate-variability coordinate.
      other <- if (rv == "cv") cbind(sv = as.numeric(sv)) else
        cbind(cv = as.numeric(sv / pars[, "v"]))
      add <- cbind(beta = as.numeric(g[, "beta"]), p = as.numeric(g[, "p"]),
                   lambda = as.numeric(g[, "lambda"]), Etot = as.numeric(g[, "Etot"]),
                   other)
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
