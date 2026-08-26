# ============================================================================
# BAwF: the ballistic accumulator with global fading of decision-relevant
# evidence.
#
#   X(u) = h_rho(u) [z + V u],   h_rho(u) = (1 + k u / rho)^(-rho),
#                                h_Inf(u) = exp(-k u),   z ~ U(0, A)
#
# The latent trace z + V u accumulates ballistically; h_rho(u) attenuates how
# strongly it remains available to the decision mechanism.  For rho = Inf this
# is dX/du = V exp(-k u) - k X, which looks like "drive decay plus leak" but
# has a single primitive rate: the clearance term is k X, so it takes its
# evidence scale from the state rather than carrying an independent one the
# way BAwD's `ell` does.
#
# The consequence, and the reason the model exists, is that the hard right
# endpoint is T_max = x_max / k with x_max fixed by the kernel shape alone
# (1 for rho = Inf, rho/(rho - 1) otherwise).  B cannot move it.  The omission
# boundary V_c(0) = k b H'(x_max) (= e k b for rho = Inf) then identifies b
# given that time scale, so the two dominant observables have distinct jobs.
#
# The only fixed option is rho, which is never estimated.  rho = 1 has no
# finite endpoint and is not a member of this family.
#
# Numerical kernels live in src/model_BAwF.h and are shared by the wrappers
# and the sampled likelihood.  Launch codes must match the BAWD_LAUNCH_*
# constants (BAwF shares them) and the adapter suffixes.
# ============================================================================

# 0 = truncated normal launch (v, sv); 1 = lognormal launch (mu, sigma);
# 2 = median-parameterised continuous split-lognormal launch (mu, sigma,
# delta); 3 = Weibull launch (shape, scale). Must match BAWF_LAUNCH_* in
# src/model_BAwF.h and the value the adapter puts in ctx->bawd_launch from the
# c_name suffix.

# Fixed fading-kernel shape.  rho = 1 is excluded deliberately: V* then falls
# monotonically to k b, so there is no finite endpoint and none of the
# identification architecture this model is built around survives.
.bawf_rho_values <- c(2, 4, Inf)
.bawf_check_rho <- function(rho) {
  if (length(rho) != 1L || is.na(rho))
    stop("BAwF rho must be one of 2, 4, Inf; got ",
         paste(rho, collapse = ", "))
  if (is.infinite(rho)) {
    if (rho < 0) stop("BAwF rho must be one of 2, 4, Inf; got ", rho)
    return(Inf)
  }
  hit <- which(abs(rho - c(2, 4)) < 1e-12)
  if (length(hit) != 1L)
    stop("BAwF rho must be one of 2, 4, Inf (rho = 1 has no finite endpoint ",
         "and is not a member of this family); got ",
         paste(rho, collapse = ", "))
  c(2, 4)[hit]
}
.bawf_rho_suffix <- function(rho) {
  rho <- .bawf_check_rho(rho)
  if (!is.finite(rho)) "" else paste0("_RHO", rho)
}


.bawf_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwF requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwF <- function(rt, pars, launch = 1L, posdrift = TRUE, rho = Inf) {
  rho <- .bawf_check_rho(rho)
  nm <- .bawf_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbawf(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     k = pars[ok, "k"], launch = as.integer(launch),
                     posdrift = posdrift, rho = rho,
                     delta = if (launch == 2L) pars[ok, "delta"] else 0) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBAwF <- function(rt, pars, launch = 1L, posdrift = TRUE, rho = Inf) {
  rho <- .bawf_check_rho(rho)
  nm <- .bawf_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  # rt = Inf is deliberately kept: the CDF there is F_max (the complement of
  # the never-finish mass), not one.  The compiled kernel returns the frozen
  # branch.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbawf(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     k = pars[ok, "k"], launch = as.integer(launch),
                     posdrift = posdrift, rho = rho,
                     delta = if (launch == 2L) pars[ok, "delta"] else 0)
  }
  out
}

# First crossing of X(u) = h_rho(u)[z + V u] through b, or Inf if the peak
# falls short.  Mirrors bawf_hit_time_r() in src/model_rng.h -- including the
# fact that b and z enter separately, because the fading multiplies the start
# point as well as the accumulated evidence.
.bawf_hit_time <- function(V, b, z, k, rho = Inf) {
  rho <- .bawf_check_rho(rho)
  if (!isTRUE(b > z)) return(0)
  if (is.na(V) || !isTRUE(V > 0)) return(Inf)
  if (k <= 1e-10) return((b - z) / V)
  x_max <- if (is.finite(rho)) rho / (rho - 1) else 1
  VK <- V / k
  Hf <- function(x) if (is.finite(rho)) exp(rho * log1p(x / rho)) else exp(x)
  Hp <- function(x) if (is.finite(rho)) exp((rho - 1) * log1p(x / rho)) else exp(x)
  Fx <- function(x) b * Hf(x) - z - VK * x
  # Convexity plus F(0) = b - z > 0 means a crossing exists exactly when
  # F(x_max) <= 0, and then the root lies in (0, x_max].
  if (Fx(x_max) > 0) return(Inf)
  lo <- 0
  hi <- x_max
  x <- 0.5 * x_max
  for (it in seq_len(100)) {
    f <- Fx(x)
    if (f > 0) lo <- x else hi <- x
    d1 <- b * Hp(x) - VK
    xn <- if (is.finite(d1) && d1 < 0) x - f / d1 else 0.5 * (lo + hi)
    if (!isTRUE(xn > lo) || !isTRUE(xn < hi) || !is.finite(xn))
      xn <- 0.5 * (lo + hi)
    done <- abs(xn - x) <= 1e-14 * max(1, xn)
    x <- xn
    if (done) break
  }
  x / k
}

# Pure-R reference simulator; the C++ path uses the identical root solve.
# Retained for the package's explicit emc2.cpp_rfun = FALSE fallback.
rBAwF <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                  posdrift = TRUE, rho = Inf) {
  rho <- .bawf_check_rho(rho)
  nm <- .bawf_check_cols(pars, launch)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  idx <- which(ok)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    V <- if (launch == 3L) {
      rweibull(nrow(p), p[, "shape"], p[, "scale"])
    } else if (launch == 1L) {
      rlnorm(nrow(p), p[, "mu"], p[, "sigma"])
    } else if (launch == 2L) {
      .bawd_split_rlnorm(p[, "mu"], p[, "sigma"], p[, "delta"])
    } else {
      msm::rtnorm(nrow(p), p[, "v"], p[, "sv"],
                  lower = if (posdrift) 0 else -Inf)
    }
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawf_hit_time, V, p[, "b"], z, p[, "k"],
                  MoreArgs = list(rho = rho))
    hit <- .tw_inv(hit, .tw_eta(p))
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  dt <- dt + matrix(t0, nrow = nr)

  bad_col <- colSums(!is.infinite(dt)) == 0L
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, seq_len(ncol(dt)))
  rt <- dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  ok <- matrix(ok, nrow = nr)[1, ]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt[ok] <- rt[ok]
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

#' The Ballistic Accumulator with Global Fading (BAwF)
#'
#' A ballistic race in which the whole evidence trace, start point included,
#' is attenuated by a global temporal gain:
#' \verb{X(u) = h_rho(u) [z + V u]}, with
#' \verb{h_rho(u) = (1 + k u / rho)^(-rho)} and
#' \verb{h_Inf(u) = exp(-k u)}, \verb{z ~ U(0, A)} and \verb{b = B + A}.
#'
#' `X(u) = b` is the response criterion and nothing but `B` and `A` appears in
#' it, so `B` remains a pure caution parameter. The required launch strength
#' \verb{V*(u, z) = (b H_rho(k u) - z) / u}, with
#' \verb{H_rho(x) = (1 + x/rho)^rho}, falls and then rises, so weak launches
#' intrinsically omit and higher start points saturate earlier.
#'
#' Two features distinguish it from [BAwD]. First, the hard right endpoint
#' \verb{T_max = x_max / k} is set by the fading rate and the fixed kernel
#' shape alone (\verb{x_max = 1} for `rho = Inf` and `rho/(rho - 1)`
#' otherwise): `B` cannot counterfeit it. Second, there is no independent
#' clearance scale -- the implied clearance term is `k X`, proportional to the
#' state -- so the `(B, k, ell)` geometry that makes BAwD's endpoint and
#' omission boundary two consequences of one entangled triple does not arise.
#' The omission boundary at the lowest start point is
#' \verb{V_c(0) = k b H'(x_max)}, which is \verb{e k b} for `rho = Inf`;
#' given the time scale, it identifies `b`.
#'
#' The density is closed form for both launch distributions (a truncated
#' partial expectation minus a probability for the live start points, plus a
#' survivor integral over the already-saturated ones). With `A > 0` the
#' density vanishes quadratically at `T_max`, as in BAwD; a point start gives
#' linear shutdown. `k = 0` is the exact LBA limit.
#'
#' `rho` is a model option, never an estimated parameter, so it does not
#' appear in `p_types`. `rho = 1` is refused: it has no finite endpoint.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`.
#' With `"splitlognormal"`, `log V` is continuous split normal with
#' `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`;
#' `mu` is the exact median, the split point is derived from that condition,
#' and `delta` is unbounded.  `delta = 0` reduces exactly to lognormal.
#' With `"weibull"`, `V ~ Weibull(shape, scale)`; both parameters are positive
#' and the BAwF likelihood remains closed form, including its frozen branch.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median/mean log launch location. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean log width. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *scale* | log | \[0, Inf\] | log(1) | | Weibull scale (Weibull launch only). |
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength (normal launch only). |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | *T_max* = *x_max*/*k* | Fading rate. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' For `drift_distribution = "normal"`, use `v` and `sv` instead of the
#' lognormal launch rows; the normal launch is truncated at zero by default.
#' The split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' For `drift_distribution = "weibull"`, use `shape` and `scale` instead of
#' the lognormal launch rows.
#' Optional fitting parameters: `pContaminant` is the omission probability and
#' `pGuess` is the uniform-outlier probability.
#'
#' **Operational-time warp.** `eta` is a trailing free parameter (identity
#' transform, default `0`, unbounded on the natural scale) and is excluded from
#' `p_types_canonical`. For physical accumulation time `u = rt - t0`, let
#' `omega = exp(eta)` and `s = ((1 + u)^omega - 1) / omega`. The CDF and
#' survivor are the parent CDF/survivor evaluated at `s`; the density is the
#' parent density times the Jacobian `c_eta'(u) = (1 + u)^(omega - 1)`.
#' Simulation maps a parent internal finish `s` back with
#' `u = (1 + omega * s)^(1 / omega) - 1` and returns `rt = t0 + u`, so `t0`
#' remains additive. `eta = 0` is the exact identity warp.
#'
#' For BAwL, this contract applies only to the clock-free, uncorrelated
#' constructor (`erlang_type = "none"`, `correlated = FALSE`). BAwL clock or
#' correlated variants, `LogicalRulesLBA`, and all non-ballistic models reject
#' `eta`.
#' @param drift_distribution Distribution of trialwise launch strength:
#'   `"lognormal"` (default), `"splitlognormal"`, `"weibull"`, or `"normal"`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, append `IO` to the compiled model
#'   name. Only meaningful for normal launches.
#' @param rho Fixed fading-kernel shape: `2`, `4`, or `Inf` (default). `Inf`
#'   is the exponential kernel and has no suffix; finite values append
#'   `"_RHO2"` or `"_RHO4"` to the compiled model name. `rho = 1` is not a
#'   member of this family.
#' @return A model list defining the BAwF race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' # mu's intercept is left free and sigma fixes nothing: unlike BAwD there is
#' # no clearance parameter to anchor the evidence scale, so one of mu, B or A
#' # must be pinned.  Here B ~ E carries caution and mu is anchored.
#' design_BAwF <- design(data = forstmann, model = BAwF, matchfun = matchfun,
#'                       formula = list(mu ~ lM, sigma ~ 1, B ~ E, A ~ 1,
#'                                      t0 ~ 1, k ~ 1),
#'                       contrasts = list(mu = list(lM = ADmat)),
#'                       constants = c(mu = 0))
#' @export
BAwF <- function(drift_distribution = c("lognormal", "normal", "splitlognormal", "weibull"),
                 posdrift = TRUE, rho = Inf) {
  drift_distribution <- match.arg(drift_distribution)
  rho <- .bawf_check_rho(rho)

  launch <- .ba_launch_code(drift_distribution, "BAwF")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift)) {
    stop("BAwF: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal and Weibull launch strengths are positive by construction.")
  }

  if (weibull) {
    p_types <- c(shape = log(1), scale = log(1))
    transform <- c(shape = "exp", scale = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), scale = c(1e-4, Inf))
  } else if (lognormal) {
    p_types <- c("mu" = 0, "sigma" = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
    if (splitlognormal) {
      p_types <- c(p_types, delta = 0)
      transform <- c(transform, delta = "identity")
      minmax <- cbind(minmax, delta = c(-Inf, Inf))
    }
  } else {
    p_types <- c("v" = 1, "sv" = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  # Unlike BAwD there is no ell to leave out of the formula, so the evidence
  # scale must be fixed by pinning one of mu (v), B or A -- see the example.
  p_types <- c(p_types, "B" = log(1), "A" = log(0), "t0" = log(0),
               "k" = log(0))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp")
  minmax <- cbind(minmax, A = c(1e-4, Inf), B = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(1e-4, Inf))
  # k = 0 is the exact LBA limit and must stay reachable, so it is a bound
  # exception rather than clamped.
  exception <- c(A = 0, k = 0)
  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  # "_SPLIT" follows "_LOGN"; neither lognormal name contains "IO".
  c_name <- paste0("BAwF", if (weibull) "_WEIB"
                           else if (lognormal) "_LOGN"
                           else if (!posdrift) "IO" else "",
                   if (splitlognormal) "_SPLIT" else "",
                   .bawf_rho_suffix(rho))
  list(
    type = "RACE",
    c_name = c_name,
    drift_distribution = drift_distribution,
    rho = rho,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types),
                                c(.time_warp_par_name, .nuisance_par_names)),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      # Both derived diagnostics are computed by the same kernel the
      # likelihood uses.  Vcrit is the omission boundary at the lowest start
      # point; reporting it alongside Tmax is the point of the model, since
      # the two are what separately identify k and b.
      Tmax_op <- bawf_tmax_vec(pars[, "A"], b, pars[, "k"], rho)
      Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
      Vcrit <- bawf_vcrit_vec(pars[, "A"], b, pars[, "k"], rho)
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax,
            Vcrit = Vcrit)
    },
    rfun = function(data, pars) {
      .rfun_BAwF(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                 posdrift = posdrift, rho = rho)
    },
    dfun = function(rt, pars) dBAwF(rt, pars, launch = launch,
                                    posdrift = posdrift, rho = rho),
    pfun = function(rt, pars) pBAwF(rt, pars, launch = launch,
                                    posdrift = posdrift, rho = rho),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BAwF: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
