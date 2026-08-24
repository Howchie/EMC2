# ============================================================================
# BAwR: the ballistic accumulator with a ramping clearance rate.
#
#   dX/du = V - kappa u^p,   p > 0,   z ~ U(0, A)
#   X(u)  = z + V u - kappa u^(p+1) / (p + 1)
#
# The momentary drive starts at the launch strength V and deteriorates
# deterministically with PHYSICAL time.  There is no state dependence (contrast
# BAwL's leak and BAwF's global fading) and no second opposing process with its
# own evidence scale (contrast BAwD's `ell`): the entire decay law is one
# coefficient and one exponent.  The derivation calls the coefficient `a`; it is
# `kappa` here so that a formula can never confuse it with `A`.
#
# The "R" is for RAMP -- the clearance rate kappa u^p ramps up with elapsed
# time.  It was called BAwP until the exponent `p` made that name collide with
# BAwDp's proportional clearance, which is an unrelated mechanism.  See
# Math/BallisticAccumulators.tex for the derivation and the naming rationale.
#
# Because every trial loses drive at the same ABSOLUTE rate, the trajectory of
# a strong launch peaks later -- at u = (V/kappa)^(1/p) -- so strong sensory
# representations keep supporting accumulation while weak ones expire.
#
# The hard right endpoint is `T_max = [(p+1) b / (kappa p)]^(1/(p+1))`.
# Unlike BAwF's it DEPENDS ON b: decay in physical time knows nothing about
# the threshold, so a more cautious accumulator simply gets longer before
# the drive expires. That is the deliberate trade -- B stays a pure caution
# parameter and kappa a caution-free property of the stimulus representation,
# at the cost of the b-free endpoint.
#
# Numerical kernels live in src/model_BAwR.h and are shared by the wrappers and
# the sampled likelihood.  Launch codes must match the BAWR_LAUNCH_* constants
# (BAwR shares BAwD's) and the adapter suffixes.
# ============================================================================

# 0 = truncated normal launch (v, sv); 1 = lognormal launch (mu, sigma).
# Must match BAWR_LAUNCH_* in src/model_BAwR.h and the value the adapter puts
# in ctx->bawd_launch from the "_LOGN" c_name suffix.
.bawr_launch_code <- function(drift_distribution) {
  switch(drift_distribution, lognormal = 1L, normal = 0L,
         stop("Unknown BAwR drift_distribution: ", drift_distribution))
}

.bawr_par_names <- function(launch) {
  if (launch == 1L) c("mu", "sigma") else c("v", "sv")
}

.bawr_check_cols <- function(pars, launch) {
  need <- c(.bawr_par_names(launch), "b", "A", "t0", "kappa", "p")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwR requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwR <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawr_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbawr(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     kappa = pars[ok, "kappa"], pw = pars[ok, "p"],
                     launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

pBAwR <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawr_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately kept: the CDF there is F_max (the complement of
  # the never-finish mass), not one.  The compiled kernel returns the frozen
  # branch.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pbawr(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     kappa = pars[ok, "kappa"], pw = pars[ok, "p"],
                     launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

# First crossing of X(u) = z + V u - kappa u^(p+1)/(p+1) through b, or Inf if
# the trajectory peak falls short.  Mirrors bawr_hit_time_r() in
# src/model_rng.h.  Takes the distance d = b - z: the decay is a pure function
# of elapsed time, so b and z enter only through their difference (contrast
# .bawf_hit_time, where the fading multiplies the start point).
.bawr_hit_time <- function(V, d, kappa, pw) {
  if (!isTRUE(d > 0)) return(0)
  if (is.na(V) || !isTRUE(V > 0)) return(Inf)
  if (!isTRUE(pw > 0)) return(Inf)
  if (kappa <= 1e-10) return(d / V)
  u_peak <- exp((log(V) - log(kappa)) / pw)
  if (!is.finite(u_peak) || !isTRUE(u_peak > 0)) return(d / V)
  # kappa u_peak^p = V, so the maximum of the accumulated evidence is
  # u_peak V p / (p + 1); below d, the threshold is never reached.
  if (u_peak * V * pw / (pw + 1) < d) return(Inf)
  Gx <- function(u) V * u - kappa * exp((pw + 1) * log(u)) / (pw + 1)
  lo <- 0
  hi <- u_peak
  u <- 0.5 * u_peak
  for (it in seq_len(100)) {
    f <- Gx(u) - d
    if (f < 0) lo <- u else hi <- u
    d1 <- V - kappa * exp(pw * log(u))
    un <- if (is.finite(d1) && d1 > 0) u - f / d1 else 0.5 * (lo + hi)
    if (!isTRUE(un > lo) || !isTRUE(un < hi) || !is.finite(un))
      un <- 0.5 * (lo + hi)
    done <- abs(un - u) <= 1e-14 * max(1, un)
    u <- un
    if (done) break
  }
  u
}

# Pure-R reference simulator; the C++ path uses the identical root solve.
# Retained for the package's explicit emc2.cpp_rfun = FALSE fallback.
rBAwR <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                  posdrift = TRUE) {
  nm <- .bawr_check_cols(pars, launch)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  idx <- which(ok)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    V <- if (launch == 1L) {
      rlnorm(nrow(p), p[, "mu"], p[, "sigma"])
    } else {
      msm::rtnorm(nrow(p), p[, "v"], p[, "sv"],
                  lower = if (posdrift) 0 else -Inf)
    }
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawr_hit_time, V, p[, "b"] - z, p[, "kappa"], p[, "p"])
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  dt <- dt + matrix(t0, nrow = nr)

  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
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

#' The Ballistic Accumulator with a Ramping Clearance Rate (BAwR)
#'
#' A ballistic race in which evidence is cleared at a rate that ramps up with
#' elapsed physical time: \verb{dX/du = V - kappa u^p}, so
#' \verb{X(u) = z + V u - kappa u^(p+1)/(p+1)} with \verb{z ~ U(0, A)} and
#' \verb{b = B + A}.
#'
#' `X(u) = b` is the response criterion and nothing but `B` and `A` appears in
#' it, so `B` remains a pure caution parameter. The required launch strength
#' \verb{V*(u, z) = (b - z)/u + kappa u^p/(p + 1)} is U-shaped, so weak
#' launches intrinsically omit and higher start points saturate earlier.
#'
#' The mechanism is deliberately simpler than [BAwD] or [BAwF]: there is no
#' state dependence at all and no second opposing process, so no independent
#' evidence scale. The whole decay law is the coefficient `kappa` and the
#' exponent `p`. Because every trial loses drive at the same absolute rate, a
#' strong launch stays productive longer -- its trajectory peaks at
#' \verb{(V/kappa)^(1/p)} -- which reads as weak sensory representations
#' expiring quickly while strong ones keep supporting accumulation. `p = 1` is
#' the linear-decay member, where trajectories are parabolas.
#'
#' The hard right endpoint is
#' \verb{T_max = [(p + 1) b / (kappa p)]^(1/(p + 1))} and the omission boundary
#' at the lowest start point is \verb{V_c(0) = kappa T_max^p}. Note that unlike
#' [BAwF] the endpoint here **does** depend on `b`: decay in physical time knows
#' nothing about the threshold, so raising caution simply buys more time before
#' the drive expires. `Tmax`, `rt_max` and `Vcrit` are reported by the model's
#' `Ttransform`.
#'
#' The density is closed form for both launch distributions (a truncated
#' partial expectation minus a probability for the live start points, plus a
#' survivor integral over the already-saturated ones). With `A > 0` the density
#' vanishes quadratically at `T_max`, as in BAwD and BAwF; a point start gives
#' linear shutdown. `kappa = 0` is the exact LBA limit.
#'
#' There is no parameter fixing the evidence scale: \verb{(V, b, A, kappa)}
#' scale jointly with `T_max` and `p` invariant, so one of `mu` (or `v`), `B`
#' or `A` must be pinned -- see the example.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Mean of the log launch strength. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | SD of log launch strength. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *kappa* | log | \[0, Inf\] | log(0) | | Drive-decay coefficient. |
#' | *p* | log | \[0, Inf\] | log(1) | | Drive-decay exponent. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by
#' `v` and `sv`; the launch is truncated positive when `posdrift = TRUE`.
#'
#' @param drift_distribution Distribution of trialwise launch strength:
#'   `"lognormal"` (default) or `"normal"`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, append `IO` to the compiled model
#'   name. Only meaningful for normal launches.
#' @return A model list defining the BAwR race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lM
#' # p is commonly fixed rather than estimated: p = 1 is the linear-decay
#' # member.  mu is pinned because nothing else anchors the evidence scale.
#' design_BAwR <- design(data = forstmann, model = BAwR, matchfun = matchfun,
#'                       formula = list(mu ~ lM, sigma ~ 1, B ~ E, A ~ 1,
#'                                      t0 ~ 1, kappa ~ 1, p ~ 1),
#'                       contrasts = list(mu = list(lM = ADmat)),
#'                       constants = c(mu = 0, p = log(1)))
#' @export
BAwR <- function(drift_distribution = c("lognormal", "normal"),
                 posdrift = TRUE) {
  drift_distribution <- match.arg(drift_distribution)

  launch <- .bawr_launch_code(drift_distribution)
  lognormal <- (launch == 1L)
  if (lognormal && !isTRUE(posdrift)) {
    stop("BAwR: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal launch strength is positive by construction.")
  }

  if (lognormal) {
    p_types <- c("mu" = 0, "sigma" = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
  } else {
    p_types <- c("v" = 1, "sv" = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  # As for BAwF there is no clearance parameter to leave out of the formula, so
  # the evidence scale must be fixed by pinning one of mu (v), B or A.
  p_types <- c(p_types, "B" = log(1), "A" = log(0), "t0" = log(0),
               "kappa" = log(0), "p" = log(1))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", kappa = "exp",
                 p = "exp")
  minmax <- cbind(minmax, A = c(1e-4, Inf), B = c(1e-4, Inf),
                  t0 = c(0.05, Inf), kappa = c(1e-4, Inf), p = c(1e-3, Inf))
  # kappa = 0 is the exact LBA limit and must stay reachable, so it is a bound
  # exception rather than clamped. p has no such limit: p = 0 would make the
  # drive lose a constant, degenerating to an LBA with a shifted launch and an
  # infinite T_max, so its lower bound is enforced.
  exception <- c(A = 0, kappa = 0)
  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  # "_LOGN" (not "_LN": resolve_race_model_adapter dispatches by substring and
  # "LNR" is an existing key).  The IO suffix is only reachable for the normal
  # launch, and "BAwR_LOGN" deliberately contains no "IO".
  c_name <- paste0("BAwR", if (lognormal) "_LOGN"
                           else if (!posdrift) "IO" else "")
  list(
    type = "RACE",
    c_name = c_name,
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      # Both derived diagnostics are computed by the same kernel the likelihood
      # uses. Tmax is reported because it is not a sampled quantity here (it
      # is a function of kappa, p AND b), and Vcrit because it sets the
      # omission rate.
      Tmax <- bawr_tmax_vec(pars[, "A"], b, pars[, "kappa"], pars[, "p"])
      Vcrit <- bawr_vcrit_vec(pars[, "A"], b, pars[, "kappa"], pars[, "p"])
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax,
            Vcrit = Vcrit)
    },
    rfun = function(data, pars) {
      .rfun_BAwR(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                 posdrift = posdrift)
    },
    dfun = function(rt, pars) dBAwR(rt, pars, launch = launch,
                                    posdrift = posdrift),
    pfun = function(rt, pars) pBAwR(rt, pars, launch = launch,
                                    posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BAwR: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
