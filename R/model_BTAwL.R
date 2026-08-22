# ============================================================================
# BTAwL: Ballistic Transient Accumulator with Leak
#
#   dX/du = V (u/tau) exp(-u/tau) - k X,   X(0) = z, z ~ U(0, A)
#
# The compiled scalar/vector kernels in src/model_BTAwL.h evaluate the exact
# transient solution.  The launch distribution can be a (positive-truncated)
# normal or a lognormal, and the constructor follows the BAwD/BAwL race-model
# contract so it can be used directly from design().
# ============================================================================

.btawl_launch_code <- function(drift_distribution) {
  switch(drift_distribution, normal = 0L, lognormal = 1L,
         stop("Unknown BTAwL drift_distribution: ", drift_distribution))
}

.btawl_par_names <- function(launch) {
  if (launch == 1L) c("mu", "sigma") else c("v", "sv")
}

.btawl_check_cols <- function(pars, launch) {
  need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBTAwL <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbtawl(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = pars[ok, "tau"],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

pBTAwL <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately retained: pBTAwL(Inf) is the eventual response
  # probability, not one, because weak launches can fail to hit.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pbtawl(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = pars[ok, "tau"],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

.btawl_mix_check_cols <- function(pars, launch) {
  need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", "tau_t", "pi")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL mixed requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBTAwLMix <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_mix_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbtawlmix(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = pars[ok, "tau_t"], pi = pars[ok, "pi"],
                         launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

pBTAwLMix <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_mix_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pbtawlmix(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = pars[ok, "tau_t"], pi = pars[ok, "pi"],
                         launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

.btawl_h <- function(t, k, tau) {
  ifelse(k <= 1e-10,
         tau * (1 - (1 + t / tau) * exp(-t / tau)),
         ifelse(abs(k - 1 / tau) < 1e-8 / pmax(1, tau),
                0.5 * t^2 / tau * exp(-t / tau),
                (exp(-t / tau) * ((k - 1 / tau) * t - 1) + exp(-k * t)) /
                  (tau * (k - 1 / tau)^2)))
}

.btawl_hs <- function(t, k, tau) {
  ifelse(k <= 1e-10,
         t - tau * (1 - exp(-t / tau)),
         ifelse(abs(k - 1 / tau) < 1e-8 / pmax(1, tau),
                (1 - exp(-k * t)) / k - t * exp(-k * t),
                (1 - exp(-k * t)) / k -
                  (exp(-t / tau) - exp(-k * t)) / (k - 1 / tau)))
}

.btawl_hmix <- function(t, k, tau_s, tau_t, pi) {
  pi * .btawl_hs(t, k, tau_s) + (1 - pi) * .btawl_h(t, k, tau_t)
}

.btawl_tmax_mix <- function(k, tau_s, tau_t, pi) {
  if (!is.finite(k) || k <= 1e-10) return(Inf)
  hp <- function(t) pi * (1 - exp(-t / tau_s)) +
    (1 - pi) * (t / tau_t) * exp(-t / tau_t) -
    k * .btawl_hmix(t, k, tau_s, tau_t, pi)
  hi <- 2 * max(tau_s, tau_t, 1 / k)
  if (!isTRUE(hp(hi) < 0)) return(Inf)
  lo <- max(.Machine$double.eps * max(1, hi), 1e-12 * min(tau_s, tau_t, 1 / k))
  uniroot(hp, c(lo, hi), tol = 1e-12)$root
}

.btawl_tmax <- function(k, tau) {
  if (!is.finite(k) || !is.finite(tau) || k <= 1e-10 || tau <= 0) return(Inf)
  hp <- function(t) (t / tau) * exp(-t / tau) - k * .btawl_h(t, k, tau)
  hi <- 2 * max(tau, 1 / k)
  while (isTRUE(hp(hi) > 0) && hi < 1e12 * max(tau, 1 / k)) hi <- 2 * hi
  if (!isTRUE(hp(hi) <= 0)) return(Inf)
  lo <- max(.Machine$double.eps * max(1, hi), 1e-12 * min(tau, 1 / k))
  uniroot(hp, c(lo, hi), tol = 1e-12)$root
}

.btawl_hit_time <- function(V, z, b, k, tau,
                            tau_s = tau, tau_t = tau, pi = 0) {
  if (!is.finite(V) || V <= 0 || !is.finite(z) || z >= b) return(if (z >= b) 0 else Inf)
  H <- function(t) .btawl_hmix(t, k, tau_s, tau_t, pi)
  X <- function(t) z * exp(-k * t) + V * H(t)
  if (k <= 1e-10) {
    if (V <= 0 || (pi <= 1e-12 && z + V * tau_t < b)) return(Inf)
    hi <- max(tau_s, tau_t, 1)
    while (X(hi) < b && hi < 1e12 * max(tau_s, tau_t, 1)) hi <- 2 * hi
    return(if (X(hi) < b) Inf else uniroot(function(t) X(t) - b,
                                             c(0, hi), tol = 1e-11)$root)
  }
  tm <- .btawl_tmax_mix(k, tau_s, tau_t, pi)
  if (is.finite(tm)) {
    if (X(tm) < b) return(Inf)
    return(uniroot(function(t) X(t) - b, c(0, tm), tol = 1e-11)$root)
  }
  hi <- max(tau_s, tau_t, 1 / k)
  while (X(hi) < b && hi < 1e12 * max(tau_s, tau_t, 1 / k)) hi <- 2 * hi
  if (X(hi) < b) Inf else uniroot(function(t) X(t) - b, c(0, hi), tol = 1e-11)$root
}

rBTAwL <- function(lR, pars, ok = rep(TRUE, length(lR)),
                   p_types = NULL, posdrift = TRUE, .drifts = NULL,
                   launch = 0L, mixed = FALSE) {
  nm <- if (mixed) {
    need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", "tau_t", "pi")
    miss <- setdiff(need, colnames(pars))
    if (length(miss)) stop("BTAwL mixed requires parameter columns ", paste(miss, collapse = ", "))
    c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", "tau_t", "pi")
  } else .btawl_check_cols(pars, launch)
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr)
    stop("BTAwL requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  if (is.null(p_types)) {
    # The simulator receives the transformed lower-case threshold `b`; the
    # constructor also carries the pre-transform `B` column, but it is not
    # part of the ballistic kernel contract.
    p_types <- if (mixed) c(nm, "b", "A", "t0", "k", "tau_s", "tau_t", "pi") else
      c(nm, "b", "A", "t0", "k", "tau")
  }
  if (!all(p_types %in% colnames(pars)))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  if (is.null(.drifts)) {
    V <- if (launch == 1L) rlnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]]) else
      if (posdrift) msm::rtnorm(nrow(pars), mean = pars[, nm[1]], sd = pars[, nm[2]], lower = 0) else
        rnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  } else {
    if (length(.drifts) != nrow(pars_all))
      stop(".drifts must have one value per row of the original parameter matrix.")
    V <- .drifts[ok]
  }
  dt <- matrix(Inf, nr, n_trials)
  for (j in seq_len(nrow(pars))) {
    row <- ok_idx[j]
    tr <- ((row - 1L) %% nr) + 1L
    trial <- ((row - 1L) %/% nr) + 1L
    z <- pars[j, "A"] * runif(1)
    if (mixed) {
      dt[tr, trial] <- .btawl_hit_time(V[j], z, pars[j, "b"], pars[j, "k"],
                                        pars[j, "tau_t"], pars[j, "tau_s"],
                                        pars[j, "tau_t"], pars[j, "pi"])
    } else {
      dt[tr, trial] <- .btawl_hit_time(V[j], z, pars[j, "b"], pars[j, "k"], pars[j, "tau"])
    }
    dt[tr, trial] <- dt[tr, trial] + pars[j, "t0"]
  }
  bad <- apply(dt, 2, function(x) all(!is.finite(x)))
  # Each column is one trial and each row one racer; transpose before
  # max.col so the winning racer is selected independently per trial.
  win <- max.col(-t(dt), ties.method = "first")
  rt <- dt[cbind(win, seq_len(n_trials))]
  R <- factor(levels(lR)[win], levels = levels(lR))
  R[bad] <- NA; rt[bad] <- Inf
  # Trials with no active parameter rows are omissions.
  active <- matrix(ok, nrow = nr)[1, ]
  R[!active] <- NA; rt[!active] <- Inf
  data.frame(R = R, rt = rt)
}

rBTAwLMix <- function(lR, pars, ok = rep(TRUE, length(lR)),
                      posdrift = TRUE, .drifts = NULL, launch = 0L) {
  rBTAwL(lR, pars, ok = ok, posdrift = posdrift, .drifts = .drifts,
         launch = launch, mixed = TRUE)
}

#' The Ballistic Transient Accumulator with Leak (BTAwL)
#'
#' BTAwL keeps Smith's onset-evoked transient drive while replacing the
#' within-trial Wiener noise by trialwise launch variability:
#' `dX/du = V (u/tau) exp(-u/tau) - k X`, with `X(0) ~ U(0,A)`.
#' For `k > 0` the transient solution has a finite intrinsic response window;
#' for `k = 0` it approaches a finite asymptote and has no hard endpoint.
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"` (`V ~ N(v, sv^2)`) or
#'   `"lognormal"` (`log V ~ N(mu, sigma^2)`).
#' @return A BTAwL race-model specification.
#' @export
BTAwL <- function(posdrift = TRUE,
                  drift_distribution = c("normal", "lognormal"),
                  mixture = FALSE) {
  if (isTRUE(mixture) || identical(mixture, "shared") || identical(mixture, "mixed"))
    return(BTAwL_mixed(posdrift = posdrift,
                       drift_distribution = drift_distribution))
  drift_distribution <- match.arg(drift_distribution)
  launch <- .btawl_launch_code(drift_distribution)
  lognormal <- launch == 1L
  if (lognormal && !isTRUE(posdrift))
    stop("BTAwL: posdrift only applies to drift_distribution = \"normal\".")
  base_name <- paste0("BTAwL", if (lognormal) "_LOGN" else "",
                      if (!lognormal && !posdrift) "_IO" else "")
  if (lognormal) {
    p_types <- c(mu = 0, sigma = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
  } else {
    p_types <- c(v = 1, sv = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0),
               k = log(0), tau = log(1))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp",
                 k = "exp", tau = "exp")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(0, Inf), tau = c(1e-4, Inf))
  exception <- c(A = 0, k = 0)
  launch_pars <- .btawl_par_names(launch)
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  list(
    type = "RACE",
    c_name = base_name,
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = c(launch_pars, "B", "A", "t0", "k", "tau"),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      lead <- c(launch_pars, "B", "A", "t0", "k", "tau")
      extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
      cbind(pars[, lead, drop = FALSE], extra,
            b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data, pars)
      rBTAwL(data$lR, pars, ok = attr(pars, "ok"),
             posdrift = posdrift, launch = launch),
    dfun = function(rt, pars)
      dBTAwL(rt, pars, launch = launch, posdrift = posdrift),
    pfun = function(rt, pars)
      pBTAwL(rt, pars, launch = launch, posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm,
                                   model = model, min_ll = min_ll)
  )
}

#' Shared-strength sustained/transient BTAwL
#'
#' This is the nested two-channel extension of [BTAwL()].  One trialwise
#' launch strength multiplies both Smith channels:
#' `V_S = pi * V`, `V_T = (1 - pi) * V`.  `pi = 0` is exactly the transient-only
#' BTAwL kernel (with `tau_t` playing the old `tau` role), so the original model
#' is a literal submodel rather than a limiting approximation.
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"` or `"lognormal"`.
#' @return A shared-strength BTAwL race-model specification.
#' @export
BTAwL_mixed <- function(posdrift = TRUE,
                        drift_distribution = c("normal", "lognormal")) {
  drift_distribution <- match.arg(drift_distribution)
  launch <- .btawl_launch_code(drift_distribution)
  lognormal <- launch == 1L
  if (lognormal && !isTRUE(posdrift))
    stop("BTAwL_mixed: posdrift only applies to drift_distribution = \"normal\".")
  if (lognormal) {
    p_types <- c(mu = 0, sigma = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
  } else {
    p_types <- c(v = 1, sv = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0), k = log(0),
               tau_s = log(1), tau_t = log(1), pi = qnorm(.5))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                 tau_s = "exp", tau_t = "exp", pi = "pnorm")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(0, Inf),
                  tau_s = c(1e-4, Inf), tau_t = c(1e-4, Inf), pi = c(0, 1))
  # pi = 0 is deliberately an exception so the nested transient-only member
  # is reachable exactly on the natural scale (and can be used as a constant).
  exception <- c(A = 0, k = 0, pi = 0)
  launch_pars <- .btawl_par_names(launch)
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  list(
    type = "RACE",
    c_name = paste0("BTAwL_MIX", if (lognormal) "_LOGN" else "",
                    if (!lognormal && !posdrift) "_IO" else ""),
    drift_distribution = drift_distribution,
    shared_strength = TRUE,
    p_types = p_types,
    p_types_canonical = c(launch_pars, "B", "A", "t0", "k", "tau_s", "tau_t", "pi"),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      lead <- c(launch_pars, "B", "A", "t0", "k", "tau_s", "tau_t", "pi")
      extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
      cbind(pars[, lead, drop = FALSE], extra,
            b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data, pars)
      rBTAwLMix(data$lR, pars, ok = attr(pars, "ok"),
                posdrift = posdrift, launch = launch),
    dfun = function(rt, pars)
      dBTAwLMix(rt, pars, launch = launch, posdrift = posdrift),
    pfun = function(rt, pars)
      pBTAwLMix(rt, pars, launch = launch, posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm,
                                   model = model, min_ll = min_ll)
  )
}

# Short alias matching the c_name suffix, retained for interactive use.
BTAwLMix <- BTAwL_mixed
BTAwLTransient <- function(posdrift = TRUE,
                           drift_distribution = c("normal", "lognormal"))
  BTAwL(posdrift = posdrift, drift_distribution = drift_distribution,
        mixture = FALSE)
