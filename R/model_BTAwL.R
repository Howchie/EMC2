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
  clear <- if ("Ttrans" %in% colnames(pars)) "Ttrans" else "tau"
  need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", clear)
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL requires parameter columns ", paste(missing, collapse = ", "))
  need
}

.btawl_tau_col <- function(pars) {
  if ("Ttrans" %in% colnames(pars))
    btawl_tau_vec(pars[, "k"], pars[, "Ttrans"])
  else pars[, "tau"]
}

dBTAwL <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau <- .btawl_tau_col(pars)
    out[ok] <- dbtawl(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = tau[ok],
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
    tau <- .btawl_tau_col(pars)
    out[ok] <- pbtawl(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = tau[ok],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

.btawl_mix_check_cols <- function(pars, launch) {
  clear <- if ("Ttrans" %in% colnames(pars)) "Ttrans" else "tau_t"
  need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", clear, "pi")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL mixed requires parameter columns ", paste(missing, collapse = ", "))
  need
}

.btawl_mix_tau_t_col <- function(pars) {
  if ("Ttrans" %in% colnames(pars))
    btawl_tau_vec(pars[, "k"], pars[, "Ttrans"])
  else pars[, "tau_t"]
}

dBTAwLMix <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_mix_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau_t <- .btawl_mix_tau_t_col(pars)
    out[ok] <- dbtawlmix(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = tau_t[ok], pi = pars[ok, "pi"],
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
    tau_t <- .btawl_mix_tau_t_col(pars)
    out[ok] <- pbtawlmix(dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = tau_t[ok], pi = pars[ok, "pi"],
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
  btawl_tmax_vec(k, tau)[1L]
}

.btawl_hit_time_transient <- function(V, z, b, k, tau) {
  if (!is.finite(V) || V <= 0 || !is.finite(z) || z >= b) return(if (z >= b) 0 else Inf)
  X <- function(t) z * exp(-k * t) + V * .btawl_h(t, k, tau)
  if (k <= 1e-10) {
    if (z + V * tau < b) return(Inf)
    hi <- max(tau, 1)
    while (X(hi) < b && hi < 1e12 * max(tau, 1)) hi <- 2 * hi
    return(if (X(hi) < b) Inf else uniroot(function(t) X(t) - b, c(0, hi), tol = 1e-11)$root)
  }
  tm <- .btawl_tmax(k, tau)
  if (is.finite(tm)) {
    if (X(tm) < b) return(Inf)
    return(uniroot(function(t) X(t) - b, c(0, tm), tol = 1e-11)$root)
  }
  hi <- max(tau, 1 / k)
  while (X(hi) < b && hi < 1e12 * max(tau, 1 / k)) hi <- 2 * hi
  if (X(hi) < b) Inf else uniroot(function(t) X(t) - b, c(0, hi), tol = 1e-11)$root
}

.btawl_hit_time_sustained <- function(V, z, b, k, tau_s) {
  if (!is.finite(V) || V <= 0 || !is.finite(z) || z >= b) return(if (z >= b) 0 else Inf)
  X <- function(t) z * exp(-k * t) + V * .btawl_hs(t, k, tau_s)
  if (k > 1e-10 && V <= k * b) return(Inf)
  hi <- max(tau_s, if (k > 1e-10) 1 / k else 1)
  while (X(hi) < b && hi < 1e12 * max(tau_s, 1)) hi <- 2 * hi
  if (X(hi) < b) Inf else uniroot(function(t) X(t) - b, c(0, hi), tol = 1e-11)$root
}

.btawl_hit_time <- function(V, z, b, k, tau,
                            tau_s = tau, tau_t = tau, pi = 0) {
  if (pi <= 1e-14) return(.btawl_hit_time_transient(V, z, b, k, tau_t))
  if (pi >= 1 - 1e-14) return(.btawl_hit_time_sustained(V, z, b, k, tau_s))
  t_T <- .btawl_hit_time_transient(V * (1 - pi), z, b, k, tau_t)
  t_S <- .btawl_hit_time_sustained(V * pi, z, b, k, tau_s)
  min(t_T, t_S)
}

rBTAwL <- function(lR, pars, ok = rep(TRUE, length(lR)),
                   p_types = NULL, posdrift = TRUE, .drifts = NULL,
                   launch = 0L, mixed = FALSE) {
  nm <- if (mixed) {
    clear <- if ("Ttrans" %in% colnames(pars)) "Ttrans" else "tau_t"
    need <- c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", clear, "pi")
    miss <- setdiff(need, colnames(pars))
    if (length(miss)) stop("BTAwL mixed requires parameter columns ", paste(miss, collapse = ", "))
    c(.btawl_par_names(launch), "b", "A", "t0", "k", "tau_s", clear, "pi")
  } else .btawl_check_cols(pars, launch)
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr)
    stop("BTAwL requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  if (is.null(p_types)) {
    # The simulator receives the transformed lower-case threshold `b`; the
    # constructor also carries the pre-transform `B` column, but it is not
    # part of the ballistic kernel contract.
    p_types <- c(nm, "b", "A", "t0", "k")
  }
  if (!all(p_types %in% colnames(pars)))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  tau_all <- if (mixed) NULL else .btawl_tau_col(pars_all)
  tau_t_all <- if (mixed) .btawl_mix_tau_t_col(pars_all) else NULL
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
    if (mixed) {
      pi_val <- pars[j, "pi"]
      if (pi_val <= 1e-14) {
        z <- pars[j, "A"] * runif(1)
        dt[tr, trial] <- .btawl_hit_time_transient(V[j], z, pars[j, "b"], pars[j, "k"],
                                                  tau_t_all[ok_idx[j]])
      } else if (pi_val >= 1 - 1e-14) {
        z <- pars[j, "A"] * runif(1)
        dt[tr, trial] <- .btawl_hit_time_sustained(V[j], z, pars[j, "b"], pars[j, "k"],
                                                  pars[j, "tau_s"])
      } else {
        z_T <- pars[j, "A"] * runif(1)
        z_S <- pars[j, "A"] * runif(1)
        if (!is.null(.drifts)) {
          V_S <- V[j] * pi_val
          V_T <- V[j] * (1 - pi_val)
        } else {
          if (launch == 1L) {
            p1_T <- pars[j, nm[1]] + log(1 - pi_val)
            p1_S <- pars[j, nm[1]] + log(pi_val)
            p2_T <- pars[j, nm[2]]
            p2_S <- pars[j, nm[2]]
            V_T <- rlnorm(1, p1_T, p2_T)
            V_S <- rlnorm(1, p1_S, p2_S)
          } else {
            p1_T <- pars[j, nm[1]] * (1 - pi_val)
            p1_S <- pars[j, nm[1]] * pi_val
            p2_T <- pars[j, nm[2]] * (1 - pi_val)
            p2_S <- pars[j, nm[2]] * pi_val
            if (posdrift) {
              V_T <- msm::rtnorm(1, mean = p1_T, sd = p2_T, lower = 0)
              V_S <- msm::rtnorm(1, mean = p1_S, sd = p2_S, lower = 0)
            } else {
              V_T <- rnorm(1, p1_T, p2_T)
              V_S <- rnorm(1, p1_S, p2_S)
            }
          }
        }
        t_T <- .btawl_hit_time_transient(V_T, z_T, pars[j, "b"], pars[j, "k"],
                                         tau_t_all[ok_idx[j]])
        t_S <- .btawl_hit_time_sustained(V_S, z_S, pars[j, "b"], pars[j, "k"],
                                         pars[j, "tau_s"])
        dt[tr, trial] <- min(t_T, t_S)
      }
    } else {
      z <- pars[j, "A"] * runif(1)
      dt[tr, trial] <- .btawl_hit_time_transient(V[j], z, pars[j, "b"], pars[j, "k"],
                                                tau_all[ok_idx[j]])
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
#' @param mixture If `TRUE`, `"shared"`, or `"mixed"`, return the shared-strength
#'   sustained/transient extension.
#' @param chart Either `"endpoint"` (sample `Ttrans`) or `"rate"` (sample `tau`).
#' @return A BTAwL race-model specification.
#' @export
BTAwL <- function(posdrift = TRUE,
                  drift_distribution = c("normal", "lognormal"),
                  mixture = FALSE,
                  chart = c("endpoint", "rate")) {
  chart <- match.arg(chart)
  if (isTRUE(mixture) || identical(mixture, "shared") || identical(mixture, "mixed"))
    return(BTAwL_mixed(posdrift = posdrift,
                       drift_distribution = drift_distribution,
                       chart = chart))
  drift_distribution <- match.arg(drift_distribution)
  launch <- .btawl_launch_code(drift_distribution)
  lognormal <- launch == 1L
  if (lognormal && !isTRUE(posdrift))
    stop("BTAwL: posdrift only applies to drift_distribution = \"normal\".")
  base_name <- paste0("BTAwL", if (lognormal) "_LOGN" else "",
                      if (chart == "rate") "_RATE" else "",
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
  clear_name <- if (chart == "endpoint") "Ttrans" else "tau"
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0),
               k = log(0), setNames(log(1), clear_name))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp",
                 k = "exp", setNames("exp", clear_name))
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(1e-4, Inf))
  minmax <- cbind(minmax, c(1e-4, Inf))
  colnames(minmax)[ncol(minmax)] <- clear_name
  exception <- if (chart == "endpoint") c(A = 0) else c(A = 0, k = 0)
  launch_pars <- .btawl_par_names(launch)
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  list(
    type = "RACE",
    c_name = base_name,
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      lead <- c(launch_pars, "B", "A", "t0", "k", clear_name)
      extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
      b <- pars[, "B"] + pars[, "A"]
      tau <- if (chart == "endpoint") btawl_tau_vec(pars[, "k"], pars[, "Ttrans"])
             else pars[, "tau"]
      Tmax <- btawl_tmax_vec(pars[, "k"], tau)
      Vcrit <- btawl_vcrit_vec(pars[, "k"], tau, b)
      out <- cbind(pars[, lead, drop = FALSE], extra, b = b, tau = tau,
                   Tmax = Tmax, rt_max = pars[, "t0"] + Tmax, Vcrit = Vcrit)
      if (chart == "endpoint") out[, "Ttrans"] <- Tmax
      out
    },
    rfun = function(data, pars)
      rBTAwL(data$lR, pars, ok = attr(pars, "ok"),
             posdrift = posdrift, launch = launch),
    dfun = function(rt, pars)
      dBTAwL(rt, pars, launch = launch, posdrift = posdrift),
    pfun = function(rt, pars)
      pBTAwL(rt, pars, launch = launch, posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BTAwL: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}

#' Shared-strength sustained/transient BTAwL (within-accumulator race)
#'
#' This is the two-channel extension of [BTAwL()]. The transient and steady-state
#' processes act as independent sub-racers toward the threshold `b`. The
#' launch strengths for the sustained and transient channels are independent draws
#' parameterized by `pi` scaling the location and scale parameters.
#' `pi = 0` recovers the transient-only BTAwL kernel exactly, while `pi = 1`
#' recovers the pure sustained leaky accumulator.
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"` or `"lognormal"`.
#' @param chart Either `"endpoint"` (sample `Ttrans`) or `"rate"` (sample `tau_t`).
#' @return A shared-strength BTAwL race-model specification.
#' @export
BTAwL_mixed <- function(posdrift = TRUE,
                        drift_distribution = c("normal", "lognormal"),
                        chart = c("endpoint", "rate")) {
  chart <- match.arg(chart)
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
  clear_name <- if (chart == "endpoint") "Ttrans" else "tau_t"
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0), k = log(0),
               tau_s = log(1), setNames(log(1), clear_name), pi = qnorm(.5))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                 tau_s = "exp", setNames("exp", clear_name), pi = "pnorm")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = if (chart == "endpoint") c(1e-4, Inf) else c(0, Inf),
                  tau_s = c(1e-4, Inf))
  minmax <- cbind(minmax, c(1e-4, Inf))
  colnames(minmax)[ncol(minmax)] <- clear_name
  minmax <- cbind(minmax, pi = c(0, 1))
  # pi = 0 is deliberately an exception so the nested transient-only member
  # is reachable exactly on the natural scale (and can be used as a constant).
  exception <- if (chart == "endpoint") c(A = 0, pi = 0) else c(A = 0, k = 0, pi = 0)
  launch_pars <- .btawl_par_names(launch)
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  list(
    type = "RACE",
    c_name = paste0("BTAwL_MIX", if (lognormal) "_LOGN" else "",
                    if (chart == "rate") "_RATE" else "",
                    if (!lognormal && !posdrift) "_IO" else ""),
    drift_distribution = drift_distribution,
    shared_strength = TRUE,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      lead <- c(launch_pars, "B", "A", "t0", "k", "tau_s", clear_name, "pi")
      extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
      b <- pars[, "B"] + pars[, "A"]
      tau_t <- if (chart == "endpoint") btawl_tau_vec(pars[, "k"], pars[, "Ttrans"])
               else pars[, "tau_t"]
      pure_transient <- pars[, "pi"] <= 1e-14
      pure_transient[is.na(pure_transient)] <- FALSE
      Tmax <- btawl_tmax_vec(pars[, "k"], tau_t)
      Vcrit <- btawl_vcrit_vec(pars[, "k"], tau_t, b)
      # A sustained component removes the transient hard endpoint.  Keep the
      # transient channel's sampled clearance in the parameter columns, but
      # only expose Tmax/rt_max/Vcrit as endpoint diagnostics for pi = 0.
      Tmax[!pure_transient] <- Inf
      Vcrit[!pure_transient] <- Inf
      out <- cbind(pars[, lead, drop = FALSE], extra, b = b, tau_t = tau_t,
                   Tmax = Tmax, rt_max = pars[, "t0"] + Tmax, Vcrit = Vcrit)
      if (chart == "endpoint") out[pure_transient, "Ttrans"] <- Tmax[pure_transient]
      out
    },
    rfun = function(data, pars)
      rBTAwLMix(data$lR, pars, ok = attr(pars, "ok"),
                posdrift = posdrift, launch = launch),
    dfun = function(rt, pars)
      dBTAwLMix(rt, pars, launch = launch, posdrift = posdrift),
    pfun = function(rt, pars)
      pBTAwLMix(rt, pars, launch = launch, posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BTAwL_mixed: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}

# Short alias matching the c_name suffix, retained for interactive use.
BTAwLMix <- BTAwL_mixed
BTAwLTransient <- function(posdrift = TRUE,
                           drift_distribution = c("normal", "lognormal"))
  BTAwL(posdrift = posdrift, drift_distribution = drift_distribution,
        mixture = FALSE)
