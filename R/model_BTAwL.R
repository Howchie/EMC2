# ============================================================================
# BTAwL: Ballistic transient/sustained local race with leak
#
#   transient: dX_T/du = V_T (u/tau_t) exp(-u/tau_t) - k X_T
#   sustained: dX_S/du = V_S (1 - exp(-u/tau_s)) - k X_S
#
# The compiled scalar/vector kernels in src/model_BTAwL.h evaluate the exact
# transient and sustained local-race solutions.  The launch distribution can
# be a (positive-truncated) normal, lognormal, or Weibull, and the constructor follows
# the BAwD/BAwL race-model contract so it can be used directly from design().
# ============================================================================



.btawl_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "tau")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL requires parameter columns ", paste(missing, collapse = ", "))
  need
}

.btawl_tau_col <- function(pars) {
  pars[, "tau"]
}

dBTAwLTransient <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau <- .btawl_tau_col(pars)
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbtawl_transient(s, A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = tau[ok],
                      launch = as.integer(launch), posdrift = posdrift,
                      delta = if (launch == 2L) pars[ok, "delta"] else 0) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBTAwLTransient <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  # rt = Inf is deliberately retained: pBTAwLTransient(Inf) is the eventual response
  # probability, not one, because weak launches can fail to hit.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau <- .btawl_tau_col(pars)
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbtawl_transient(s, A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], tau = tau[ok],
                      launch = as.integer(launch), posdrift = posdrift,
                      delta = if (launch == 2L) pars[ok, "delta"] else 0)
  }
  out
}

.btawl_local_race_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "tau_s", "tau_t", "pi")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwL requires parameter columns ", paste(missing, collapse = ", "))
  need
}

.btawl_local_race_tau_t_col <- function(pars) {
  pars[, "tau_t"]
}

dBTAwL <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_local_race_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau_t <- .btawl_local_race_tau_t_col(pars)
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbtawl_local_race(s, A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = tau_t[ok], pi = pars[ok, "pi"],
                         launch = as.integer(launch), posdrift = posdrift,
                         delta = if (launch == 2L) pars[ok, "delta"] else 0) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBTAwL <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_local_race_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    tau_t <- .btawl_local_race_tau_t_col(pars)
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbtawl_local_race(s, A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = tau_t[ok], pi = pars[ok, "pi"],
                         launch = as.integer(launch), posdrift = posdrift,
                         delta = if (launch == 2L) pars[ok, "delta"] else 0)
  }
  out
}

.btawl_separate_par_names <- function(launch) {
  switch(as.character(launch),
         `0` = c("v_S", "sv_S", "v_T", "sv_T"),
         `1` = c("mu_S", "sigma_S", "mu_T", "sigma_T"),
         `2` = c("mu_S", "sigma_S", "delta_S", "mu_T", "sigma_T", "delta_T"),
         `3` = c("shape_S", "mean_S", "shape_T", "mean_T"),
         stop("invalid BTAwL launch code"))
}

.btawl_separate_check_cols <- function(pars, launch) {
  need <- c(.btawl_separate_par_names(launch),
            "b", "A", "t0", "k", "tau_s", "tau_t")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwLSeparate requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBTAwLSeparate <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_separate_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    get <- function(nm) if (nm %in% colnames(pars)) pars[ok, nm] else 0
    out[ok] <- dbtawl_local_race_separate(
      s, A = pars[ok, "A"], b = pars[ok, "b"],
      p1_S = pars[ok, nm[1]], p2_S = pars[ok, nm[2]],
      p1_T = pars[ok, nm[if (launch == 2L) 4L else 3L]],
      p2_T = pars[ok, nm[if (launch == 2L) 5L else 4L]],
      k = pars[ok, "k"], tau_s = pars[ok, "tau_s"], tau_t = pars[ok, "tau_t"],
      launch = as.integer(launch), posdrift = posdrift,
      delta_S = get("delta_S"), delta_T = get("delta_T")) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBTAwLSeparate <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_separate_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    get <- function(nm) if (nm %in% colnames(pars)) pars[ok, nm] else 0
    out[ok] <- pbtawl_local_race_separate(
      s, A = pars[ok, "A"], b = pars[ok, "b"],
      p1_S = pars[ok, nm[1]], p2_S = pars[ok, nm[2]],
      p1_T = pars[ok, nm[if (launch == 2L) 4L else 3L]],
      p2_T = pars[ok, nm[if (launch == 2L) 5L else 4L]],
      k = pars[ok, "k"], tau_s = pars[ok, "tau_s"], tau_t = pars[ok, "tau_t"],
      launch = as.integer(launch), posdrift = posdrift,
      delta_S = get("delta_S"), delta_T = get("delta_T"))
  }
  out
}

.btawl_sustained_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "tau_s")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BTAwLSustained requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBTAwLSustained <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_sustained_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    # The transient clearance and pi are inert at the sustained boundary;
    # using the full kernel here keeps the R reference exactly aligned with
    # the compiled sustained adapter.
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbtawl_local_race(s, A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = pars[ok, "tau_s"], pi = 1,
                         launch = as.integer(launch), posdrift = posdrift,
                         delta = if (launch == 2L) pars[ok, "delta"] else 0) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBTAwLSustained <- function(rt, pars, launch = 0L, posdrift = TRUE) {
  nm <- .btawl_sustained_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbtawl_local_race(s, A = pars[ok, "A"], b = pars[ok, "b"],
                         p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                         k = pars[ok, "k"], tau_s = pars[ok, "tau_s"],
                         tau_t = pars[ok, "tau_s"], pi = 1,
                         launch = as.integer(launch), posdrift = posdrift,
                         delta = if (launch == 2L) pars[ok, "delta"] else 0)
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

.rBTAwLTransient_R <- function(lR, pars, ok = rep(TRUE, length(lR)),
                               p_types = NULL, posdrift = TRUE, launch = 0L) {
  nm <- .btawl_check_cols(pars, launch)
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
  tau_all <- .btawl_tau_col(pars_all)
  eta_all <- .tw_eta(pars_all)
  V <- if (launch == 3L) {
    .rweibull_mean(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  } else if (launch == 1L) {
    rlnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  } else if (launch == 2L) {
    .bawd_split_rlnorm(pars[, "mu"], pars[, "sigma"], pars[, "delta"])
  } else if (posdrift) {
    msm::rtnorm(nrow(pars), mean = pars[, nm[1]], sd = pars[, nm[2]], lower = 0)
  } else {
    rnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  }
  dt <- matrix(Inf, nr, n_trials)
  for (j in seq_len(nrow(pars))) {
    row <- ok_idx[j]
    tr <- ((row - 1L) %% nr) + 1L
    trial <- ((row - 1L) %/% nr) + 1L
    z <- pars[j, "A"] * runif(1)
    hit <- .btawl_hit_time_transient(V[j], z, pars[j, "b"], pars[j, "k"],
                                     tau_all[ok_idx[j]])
    hit <- .tw_inv(hit, eta_all[ok_idx[j]])
    dt[tr, trial] <- hit + pars[j, "t0"]
  }
  bad <- colSums(is.finite(dt)) == 0L
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

.rBTAwL_R <- function(lR, pars, ok = rep(TRUE, length(lR)),
                      p_types = NULL, posdrift = TRUE, launch = 0L) {
  nm <- .btawl_local_race_check_cols(pars, launch)
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr)
    stop("BTAwL requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  if (is.null(p_types))
    p_types <- c(nm, "b", "A", "t0", "k")
  if (!all(p_types %in% colnames(pars)))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  tau_t_all <- .btawl_local_race_tau_t_col(pars_all)
  eta_all <- .tw_eta(pars_all)
  dt <- matrix(Inf, nr, n_trials)
  draw_launch <- function(p1, p2, delta = 0) {
    if (launch == 3L) return(.rweibull_mean(1, p1, p2))
    if (launch == 1L) return(rlnorm(1, p1, p2))
    if (launch == 2L) return(.bawd_split_rlnorm(p1, p2, delta))
    if (posdrift) msm::rtnorm(1, mean = p1, sd = p2, lower = 0) else
      rnorm(1, p1, p2)
  }
  for (j in seq_len(nrow(pars))) {
    row <- ok_idx[j]
    tr <- ((row - 1L) %% nr) + 1L
    trial <- ((row - 1L) %/% nr) + 1L
    pi_val <- pars[j, "pi"]
    # The likelihood integrates independent channel launches.  The simulator
    # therefore draws them independently from the same constrained scaling:
    # lognormal locations shift by log(pi), while normal locations and scales
    # are both scaled by pi; split widths are unchanged by this shift.
    if (pi_val <= 1e-14) {
      p1_T <- pars[j, nm[1]]; p2_T <- pars[j, nm[2]]
      V_T <- draw_launch(p1_T, p2_T,
                         if (launch == 2L) pars[j, "delta"] else 0)
      V_S <- 0
    } else if (pi_val >= 1 - 1e-14) {
      p1_S <- pars[j, nm[1]]; p2_S <- pars[j, nm[2]]
      V_T <- 0
      V_S <- draw_launch(p1_S, p2_S,
                         if (launch == 2L) pars[j, "delta"] else 0)
    } else if (launch == 1L || launch == 2L) {
      p1_T <- pars[j, nm[1]] + log1p(-pi_val)
      p1_S <- pars[j, nm[1]] + log(pi_val)
      p2_T <- p2_S <- pars[j, nm[2]]
    } else if (launch == 3L) {
      p1_T <- p1_S <- pars[j, nm[1]]
      p2_T <- pars[j, nm[2]] * (1 - pi_val)
      p2_S <- pars[j, nm[2]] * pi_val
    } else {
      p1_T <- pars[j, nm[1]] * (1 - pi_val)
      p1_S <- pars[j, nm[1]] * pi_val
      p2_T <- pars[j, nm[2]] * (1 - pi_val)
      p2_S <- pars[j, nm[2]] * pi_val
    }
    z_T <- pars[j, "A"] * runif(1)
    z_S <- pars[j, "A"] * runif(1)
    if (pi_val > 1e-14 && pi_val < 1 - 1e-14) {
      V_T <- draw_launch(p1_T, p2_T,
                         if (launch == 2L) pars[j, "delta"] else 0)
      V_S <- draw_launch(p1_S, p2_S,
                         if (launch == 2L) pars[j, "delta"] else 0)
    }
    t_T <- .btawl_hit_time_transient(V_T, z_T, pars[j, "b"], pars[j, "k"],
                                     tau_t_all[ok_idx[j]])
    t_S <- .btawl_hit_time_sustained(V_S, z_S, pars[j, "b"], pars[j, "k"],
                                     pars[j, "tau_s"])
    hit <- .tw_inv(min(t_T, t_S), eta_all[ok_idx[j]])
    dt[tr, trial] <- hit + pars[j, "t0"]
  }
  bad <- colSums(is.finite(dt)) == 0L
  win <- max.col(-t(dt), ties.method = "first")
  rt <- dt[cbind(win, seq_len(n_trials))]
  R <- factor(levels(lR)[win], levels = levels(lR))
  R[bad] <- NA; rt[bad] <- Inf
  active <- matrix(ok, nrow = nr)[1, ]
  R[!active] <- NA; rt[!active] <- Inf
  data.frame(R = R, rt = rt)
}

.rBTAwLSeparate_R <- function(lR, pars, ok = rep(TRUE, length(lR)),
                              p_types = NULL, posdrift = TRUE, launch = 0L) {
  nm <- .btawl_separate_check_cols(pars, launch)
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr)
    stop("BTAwLSeparate requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  if (is.null(p_types)) p_types <- c(nm, "b", "A", "t0", "k")
  if (!all(p_types %in% colnames(pars)))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  eta_all <- .tw_eta(pars_all)
  dt <- matrix(Inf, nr, n_trials)
  draw_launch <- function(p1, p2, delta = 0) {
    if (launch == 3L) return(.rweibull_mean(1, p1, p2))
    if (launch == 1L) return(rlnorm(1, p1, p2))
    if (launch == 2L) return(.bawd_split_rlnorm(p1, p2, delta))
    if (posdrift) msm::rtnorm(1, mean = p1, sd = p2, lower = 0) else
      rnorm(1, p1, p2)
  }
  for (j in seq_len(nrow(pars))) {
    row <- ok_idx[j]
    tr <- ((row - 1L) %% nr) + 1L
    trial <- ((row - 1L) %/% nr) + 1L
    if (launch == 2L) {
      iS <- c("mu_S", "sigma_S", "delta_S")
      iT <- c("mu_T", "sigma_T", "delta_T")
    } else if (launch == 3L) {
      iS <- c("shape_S", "mean_S")
      iT <- c("shape_T", "mean_T")
    } else if (launch == 1L) {
      iS <- c("mu_S", "sigma_S")
      iT <- c("mu_T", "sigma_T")
    } else {
      iS <- c("v_S", "sv_S")
      iT <- c("v_T", "sv_T")
    }
    V_S <- draw_launch(pars[j, iS[1]], pars[j, iS[2]],
                       if (launch == 2L) pars[j, iS[3]] else 0)
    V_T <- draw_launch(pars[j, iT[1]], pars[j, iT[2]],
                       if (launch == 2L) pars[j, iT[3]] else 0)
    z_S <- pars[j, "A"] * runif(1)
    z_T <- pars[j, "A"] * runif(1)
    t_S <- .btawl_hit_time_sustained(V_S, z_S, pars[j, "b"], pars[j, "k"],
                                     pars[j, "tau_s"])
    t_T <- .btawl_hit_time_transient(V_T, z_T, pars[j, "b"], pars[j, "k"],
                                     pars[j, "tau_t"])
    hit <- .tw_inv(min(t_S, t_T), eta_all[row])
    dt[tr, trial] <- hit + pars[j, "t0"]
  }
  bad <- colSums(is.finite(dt)) == 0L
  win <- max.col(-t(dt), ties.method = "first")
  rt <- dt[cbind(win, seq_len(n_trials))]
  R <- factor(levels(lR)[win], levels = levels(lR))
  R[bad] <- NA; rt[bad] <- Inf
  active <- matrix(ok, nrow = nr)[1, ]
  R[!active] <- NA; rt[!active] <- Inf
  data.frame(R = R, rt = rt)
}

.rBTAwLSustained_R <- function(lR, pars, ok = rep(TRUE, length(lR)),
                               p_types = NULL, posdrift = TRUE, launch = 0L) {
  nm <- .btawl_sustained_check_cols(pars, launch)
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr)
    stop("BTAwL requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  if (is.null(p_types)) p_types <- c(nm, "b", "A", "t0", "k")
  if (!all(p_types %in% colnames(pars)))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  eta_vec <- .tw_eta(pars)
  V <- if (launch == 3L) {
    .rweibull_mean(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  } else if (launch == 1L) {
    rlnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  } else if (launch == 2L) {
    .bawd_split_rlnorm(pars[, "mu"], pars[, "sigma"], pars[, "delta"])
  } else if (posdrift) {
    msm::rtnorm(nrow(pars), mean = pars[, nm[1]], sd = pars[, nm[2]], lower = 0)
  } else {
    rnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
  }
  dt <- matrix(Inf, nr, n_trials)
  for (j in seq_len(nrow(pars))) {
    row <- ok_idx[j]
    tr <- ((row - 1L) %% nr) + 1L
    trial <- ((row - 1L) %/% nr) + 1L
    z <- pars[j, "A"] * runif(1)
    hit <- .btawl_hit_time_sustained(V[j], z, pars[j, "b"], pars[j, "k"],
                                     pars[j, "tau_s"])
    hit <- .tw_inv(hit, eta_vec[j])
    dt[tr, trial] <- hit + pars[j, "t0"]
  }
  bad <- colSums(is.finite(dt)) == 0L
  win <- max.col(-t(dt), ties.method = "first")
  rt <- dt[cbind(win, seq_len(n_trials))]
  R <- factor(levels(lR)[win], levels = levels(lR))
  R[bad] <- NA; rt[bad] <- Inf
  active <- matrix(ok, nrow = nr)[1, ]
  R[!active] <- NA; rt[!active] <- Inf
  data.frame(R = R, rt = rt)
}

.btawl_constructor <- function(mode = c("full", "transient", "sustained", "separate"),
                               posdrift = TRUE,
                               drift_distribution = c("lognormal", "normal",
                                                       "splitlognormal", "weibull")) {
  mode <- match.arg(mode)
  drift_distribution <- match.arg(drift_distribution)
  launch <- .ba_launch_code(drift_distribution, "BTAwL")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift))
    stop("BTAwL: posdrift only applies to drift_distribution = \"normal\".")
  if (weibull) {
    p_types <- c(shape = log(1), mean = log(1))
    transform <- c(shape = "exp", mean = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), mean = c(1e-4, Inf))
  } else if (lognormal) {
    p_types <- c(mu = 0, sigma = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
    if (splitlognormal) {
      p_types <- c(p_types, delta = 0)
      transform <- c(transform, delta = "identity")
      minmax <- cbind(minmax, delta = c(-Inf, Inf))
    }
  } else {
    p_types <- c(v = 1, sv = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  launch_pars <- .ba_par_names(launch)
  if (mode == "separate") {
    if (weibull) {
      p_types <- c(shape_S = log(1), mean_S = log(1),
                   shape_T = log(1), mean_T = log(1))
      transform <- c(shape_S = "exp", mean_S = "exp",
                     shape_T = "exp", mean_T = "exp")
      minmax <- cbind(shape_S = c(1e-4, Inf), mean_S = c(1e-4, Inf),
                      shape_T = c(1e-4, Inf), mean_T = c(1e-4, Inf))
      launch_pars <- c("shape_S", "mean_S", "shape_T", "mean_T")
    } else if (lognormal) {
      if (splitlognormal) {
        p_types <- c(mu_S = 0, sigma_S = log(1), delta_S = 0,
                     mu_T = 0, sigma_T = log(1), delta_T = 0)
        transform <- c(mu_S = "identity", sigma_S = "exp", delta_S = "identity",
                       mu_T = "identity", sigma_T = "exp", delta_T = "identity")
        minmax <- cbind(mu_S = c(-Inf, Inf), sigma_S = c(1e-4, Inf),
                        delta_S = c(-Inf, Inf), mu_T = c(-Inf, Inf),
                        sigma_T = c(1e-4, Inf), delta_T = c(-Inf, Inf))
        launch_pars <- c("mu_S", "sigma_S", "delta_S",
                         "mu_T", "sigma_T", "delta_T")
      } else {
        p_types <- c(mu_S = 0, sigma_S = log(1),
                     mu_T = 0, sigma_T = log(1))
        transform <- c(mu_S = "identity", sigma_S = "exp",
                       mu_T = "identity", sigma_T = "exp")
        minmax <- cbind(mu_S = c(-Inf, Inf), sigma_S = c(1e-4, Inf),
                        mu_T = c(-Inf, Inf), sigma_T = c(1e-4, Inf))
        launch_pars <- c("mu_S", "sigma_S", "mu_T", "sigma_T")
      }
    } else {
      p_types <- c(v_S = 1, sv_S = log(1), v_T = 1, sv_T = log(1))
      transform <- c(v_S = "identity", sv_S = "exp",
                     v_T = "identity", sv_T = "exp")
      minmax <- cbind(v_S = c(-Inf, Inf), sv_S = c(1e-4, Inf),
                      v_T = c(-Inf, Inf), sv_T = c(1e-4, Inf))
      launch_pars <- c("v_S", "sv_S", "v_T", "sv_T")
    }
    p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0), k = log(0),
                 tau_s = log(1), tau_t = log(1))
    transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                   tau_s = "exp", tau_t = "exp")
    minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                    t0 = c(0.05, Inf), k = c(0, Inf),
                    tau_s = c(1e-4, Inf), tau_t = c(1e-4, Inf))
    exception <- c(A = 0, k = 0)
    .tw <- add_time_warp_par(p_types, transform, minmax, exception)
    p_types <- .tw$p_types; transform <- .tw$transform
    minmax <- .tw$minmax; exception <- .tw$exception
    .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
    p_types <- .nuis$p_types; transform <- .nuis$transform
    minmax <- .nuis$minmax; exception <- .nuis$exception
    c_name <- paste0("BTAwL_SEPARATE", if (weibull) "_WEIB" else if (lognormal) "_LOGN" else "",
                     if (splitlognormal) "_SPLIT" else "", "_RATE",
                     if (!lognormal && !posdrift) "_IO" else "")
    return(list(
      type = "RACE", c_name = c_name, drift_distribution = drift_distribution,
      p_types = p_types,
      p_types_canonical = setdiff(names(p_types),
                                  c(.time_warp_par_name, .nuisance_par_names)),
      transform = list(func = transform),
      bound = list(minmax = minmax, exception = exception),
      Ttransform = function(pars, dadm) {
        lead <- c(launch_pars, "B", "A", "t0", "k", "tau_s", "tau_t")
        extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
        b <- pars[, "B"] + pars[, "A"]
        cbind(pars[, lead, drop = FALSE], extra, b = b)
      },
      rfun = function(data, pars)
        .rfun_BTAwL(data$lR, pars, ok = attr(pars, "ok"),
                    posdrift = posdrift, launch = launch,
                    mode = "separate", separate = TRUE),
      dfun = function(rt, pars)
        dBTAwLSeparate(rt, pars, launch = launch, posdrift = posdrift),
      pfun = function(rt, pars)
        pBTAwLSeparate(rt, pars, launch = launch, posdrift = posdrift),
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("BTAwLSeparate: the likelihood is implemented in the compiled race path; ",
             "the R likelihood route is not supported.")
      }
    ))
  }
  if (mode == "transient") {
    p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0),
                 k = log(0), tau = log(1))
    transform <- c(transform, B = "exp", A = "exp", t0 = "exp",
                   k = "exp", tau = "exp")
    minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                    t0 = c(0.05, Inf), k = c(1e-4, Inf),
                    tau = c(1e-4, Inf))
    exception <- c(A = 0, k = 0)
    c_name <- paste0("BTAwL_TRANSIENT", if (weibull) "_WEIB" else if (lognormal) "_LOGN" else "",
                     if (splitlognormal) "_SPLIT" else "", "_RATE",
                     if (!lognormal && !posdrift) "_IO" else "")
    .tw <- add_time_warp_par(p_types, transform, minmax, exception)
    p_types <- .tw$p_types; transform <- .tw$transform
    minmax <- .tw$minmax; exception <- .tw$exception
    .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
    p_types <- .nuis$p_types; transform <- .nuis$transform
    minmax <- .nuis$minmax; exception <- .nuis$exception
    return(list(
      type = "RACE", c_name = c_name, drift_distribution = drift_distribution,
      p_types = p_types,
      p_types_canonical = setdiff(names(p_types),
                                  c(.time_warp_par_name, .nuisance_par_names)),
      transform = list(func = transform),
      bound = list(minmax = minmax, exception = exception),
      Ttransform = function(pars, dadm) {
        lead <- c(launch_pars, "B", "A", "t0", "k", "tau")
        extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
        b <- pars[, "B"] + pars[, "A"]
        tau <- pars[, "tau"]
        Tmax_op <- btawl_tmax_vec(pars[, "k"], tau)
        Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
        Vcrit <- btawl_vcrit_vec(pars[, "k"], tau, b)
        out <- cbind(pars[, lead, drop = FALSE], extra, b = b,
                     Tmax = Tmax, rt_max = pars[, "t0"] + Tmax, Vcrit = Vcrit)
        out
      },
      rfun = function(data, pars)
        .rfun_BTAwL(data$lR, pars, ok = attr(pars, "ok"),
                    posdrift = posdrift, launch = launch,
                    mode = "transient"),
      dfun = function(rt, pars)
        dBTAwLTransient(rt, pars, launch = launch, posdrift = posdrift),
      pfun = function(rt, pars)
        pBTAwLTransient(rt, pars, launch = launch, posdrift = posdrift),
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("BTAwLTransient: the likelihood is implemented in the compiled race path; ",
             "the R likelihood route is not supported.")
      }
    ))
  }

  if (mode == "sustained") {
    p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0), k = log(0),
                 tau_s = log(1))
    transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                   tau_s = "exp")
    minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                    t0 = c(0.05, Inf), k = c(0, Inf), tau_s = c(1e-4, Inf))
    exception <- c(A = 0, k = 0)
    c_name <- paste0("BTAwL_SUSTAINED", if (weibull) "_WEIB" else if (lognormal) "_LOGN" else "",
                     if (splitlognormal) "_SPLIT" else "",
                     if (!lognormal && !posdrift) "_IO" else "")
    .tw <- add_time_warp_par(p_types, transform, minmax, exception)
    p_types <- .tw$p_types; transform <- .tw$transform
    minmax <- .tw$minmax; exception <- .tw$exception
    .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
    p_types <- .nuis$p_types; transform <- .nuis$transform
    minmax <- .nuis$minmax; exception <- .nuis$exception
    return(list(
      type = "RACE", c_name = c_name, drift_distribution = drift_distribution,
      p_types = p_types,
      p_types_canonical = setdiff(names(p_types),
                                  c(.time_warp_par_name, .nuisance_par_names)),
      transform = list(func = transform),
      bound = list(minmax = minmax, exception = exception),
      Ttransform = function(pars, dadm) {
        lead <- c(launch_pars, "B", "A", "t0", "k", "tau_s")
        extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
        b <- pars[, "B"] + pars[, "A"]
        cbind(pars[, lead, drop = FALSE], extra, b = b)
      },
      rfun = function(data, pars)
        .rfun_BTAwL(data$lR, pars, ok = attr(pars, "ok"),
                    posdrift = posdrift, launch = launch,
                    mode = "sustained"),
      dfun = function(rt, pars)
        dBTAwLSustained(rt, pars, launch = launch, posdrift = posdrift),
      pfun = function(rt, pars)
        pBTAwLSustained(rt, pars, launch = launch, posdrift = posdrift),
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("BTAwLSustained: the likelihood is implemented in the compiled race path; ",
             "the R likelihood route is not supported.")
      }
    ))
  }

  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0), k = log(0),
               tau_s = log(1), tau_t = log(1), pi = qnorm(.5))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                 tau_s = "exp", tau_t = "exp", pi = "pnorm")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(0, Inf),
                  tau_s = c(1e-4, Inf))
  minmax <- cbind(minmax, tau_t = c(1e-4, Inf), pi = c(0, 1))
  # `pi` is permitted to sit on either endpoint so the two nested submodels are
  # reachable from the full kernel: pi = 0 is a pure transient race (tau_s is
  # then inert) and pi = 1 a pure sustained one (tau_t is inert).  The compiled
  # kernel dispatches both exactly, so each matches its dedicated wrapper.
  exception <- c(A = 0, k = 0, pi = 0, pi = 1)
  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  list(
    type = "RACE",
    c_name = paste0("BTAwL", if (weibull) "_WEIB" else if (lognormal) "_LOGN" else "",
                    if (splitlognormal) "_SPLIT" else "", "_RATE",
                    if (!lognormal && !posdrift) "_IO" else ""),
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types),
                                c(.time_warp_par_name, .nuisance_par_names)),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      lead <- c(launch_pars, "B", "A", "t0", "k", "tau_s", "tau_t", "pi")
      extra <- pars[, setdiff(colnames(pars), lead), drop = FALSE]
      b <- pars[, "B"] + pars[, "A"]
      tau_t <- pars[, "tau_t"]
      pure_transient <- pars[, "pi"] <= 1e-14
      pure_transient[is.na(pure_transient)] <- FALSE
      Tmax_op <- btawl_tmax_vec(pars[, "k"], tau_t)
      Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
      Vcrit <- btawl_vcrit_vec(pars[, "k"], tau_t, b)
      # A sustained component removes the transient hard endpoint.  Keep the
      # transient channel's sampled clearance in the parameter columns, but
      # only expose Tmax/rt_max/Vcrit as endpoint diagnostics for pi = 0.
      Tmax[!pure_transient] <- Inf
      Vcrit[!pure_transient] <- Inf
      out <- cbind(pars[, lead, drop = FALSE], extra, b = b,
                   Tmax = Tmax, rt_max = pars[, "t0"] + Tmax, Vcrit = Vcrit)
      out
    },
    rfun = function(data, pars)
      .rfun_BTAwL(data$lR, pars, ok = attr(pars, "ok"),
                  posdrift = posdrift, launch = launch,
                  mode = "full"),
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

#' Ballistic sustained/transient local race with leak.
#'
#' The full BTAwL model is a local race between independent sustained and
#' transient processes.  `pi` allocates one base launch distribution between
#' those processes; it is not a shared trialwise launch draw.  Normal launches
#' scale both location and spread, lognormal launches shift the log location,
#' and Weibull launches use `V ~ Weibull(shape, scale)` with the arithmetic
#' mean supplied publicly (the conventional scale is derived internally).
#'
#' The model uses the following parameter matrix.  `B` is the distance from the
#' upper end of the start-point range to the threshold, so the kernel receives
#' `b = B + A`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median log launch location for a lognormal launch. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean log width. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength (normal launch only). |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate. |
#' | *tau_s* | log | \[0, Inf\] | log(1) | | Sustained-drive time constant. |
#' | *tau_t* | log | \[0, Inf\] | log(1) | | Transient-drive time constant. |
#' | *pi* | probit | \[0, 1\] | qnorm(.5) | | Probability allocated to the sustained process. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by `v`
#' and `sv`; the launch is truncated positive when `posdrift = TRUE`.
#' With `"splitlognormal"`, `log V` is continuous split normal with
#' `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`;
#' `mu` is the exact median, the split point is derived from that condition,
#' and `delta` is unbounded.  `delta = 0` reduces exactly to lognormal.
#' The split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' With `drift_distribution = "weibull"`, `shape` and arithmetic `mean` are
#' positive; for the local race, `pi` multiplies the mean of the sustained
#' channel and `1 - pi` multiplies the mean of the transient channel.
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
#'
#' The evidence scale is not intrinsic: multiplying the launch strength,
#' threshold distance, and start-point range by the same positive constant
#' leaves response times unchanged.  Fix one scale parameter (normally `sv` for
#' a normal launch, or one intercept among `mu`, `B`, and `A` for a lognormal
#' launch) when specifying a design.  Anchoring on a `B` intercept is usually
#' better conditioned than anchoring on `mu`, which tends to drive `B` onto its
#' lower bound.
#'
#' # Nested submodels
#'
#' Both single-channel models are reachable from the full kernel by pinning
#' `pi`, so a nested model comparison needs no change of `model`:
#'
#' \describe{
#'   \item{`pi = qnorm(0)`}{Pure transient race, identical to
#'     [BTAwLTransient()].}
#'   \item{`pi = qnorm(1)`}{Pure sustained race, identical to
#'     [BTAwLSustained()].}
#' }
#'
#' At a pinned `pi` the other channel's time constant no longer enters the
#' likelihood, so it must be supplied as a constant --- left free it would
#' sample against a flat likelihood and never converge.  Any value inside its
#' bounds gives exactly the same fit; `log(1)` is the conventional choice:
#'
#' ```
#' # pure transient                       # pure sustained
#' constants = c(pi = qnorm(0),           constants = c(pi = qnorm(1),
#'               tau_s = log(1))                        tau_t = log(1))
#' ```
#'
#' `tau_s` and `tau_t` are the direct transient/sustained rate parameters;
#' there is no endpoint parameterisation to switch between.  Do not try to
#' switch a channel off with `tau_s = log(0)` or `tau_t = log(0)`.  Both are
#' outside the bounds, and neither means what it looks like: as `tau_s` falls
#' to zero the sustained drive becomes a *step* at full strength (the fastest
#' sustained input, not the absence of one), while a vanishing `tau_t`
#' delivers zero total transient impulse.  `pi` is the switch.
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"`, `"lognormal"`, or
#'   `"splitlognormal"` (continuous median-parameterised split-lognormal), or
#'   `"weibull"` (`V ~ Weibull(shape, scale)` with public arithmetic `mean`).
#' @return A BTAwL race-model specification.
#' @seealso [BTAwLTransient()], [BTAwLSustained()]
#' @export
BTAwL <- function(posdrift = TRUE,
                  drift_distribution = c("lognormal", "normal", "splitlognormal", "weibull"))
  .btawl_constructor("full", posdrift, drift_distribution)

#' BTAwL local race with independent sustained and transient launch scales.
#'
#' `BTAwLSeparate` defines a local race with independent sustained and
#' transient launch distributions.  It has channel-specific launch location
#' and dispersion parameters (`mu_S`, `mu_T` and `sigma_S`, `sigma_T`, or the
#' corresponding normal/Weibull parameters).  The closed-form survivor is the
#' product of the channel survivors.
#'
#' The model uses the following parameter matrix.  `B` is the distance from the
#' upper end of the start-point range to the threshold, so the kernel receives
#' `b = B + A`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu_S* | identity | \[-Inf, Inf\] | 0 | | Median log launch location for the sustained channel (lognormal launch only). |
#' | *sigma_S* | log | \[0, Inf\] | log(1) | | Geometric-mean log width for the sustained channel (lognormal launch only). |
#' | *mu_T* | identity | \[-Inf, Inf\] | 0 | | Median log launch location for the transient channel (lognormal launch only). |
#' | *sigma_T* | log | \[0, Inf\] | log(1) | | Geometric-mean log width for the transient channel (lognormal launch only). |
#' | *delta_S* | identity | \[-Inf, Inf\] | 0 | | Sustained split-lognormal shape parameter (split-lognormal launch only). |
#' | *delta_T* | identity | \[-Inf, Inf\] | 0 | | Transient split-lognormal shape parameter (split-lognormal launch only). |
#' | *shape_S* | log | \[0, Inf\] | log(1) | | Weibull shape for the sustained channel (Weibull launch only). |
#' | *mean_S* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean for the sustained channel (Weibull launch only). |
#' | *shape_T* | log | \[0, Inf\] | log(1) | | Weibull shape for the transient channel (Weibull launch only). |
#' | *mean_T* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean for the transient channel (Weibull launch only). |
#' | *v_S* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength for the sustained channel (normal launch only). |
#' | *sv_S* | log | \[0, Inf\] | log(1) | | SD of normal launch strength for the sustained channel (normal launch only). |
#' | *v_T* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength for the transient channel (normal launch only). |
#' | *sv_T* | log | \[0, Inf\] | log(1) | | SD of normal launch strength for the transient channel (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate. |
#' | *tau_s* | log | \[0, Inf\] | log(1) | | Sustained-drive time constant. |
#' | *tau_t* | log | \[0, Inf\] | log(1) | | Transient-drive time constant. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' With `drift_distribution = "normal"`, use `v_S`, `sv_S`, `v_T`, and
#' `sv_T`; the launches are truncated positive when `posdrift = TRUE`.
#' With `"lognormal"`, use the `mu_S`, `sigma_S`, `mu_T`, and `sigma_T` rows.
#' With `"splitlognormal"`, each channel has its own `delta` parameter;
#' `delta = 0` reduces that channel exactly to the lognormal launch.  With
#' `drift_distribution = "weibull"`, use the positive `shape` and `mean`
#' rows.  Optional fitting parameters are `pContaminant` (omission
#' probability) and `pGuess` (uniform-outlier probability).
#'
#' **Operational-time warp.** `eta` is a trailing free parameter (identity
#' transform, default `0`, unbounded on the natural scale) and is excluded from
#' `p_types_canonical`.  For physical accumulation time `u = rt - t0`, let
#' `omega = exp(eta)` and `s = ((1 + u)^omega - 1) / omega`.  The CDF and
#' survivor are the parent CDF/survivor evaluated at `s`; the density is the
#' parent density times the Jacobian `c_eta'(u) = (1 + u)^(omega - 1)`.
#' Simulation maps a parent internal finish `s` back with
#' `u = (1 + omega * s)^(1 / omega) - 1` and returns `rt = t0 + u`, so `t0`
#' remains additive.  `eta = 0` is the exact identity warp.
#'
#' For BAwL, this contract applies only to the clock-free, uncorrelated
#' constructor; all non-ballistic models reject `eta`.
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Launch family: `"normal"`, `"lognormal"`,
#'   `"splitlognormal"`, or `"weibull"`.
#' @return A BTAwL race-model specification.
#' @seealso [BTAwL()]
#' @export
BTAwLSeparate <- function(posdrift = TRUE,
                           drift_distribution = c("lognormal", "normal", "splitlognormal", "weibull"))
  .btawl_constructor("separate", posdrift, drift_distribution)

#' Pure transient BTAwL wrapper.
#'
#' The transient wrapper uses the following parameter matrix:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median log launch location. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean log width. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength (normal launch only). |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate. |
#' | *tau* | log | \[0, Inf\] | log(1) | | Transient time constant. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by `v`
#' and `sv`; the launch is truncated positive when `posdrift = TRUE`.
#' With `"splitlognormal"`, `log V` is continuous split normal with
#' `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`;
#' `mu` is the exact median and `delta` is unbounded.  The split point is
#' derived from the median condition, and `delta = 0` is exactly lognormal.
#' The split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' With `drift_distribution = "weibull"`, use positive `shape` and `mean`
#' instead of the lognormal launch rows.
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
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"`, `"lognormal"`, or
#'   `"splitlognormal"` (continuous median-parameterised split-lognormal), or
#'   `"weibull"` (`V ~ Weibull(shape, scale)` with public arithmetic `mean`).
#' @return A transient-only BTAwL race-model specification.
#' @export
BTAwLTransient <- function(posdrift = TRUE,
                           drift_distribution = c("lognormal", "normal",
                                                   "splitlognormal", "weibull"))
  .btawl_constructor("transient", posdrift, drift_distribution)
#' Pure sustained BTAwL wrapper.
#'
#' The sustained wrapper uses the following parameter matrix:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median log launch location. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean log width. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength (normal launch only). |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate. |
#' | *tau_s* | log | \[0, Inf\] | log(1) | | Sustained-drive time constant. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by `v`
#' and `sv`; the launch is truncated positive when `posdrift = TRUE`.
#' With `"splitlognormal"`, `log V` is continuous split normal with
#' `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`;
#' `mu` is the exact median and `delta` is unbounded.  The split point is
#' derived from the median condition, and `delta = 0` is exactly lognormal.
#' The split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' With `drift_distribution = "weibull"`, use positive `shape` and `mean`
#' instead of the lognormal launch rows.
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
#'
#' @param posdrift Logical. For a normal launch, truncate `V` below zero.
#' @param drift_distribution Either `"normal"`, `"lognormal"`, or
#'   `"splitlognormal"` (continuous median-parameterised split-lognormal), or
#'   `"weibull"` (`V ~ Weibull(shape, scale)` with public arithmetic `mean`).
#' @return A sustained-only BTAwL race-model specification.
#' @export
BTAwLSustained <- function(posdrift = TRUE,
                           drift_distribution = c("lognormal", "normal",
                                                   "splitlognormal", "weibull"))
  .btawl_constructor("sustained", posdrift, drift_distribution)
