#### C++ rfun dispatch (rfun_port_plan.md) ----
#
# Thin R-side wrappers around the C++ simulation kernels (src/model_rng.cpp)
# for LBA/LBAIO/logical-rules race models, BAwL, RDM and RDMSWTN. Each wrapper falls
# back to the pure-R reference rfun (.lba_rfun/rBAwL/rRDM/rRDMSWTN) when the
# `emc2.cpp_rfun` option is turned off. The kernels are distributionally,
# not stream-, equivalent to the R rfuns -- set.seed()-reproduced simulated
# datasets will differ trial-by-trial between the two paths.

# The R reference rfuns size their output data.frame by length(lR)/n_acc but
# do their internal race computation over nrow(pars)/n_acc trials; when a
# caller passes fewer pars rows than lR (only ever exercised by hand-built
# test inputs, not by make_data()/RACE_rfun), R's `out$R <- ...` data.frame
# column assignment silently recycles the shorter result to fill the longer
# column. Replicate that here so the C++ path matches for such inputs.
.rfun_cpp_pack <- function(res, lR_levels, n_trials_target) {
  n <- length(res$R)
  if (n_trials_target != n) {
    if (n == 0L || n_trials_target %% n != 0L)
      stop("replacement has ", n, " rows, data has ", n_trials_target)
    idx <- rep_len(seq_len(n), n_trials_target)
    res$R <- res$R[idx]
    res$rt <- res$rt[idx]
    if (!is.null(res$isTime)) res$isTime <- res$isTime[idx]
  }
  out <- data.frame(R = factor(res$R, levels = seq_along(lR_levels), labels = lR_levels),
                    rt = res$rt)
  if (!is.null(res$isTime)) out$isTime <- res$isTime
  out
}

.use_cpp_rfun <- function() isTRUE(getOption("emc2.cpp_rfun", TRUE))

.rfun_LBA <- function(lR, pars, ok = rep(TRUE, length(lR)), posdrift = TRUE) {
  if (.use_cpp_rfun()) {
    res <- rlba_cpp(pars, levels(lR), ok, posdrift)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  .lba_rfun(lR, pars, ok = ok, posdrift = posdrift)
}

.rfun_RDM <- function(lR, pars, ok = rep(TRUE, dim(pars)[1])) {
  if (.use_cpp_rfun()) {
    res <- rrdm_cpp(pars, levels(lR), ok)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rRDM(lR, pars, ok = ok)
}

.rfun_ROU <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                      par = "rate") {
  if (.use_cpp_rfun()) {
    res <- rrou_cpp(
      pars, levels(lR), ok, kind,
      dt = getOption("emc2.rou_sim_dt", 1e-3),
      t_max = getOption("emc2.rou_sim_tmax", 30),
      par_kind = .ROU_PAR[[par]]
    )
    out <- .rfun_cpp_pack(
      res, levels(lR), length(lR) / length(levels(lR))
    )
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  .rfun_ROU_R(lR, pars, ok = ok, kind = kind, par = par)
}

.rfun_ROUp <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                       par = "rate", pooling = "coactive") {
  pooling <- match.arg(pooling, names(.ROUp_POOLING))
  if (.use_cpp_rfun()) {
    res <- rroup_cpp(
      pars, levels(lR), ok, kind,
      dt = getOption("emc2.roup_sim_dt", 1e-3),
      t_max = getOption("emc2.roup_sim_tmax", 30),
      par_kind = .ROUp_PAR[[par]],
      pooling = .ROUp_POOLING[[pooling]]
    )
    out <- .rfun_cpp_pack(
      res, levels(lR), length(lR) / length(levels(lR))
    )
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  .rfun_ROUp_R(lR, pars, ok = ok, kind = kind, par = par, pooling = pooling)
}



.rfun_BAwL <- function(lR, pars, ok = rep(TRUE, length(lR)), posdrift = TRUE,
                       erlang = 1L, guess = FALSE, global = FALSE, launch = 0L) {
  if (.use_cpp_rfun()) {
    res <- rbawl_cpp(pars, levels(lR), ok, posdrift, as.integer(erlang), guess, global,
                     as.integer(launch))
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rBAwL(lR, pars, ok = ok, posdrift = posdrift, erlang = erlang, guess = guess,
        global = global, launch = launch)
}

# `launch` must come from the same .bawd_launch_code() call that produced the
# model's c_name suffix, or the simulator and the likelihood describe different
# models (see R/model_BAwD.R).
.rfun_BAwD <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                       posdrift = TRUE, gamma = 0, rho = Inf) {
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)
  if (.use_cpp_rfun()) {
    res <- rbawd_cpp(pars, levels(lR), ok, as.integer(launch), posdrift, gamma,
                     rho)
    out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  rBAwD(lR, pars, ok = ok, launch = launch, posdrift = posdrift, gamma = gamma,
        rho = rho)
}

.rfun_BAwF <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                       posdrift = TRUE, rho = Inf) {
  rho <- .bawf_check_rho(rho)
  if (.use_cpp_rfun()) {
    res <- rbawf_cpp(pars, levels(lR), ok, as.integer(launch), posdrift, rho)
    out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  rBAwF(lR, pars, ok = ok, launch = launch, posdrift = posdrift, rho = rho)
}

.rfun_BAwR <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                       posdrift = TRUE) {
  if (.use_cpp_rfun()) {
    res <- rbawr_cpp(pars, levels(lR), ok, as.integer(launch), posdrift)
    out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  rBAwR(lR, pars, ok = ok, launch = launch, posdrift = posdrift)
}

.rfun_BAwDp <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                        posdrift = TRUE) {
  if (.use_cpp_rfun()) {
    res <- rbawdp_cpp(pars, levels(lR), ok, as.integer(launch), posdrift)
    out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  rBAwDp(lR, pars, ok = ok, launch = launch, posdrift = posdrift)
}

.rfun_FRQ <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  if (.use_cpp_rfun()) {
    res <- rfrq_cpp(pars, levels(lR), ok)
    out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
    return(.apply_timed_guess_winner(out, levels(lR)))
  }
  rFRQ(lR, pars, ok = ok)
}

.rfun_BAwL_corr <- function(lR, pars, ok = rep(TRUE, length(lR)), posdrift = TRUE,
                            erlang = 1L, guess = FALSE, global = FALSE) {
  if (.use_cpp_rfun()) {
    res <- rbawl_corr_cpp(pars, levels(lR), ok, posdrift, as.integer(erlang), guess, global)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rBAwL_corr(lR, pars, ok = ok, posdrift = posdrift, erlang = erlang,
             guess = guess, global = global)
}

.rfun_RDMSWTN <- function(lR, pars, ok = rep(TRUE, dim(pars)[1]), erlang_shape = 1L,
                          erlang_type = "none", posdrift = TRUE,
                          correlated = FALSE, correlate = "times") {
  drift_corr <- correlated && correlate == "drifts"
  if (.use_cpp_rfun()) {
    res <- if (drift_corr) {
      rrdmswtn_drift_corr_cpp(pars, levels(lR), ok, posdrift)
    } else if (correlated) {
      rrdmswtn_corr_cpp(pars, levels(lR), ok, posdrift)
    } else {
      rrdmswtn_cpp(pars, levels(lR), ok, as.integer(erlang_shape),
                   erlang_type, posdrift)
    }
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  if (drift_corr) {
    rRDMSWTN(lR, pars, ok = ok, erlang_shape = erlang_shape,
             erlang_type = erlang_type, posdrift = posdrift,
             drift_override = .draw_correlated_drifts(pars, lR, ok, posdrift))
  } else if (correlated) {
    rRDMSWTN_corr(lR, pars, ok = ok, posdrift = posdrift)
  } else {
    rRDMSWTN(lR, pars, ok = ok, erlang_shape = erlang_shape,
             erlang_type = erlang_type, posdrift = posdrift)
  }
}

.rfun_RDMSWTN_TT <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                             posdrift = TRUE, correlated = FALSE,
                             correlate = "times") {
  drift_corr <- correlated && correlate == "drifts"
  if (.use_cpp_rfun()) {
    res <- if (drift_corr) {
      rrdmswtn_tt_drift_corr_cpp(pars, levels(lR), ok, posdrift)
    } else if (correlated) {
      rrdmswtn_tt_corr_cpp(pars, levels(lR), ok, posdrift)
    } else {
      rrdmswtn_tt_cpp(pars, levels(lR), ok, posdrift)
    }
    return(.rfun_cpp_pack(
      res, levels(lR), length(lR) / length(levels(lR))
    ))
  }
  if (drift_corr) {
    rRDMSWTN_TT(lR, pars, ok = ok, posdrift = posdrift,
                drift_override = .draw_correlated_drifts(pars, lR, ok, posdrift))
  } else if (correlated) {
    rRDMSWTN_TT_corr(lR, pars, ok = ok, posdrift = posdrift)
  } else {
    rRDMSWTN_TT(lR, pars, ok = ok, posdrift = posdrift)
  }
}

# R counterpart of drift_factor_draw_correlated() in src/drift_factor.h: one
# shared standard-normal factor per trial, an independent residual per row, and
# joint rejection when posdrift restricts the draws to the positive orthant.
.draw_correlated_drifts <- function(pars, lR, ok, posdrift, max_iter = 100000) {
  n_acc <- length(levels(lR))
  n_rows <- nrow(pars)
  if (n_acc <= 0L || n_rows %% n_acc != 0L) {
    stop("Correlated drift draws need whole trials of accumulator rows.")
  }
  rho <- pars[, "rho"]
  if (any(!is.finite(rho[ok])) || any(abs(rho[ok]) > 1)) {
    stop("Correlated drift draws require finite rho in [-1, 1].")
  }
  sv <- pars[, "sv"]
  magnitude <- abs(rho)
  slope <- ifelse(rho < 0, -1, 1) * sv * sqrt(magnitude)
  sv_res <- sv * pmax(sqrt(pmax(0, 1 - magnitude)), 1e-12)
  out <- rep(Inf, n_rows)
  for (tr in seq_len(n_rows / n_acc)) {
    rows <- ((tr - 1L) * n_acc + 1L):(tr * n_acc)
    active <- rows[ok[rows]]
    if (!length(active)) next
    accepted <- FALSE
    for (iter in seq_len(max_iter)) {
      z <- rnorm(1)
      draw <- rnorm(length(active), pars[active, "v"] + slope[active] * z,
                    sv_res[active])
      out[active] <- draw
      if (!posdrift || all(draw > 0)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      stop("Correlated drift draws: jointly positive rejection exceeded ",
           max_iter, " attempts; check that the drift means are not far ",
           "below zero.")
    }
  }
  out
}
