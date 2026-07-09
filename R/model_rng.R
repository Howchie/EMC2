#### C++ rfun dispatch (rfun_port_plan.md) ----
#
# Thin R-side wrappers around the C++ simulation kernels (src/model_rng.cpp)
# for LBA/LBAIO/LogicalRulesLBA, BAwL, RDM and RDMSWTN. Each wrapper falls
# back to the pure-R reference rfun (rLBA/rBAwL/rRDM/rRDMSWTN) when the
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
  rLBA(lR, pars, ok = ok, posdrift = posdrift)
}

.rfun_RDM <- function(lR, pars, ok = rep(TRUE, dim(pars)[1])) {
  if (.use_cpp_rfun()) {
    res <- rrdm_cpp(pars, levels(lR), ok)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rRDM(lR, pars, ok = ok)
}

.rfun_BAwL <- function(lR, pars, ok = rep(TRUE, length(lR)), posdrift = TRUE,
                       erlang = 1L, guess = FALSE, global = FALSE) {
  if (.use_cpp_rfun()) {
    res <- rbawl_cpp(pars, levels(lR), ok, posdrift, as.integer(erlang), guess, global)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rBAwL(lR, pars, ok = ok, posdrift = posdrift, erlang = erlang, guess = guess, global = global)
}

.rfun_RDMSWTN <- function(lR, pars, ok = rep(TRUE, dim(pars)[1]), erlang_shape = 1L,
                          erlang_type = "none", posdrift = TRUE) {
  if (.use_cpp_rfun()) {
    res <- rrdmswtn_cpp(pars, levels(lR), ok, as.integer(erlang_shape), erlang_type, posdrift)
    return(.rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR))))
  }
  rRDMSWTN(lR, pars, ok = ok, erlang_shape = erlang_shape, erlang_type = erlang_type, posdrift = posdrift)
}
