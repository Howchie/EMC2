lib <- Sys.getenv("LIB"); suppressMessages(library(EMC2, lib.loc = lib)); ns <- asNamespace("EMC2")
load("Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
stage <- "adapt"; B <- 100; p <- emc[[1]]$n_pars
Vbar <- lapply(emc, function(ch) { k <- ns$.emc_alloc_draws(ch, max(2L, ch$alloc$window$n)); apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean) })
for (j in 1:3) {
  hst <- tryCatch(ns$.emc_hmc_plan_chain(emc[[j]], stage, Vbar[[j]]), error = function(e) NULL)
  ms <- ns$.emc_alloc_measure(emc, j, hst, Vbar[[j]], stage)
  u <- ns$.emc_alloc_units(ms, hst, ms$shear, optimistic = FALSE, B = B)
  sol <- ns$.emc_alloc_solve_units(u)
  w <- vapply(u$units[1:(2*p)], function(x) if (is.null(x$w)) NA else x$w, 0)
  cat(sprintf("chain %d: s0 %s | choice %s | options %s | tau %.1f\n", j, "", sol$choice,
      paste(names(sol$options), round(sol$options, 1), collapse = " "), sol$tau))
  if (!is.null(ms$g_lag)) cat("  lag1 lsd:", round(ms$g_lag$r1[p + 1:p], 2), "\n  lag2 lsd:", round(ms$g_lag$r2[p + 1:p], 2),
     "\n  s lsd:", round(ms$g_lag$s[p + 1:p], 2), "\n  s mu:", round(ms$g_lag$s[1:p], 2), "\n  w mu:", round(w[1:p], 2), "\n  w lsd:", round(w[p + 1:p], 2), "\n")
}
