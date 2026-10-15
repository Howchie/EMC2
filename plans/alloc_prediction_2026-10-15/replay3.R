lib <- Sys.getenv("LIB"); suppressMessages(library(EMC2, lib.loc = lib)); ns <- asNamespace("EMC2")
load("Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
Vbar <- lapply(emc, function(ch) { k <- ns$.emc_alloc_draws(ch, max(2L, ch$alloc$window$n)); apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean) })
j <- 1
hst <- ns$.emc_hmc_plan_chain(emc[[j]], "adapt", Vbar[[j]])
ms <- ns$.emc_alloc_measure(emc, j, hst, Vbar[[j]], "adapt")
u <- ns$.emc_alloc_units(ms, hst, ms$shear, optimistic = FALSE, B = 100)
tau <- ns$.emc_alloc_tau
for (ch in c("0", "1")) {
  uu <- u; uu$regime_fixed <- ch
  sol <- ns$.emc_alloc_solve_units(uu)
  tg <- tau(unlist(sol$pred[1:(2 * p)]))
  cat(sprintf("regime %s: objective %.1f cost %.2f | max tau group %.1f (%s) | lsd taus: %s\n", ch, sol$objective, sol$cost,
    max(tg), c(paste0("mu", 1:p), paste0("lsd.", pn))[which.max(tg)], paste(round(tg[p + 1:p], 1), collapse = " ")))
}
