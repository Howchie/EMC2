lib <- Sys.getenv("LIB"); suppressMessages(library(EMC2, lib.loc = lib)); ns <- asNamespace("EMC2")
load("Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names; tau <- ns$.emc_alloc_tau
Vbar <- lapply(emc, function(ch) { k <- ns$.emc_alloc_draws(ch, max(2L, ch$alloc$window$n)); apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean) })
out <- NULL
for (j in 1:3) {
  hst <- ns$.emc_hmc_plan_chain(emc[[j]], "adapt", Vbar[[j]])
  ms <- ns$.emc_alloc_measure(emc, j, hst, Vbar[[j]], "adapt")
  for (use in c(FALSE, TRUE)) { m <- ms; if (!use) m$g_kappa <- rep(1, 2 * p)
    sol <- ns$.emc_alloc_solve_units(ns$.emc_alloc_units(m, hst, m$shear, optimistic = FALSE, B = 100))
    out <- rbind(out, data.frame(ch = j, kappa_on = use, coord = 1:(2 * p), rho = ms$g_rho, kappa = ms$g_kappa, g = ms$g_fac, x = sol$x[1:(2 * p)],
      pred = unlist(sol$pred[1:(2 * p)]))) } }
saveRDS(out, file.path(Sys.getenv("SP"), "replay_k2.rds"))
a <- aggregate(cbind(rho, kappa, g, x) ~ coord + kappa_on, out, mean)
w <- merge(a[!a$kappa_on, c("coord", "rho", "kappa", "g", "x")], a[a$kappa_on, c("coord", "x")], by = "coord", suffixes = c("_off", "_on"))
w$lab <- c(paste0("mu", 1:p), paste0("lsd.", pn))[w$coord]
print(head(w[order(-w$x_off - w$x_on), c("lab", "rho", "kappa", "g", "x_off", "x_on")], 12), digits = 2, row.names = FALSE)
cat("total rounds (mean over chains) off", round(sum(w$x_off), 2), " on", round(sum(w$x_on), 2), "\n")
