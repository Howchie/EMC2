# The allocator's inputs for the ridge coordinate at a saved state: window variance, rest-of-iteration jump, shear record E/n, path displacement.
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
emc <- readRDS(commandArgs(TRUE)[1])
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
u <- readRDS("/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/handoff_2026-10-15/scripts/ridge_u.rds")
for (j in 1:3) {
  ch <- emc[[j]]; w <- ch$alloc$window; idx <- ch$samples$idx
  dirs <- ch$shear$plan$dirs; Dinv <- solve(dirs)
  # which shear direction is the ridge
  cosang <- abs(drop(t(dirs) %*% u)) / sqrt(colSums(dirs^2) * sum(u^2))
  k <- which.max(cosang)
  draws <- (idx - w$n + 1):idx
  X <- do.call(rbind, lapply(emc, function(x) t(Dinv %*% x$samples$theta_mu[, draws])))
  V <- apply(X, 2, var); Wv <- mean(sapply(emc, function(x) apply(t(Dinv %*% x$samples$theta_mu[, draws]), 2, var)[k]))
  T_g <- rbind(cbind(Dinv, matrix(0, p, p)), cbind(matrix(0, p, p), diag(p))); rest_esjd <- (T_g %*% w$G %*% t(T_g))[k, k] / w$nG
  rec <- ns$.emc_shear_record(ch$shear$plan$last, p)
  last <- ch$alloc$last
  cat(sprintf("ch%d: ridge = shear dir %d (|cos| %.2f); window n %d; V pooled %.3g within %.3g (B/V %.2f)\n", j, k, cosang[k], w$n, V[k], Wv, 1 - Wv / V[k]))
  cat(sprintf("   rest-of-iteration ESJD %.3g -> rho %.3f (recorded group_rho %.3f); shear E/n per round %.3g over %.0f rounds (moved %.0f) -> g_obs %.3f, recorded g %.3f a %.3f; rounds planned %.2f; N/K %s\n",
    rest_esjd, 1 - rest_esjd / (2 * V[k]), last$group_rho[k], rec$E[k] / rec$n[k], rec$n[k], rec$moved[k], 1 - rec$E[k] / rec$n[k] / (2 * V[k]), last$group_g[k], last$group_a[k], last$rounds[k],
    paste(round(rec$N[k] / max(1, sum(rec$K[k, ])), 4), collapse = " ")))
  # whole-iteration lag-1 and lag-10 jumps of stored draws in this coordinate
  x <- drop(Dinv[k, ] %*% ch$samples$theta_mu[, draws])
  cat(sprintf("   stored draws: lag-1 jump %.3g, lag-10 jump/10 %.3g, lag-25 jump/25 %.3g (within-chain)\n", mean(diff(x)^2), mean(diff(x, 10)^2) / 10, mean(diff(x, 25)^2) / 25))
  cat(sprintf("   predicted tau for ridge under plan %.1f; tau0 %.1f; plan choice %s\n", (1 + last$pred_group[k]) / (1 - last$pred_group[k]), last$tau0, last$choice))
}
