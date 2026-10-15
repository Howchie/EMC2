# What each arm's chains bought at each rebuild, and the plan's own predicted tau against the realized one.
suppressPackageStartupMessages(library(posterior))
A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15"
files <- c(old_1 = file.path(A, "arms/old_1.rds"), new_1 = file.path(A, "arms/new_1.rds"),
  new_2 = file.path(A, "handoff_2026-10-15/arms/pair_new_2.rds"), new_3 = file.path(A, "handoff_2026-10-15/arms/pair_new_3.rds"),
  point_1 = file.path(A, "handoff_2026-10-15/arms/arm_point_1.rds"), point_2 = file.path(A, "handoff_2026-10-15/arms/pair_point_2.rds"),
  point_3 = file.path(A, "handoff_2026-10-15/arms/pair_point_3.rds"), kappa_1 = file.path(A, "handoff_2026-10-15/arms/arm_kappa_1.rds"),
  transport_1 = file.path(A, "eval_transport/arm_transport_1.rds"))
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-9); (1 + r) / (1 - r) }
for (nm in names(files)) {
  r <- readRDS(files[nm]); G <- r$draws
  run <- apply(simplify2array(lapply(G, function(g) g[101:400, ])), 2, function(x) nrow(x) * ncol(x) / ess_basic(x))
  slow <- order(-run)[1:6]
  cat(sprintf("\n==== %s: cpu %.0f  ridge ESS/kCPU %.2f  median ESS/kCPU %.2f  moves(hmc/shear/coll per chain): %s\n", nm, r$cpu,
    1000 * min(r$ess[slow]) / r$cpu, 1000 * median(r$ess) / r$cpu, paste(apply(r$moves, 2, paste, collapse = "/"), collapse = "  ")))
  for (j in 1:3) for (h in r$plans[[j]]) {
    it <- h$iteration; if (is.null(h$options)) next
    cat(sprintf(" ch%d it%d %-5s obj %6.1f (none %6.1f) tau %5.1f tau0 %5.1f cost %5.2f C_rest %4.2f | rounds %4.2f coll %4.2f | hmc mean %.2f (>0: %3d) q mean %.2f (<1: %3d) | pred tau slow6: %s | subj rho max %.2f\n",
      j, it, h$choice, h$objective, h$options[["none"]], h$tau, h$tau0, h$cost, h$C_rest, sum(h$rounds), h$collapsed_rounds,
      mean(h$rate), sum(h$rate > 0), mean(h$particle_shared), sum(h$particle_shared < 1),
      paste(round(tau1(h$pred_group[slow])), collapse = " "), max(apply(h$pred_subject, 1, max))))
  }
  cat(sprintf("   realized tau slow6: %s   (%s)\n", paste(round(run[slow]), collapse = " "), paste(names(r$ess)[slow], collapse = " ")))
}
