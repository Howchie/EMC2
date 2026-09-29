## How well does the frozen ensemble proposal q_s cover each subject's recent
## posterior draws? Squared Mahalanobis distance of the chain's own last-200
## draws under q_s's normal component: E = p (17) if q matches the target.
suppressMessages(library(EMC2))
chk <- function(f, stage_pick) {
  load(f); st <- emc[[1]]$samples$stage
  idx <- tail(which(st == stage_pick), 200)
  out <- sapply(seq_along(emc), function(k) {
    q <- emc[[k]]$ensemble_q
    sapply(seq_along(q$mu), function(s) {
      a <- emc[[k]]$samples$alpha[, s, idx]
      z <- forwardsolve(t(chol(q$var[[s]])), a - q$mu[[s]])
      c(m2 = mean(colSums(z^2)),
        # posterior spread relative to q along q's axes: >1 means q too narrow
        ratio = max(eigen(cov2cor(cov(t(z))) * 0 + cov(t(z)), only.values = TRUE)$values),
        shift = sqrt(sum(rowMeans(z)^2)))
    })
  }, simplify = "array")      # 3 x subj x chain
  cat(f, "(", stage_pick, "last 200 )\n")
  for (nm in dimnames(out)[[1]]) cat(sprintf("  %-5s median %.2f  90%% %.2f  max %.2f\n", nm,
      median(out[nm, , ]), quantile(out[nm, , ], .9), max(out[nm, , ])))
}
chk("/data/work/EMC2_dev_oo/plans/convergence/samples_lba_rates_d45_snapshot_0218.RData", "sample")
chk("/data/work/PM/NirvanaHons_Nback/Fits/Match_Control_Exp1_rdm7s_mass_massR_noA_rates_t0R_UT0_vs2_mt_d45_ens.RData", "sample")
