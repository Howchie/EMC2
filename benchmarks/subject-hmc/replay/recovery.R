# Early-failure recovery: a step size far past the acceptance cliff, then the
# rebuild rule applied to windows of 20 attempts from the stuck state.
source("benchmarks/subject-hmc/replay/setup11.R")
W <- .emc_hmc_ref_whiten(H, V)
for (seed in 1:2) {
  set.seed(seed); a <- a_stuck; ll <- ll_stuck; sj <- list(eps = 0.7, n_tuned = 0L); path <- sj$eps; acc <- NULL
  for (w in 1:8) {
    ev <- .emc_hmc_new_evidence(p, "k")
    for (k in 1:20) { tr <- .emc_hmc_transition(a, ll, fr, ll_fn, sj$eps, 1e-4, floor_s)
      tr$score <- tr$acc_prob * mean(drop(W %*% tr$disp)^2); tr$jump <- tr$alpha - a
      ev <- .emc_hmc_add_evidence(ev, tr); a <- tr$alpha; ll <- tr$ll }
    acc <- c(acc, mean(ev$accp)); sj <- .emc_hmc_update_eps(sj, ev); path <- c(path, sj$eps)
  }
  cat(sprintf("seed %d eps path: %s\n        window acceptance: %s | final A %.3f\n", seed,
      paste(sprintf("%.3f", path), collapse = " "), paste(sprintf("%.2f", acc), collapse = " "), a[iA]))
}
