source("benchmarks/subject-hmc/replay/setup11.R")
library(parallel)
lap <- ns$.emc_laplace_conditional(structure(numeric(p), lik = ch$lik_laplace[[s]]), list(mu = mu, var = V))
W <- .emc_hmc_ref_whiten(H, V)
sj <- list(eps = NA_real_, n_tuned = 0L)
# (1) Stochastic eps rule vs grid: windows of NW attempts from a moving chain.
rule_run <- function(seed, NW, n_win = 8) {
  set.seed(seed); a <- a_stuck; ll <- ll_stuck; sj$eps <- 0.172; path <- sj$eps
  for (w in 1:n_win) {
    ev <- .emc_hmc_new_evidence(p, "k")
    for (k in 1:NW) { tr <- .emc_hmc_transition(a, ll, fr, ll_fn, sj$eps, 1e-4, floor_s)
      tr$score <- tr$acc_prob * mean(drop(W %*% tr$disp)^2); tr$jump <- tr$alpha - a
      ev <- .emc_hmc_add_evidence(ev, tr); a <- tr$alpha; ll <- tr$ll }
    sj <- .emc_hmc_update_eps(sj, ev); path <- c(path, sj$eps)
  }
  path
}
# (2) Pilot batch size: allocation-relevant estimate (min-coordinate rho_H).
pilot_run <- function(seed, n) {
  set.seed(seed); a <- a_stuck; ll <- ll_stuck; ev <- .emc_hmc_new_evidence(p, "k")
  for (k in 1:n) { tr <- .emc_hmc_transition(a, ll, fr, ll_fn, 0.172, 1e-4, floor_s)
    tr$score <- tr$acc_prob * mean(drop(W %*% tr$disp)^2); tr$jump <- tr$alpha - a
    ev <- .emc_hmc_add_evidence(ev, tr); a <- tr$alpha; ll <- tr$ll }
  rho <- .emc_hmc_rho(.emc_hmc_esjd(ev$M, ev$n, W))
  c(n = n, rhoH_max = max(rho), rhoH_med = median(rho), acc = mean(ev$accp))
}
jobs <- c(lapply(1:4, function(k) list(type = "rule", seed = k, NW = 40)),
          lapply(1:6, function(k) list(type = "pilot", seed = 100 + k, n = 10)),
          lapply(1:6, function(k) list(type = "pilot", seed = 200 + k, n = 20)),
          lapply(1:6, function(k) list(type = "pilot", seed = 300 + k, n = 40)))
res <- mclapply(jobs, function(j) if (j$type == "rule") rule_run(j$seed, j$NW) else pilot_run(j$seed, j$n), mc.cores = 2)
cat("eps paths under the jitter-halves rule (40 attempts per window):\n")
for (k in which(sapply(jobs, `[[`, "type") == "rule")) cat("  ", paste(sprintf("%.3f", res[[k]]), collapse = " "), "\n")
pil <- do.call(rbind, res[sapply(jobs, `[[`, "type") == "pilot"])
print(aggregate(cbind(rhoH_max, rhoH_med, acc) ~ n, data = as.data.frame(pil), FUN = function(x) sprintf("%.3f (sd %.3f)", mean(x), sd(x))))
