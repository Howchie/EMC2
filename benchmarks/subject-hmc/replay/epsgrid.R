source("benchmarks/subject-hmc/replay/setup11.R")
library(parallel)
lap <- ns$.emc_laplace_conditional(structure(numeric(p), lik = ch$lik_laplace[[s]]), list(mu = mu, var = V))
W <- .emc_hmc_ref_whiten(H, V)
# Bracket from the stuck state and from the bulk.
set.seed(1); t0 <- proc.time()[[3]]
b1 <- .emc_hmc_bracket(a_stuck, ll_stuck, fr, ll_fn, 1e-4, floor_s, p^(-1/4))
b2 <- .emc_hmc_bracket(lap$mu, ll_fn(matrix(lap$mu, 1)), fr, ll_fn, 1e-4, floor_s, p^(-1/4))
cat(sprintf("bracket eps: from stuck %.3f (%d trials, %d rows), from Laplace mean %.3f (%d trials) | %.1f s, %.2f ms/row\n",
  b1$eps, b1$trials, b1$rows, b2$eps, b2$trials, proc.time()[[3]] - t0, 1000 * (proc.time()[[3]] - t0) / rows_used))
grid <- c(0.08, 0.12, 0.17, 0.25, 0.35, 0.5, 0.7)
res <- mclapply(grid, function(eps) {
  set.seed(100 + round(eps * 100)); rows_used <<- 0
  a <- lap$mu; ll <- ll_fn(matrix(a, 1)); n <- 120
  out <- matrix(NA, n, 6); Ahist <- numeric(n)
  for (t in 1:n) { tr <- .emc_hmc_transition(a, ll, fr, ll_fn, eps, 1e-4, floor_s)
    wd <- drop(W %*% tr$disp); out[t, ] <- c(tr$acc_prob, mean(wd^2) * tr$acc_prob, min(colMeans(matrix(wd^2, 1))), tr$rows, tr$L, as.numeric(tr$status == "ok"))
    a <- tr$alpha; ll <- tr$ll; Ahist[t] <- a[iA] }
  c(eps = eps, acc = mean(out[, 1]), esjd = mean(out[, 2]), rows = mean(out[, 4]), L = mean(out[, 5]), ok = mean(out[, 6]),
    esjd_per_krow = 1000 * mean(out[, 2]) / mean(out[, 4]), essA = ess_bulk(Ahist[-(1:20)]), essA_per_krow = 1000 * ess_bulk(Ahist[-(1:20)]) / sum(out[-(1:20), 4]),
    sdA = sd(Ahist[-(1:20)]))
}, mc.cores = length(grid))
print(round(do.call(rbind, res), 4))
