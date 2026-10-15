suppressPackageStartupMessages(library(posterior))
for (f in commandArgs(TRUE)) {
  r <- readRDS(f); G <- r$draws; n <- nrow(G[[1]]); d <- ncol(G[[1]])
  arr <- simplify2array(lapply(G, function(g) g[-(1:100), , drop = FALSE]))
  tau_run <- apply(arr, 2, function(a) nrow(a) * ncol(a) / ess_basic(a))
  pr_bm <- pr_l10 <- matrix(NA, 4, d)
  for (b in 1:4) {
    idx <- (b - 1) * 100 + 1:100; blk <- lapply(G, function(g) g[idx, , drop = FALSE])
    V <- apply(do.call(rbind, blk), 2, var)
    means <- sapply(blk, colMeans)                       # coord x chain
    pr_bm[b, ] <- pmax(1, 100 * apply(means, 1, var) / V)
    L <- 10; T <- Reduce(`+`, lapply(blk, function(g) colMeans((g[-(1:L), ] - g[1:(100 - L), ])^2))) / 3 / (L * V)
    rho <- pmax(0, 1 - T * L / 2)^(1 / L); pr_l10[b, ] <- (1 + rho) / (1 - pmin(rho, 1 - 1e-9))
  }
  slow <- order(-tau_run)[1:8]
  res <- sapply(list(L10 = pr_l10, bm = pr_bm, max = pmax(pr_l10, pr_bm)), function(pr) { tp <- apply(pr, 2, median)
    c(ratio_all = median(tp / tau_run), ratio_slow8 = median(tp[slow] / tau_run[slow]), spearman = cor(tp, tau_run, method = "spearman"),
      blockcv_slow8 = median(apply(pr[, slow], 2, function(x) sd(log(x))))) })
  cat(basename(f), "\n"); print(round(res, 3))
}
