suppressPackageStartupMessages(library(posterior))
files <- commandArgs(TRUE)
geyer_tau <- function(a) { n <- nrow(a); ess <- ess_basic(a); n * ncol(a) / ess }
for (f in files) {
  r <- readRDS(f)
  G <- r$draws; n <- nrow(G[[1]]); d <- ncol(G[[1]])
  sc <- lapply(G, function(g) g[-(1:100), , drop = FALSE])
  arr <- simplify2array(sc)                       # iter x coord x chain
  tau_run <- apply(arr, 2, function(a) geyer_tau(a))
  out <- NULL
  for (L in c(1, 10, 25, 50, 99)) {
    pr <- matrix(NA, 4, d)
    for (b in 1:4) {
      idx <- (b - 1) * 100 + 1:100 + 100 * 0
      idx <- idx[idx <= n]
      blk <- lapply(G, function(g) g[idx, , drop = FALSE])
      V <- apply(do.call(rbind, blk), 2, var)
      T <- Reduce(`+`, lapply(blk, function(g) colMeans((g[-(1:L), , drop = FALSE] - g[1:(nrow(g) - L), , drop = FALSE])^2))) / length(blk) / (L * V)
      rho <- pmax(0, 1 - T * L / 2)^(1 / L)
      pr[b, ] <- (1 + rho) / (1 - pmin(rho, 1 - 1e-9))
    }
    tau_pred <- apply(pr, 2, median)
    slow <- order(-tau_run)[1:8]
    out <- rbind(out, data.frame(L = L, med_ratio_all = median(tau_pred / tau_run),
      med_ratio_slow8 = median(tau_pred[slow] / tau_run[slow]),
      spearman = cor(tau_pred, tau_run, method = "spearman")))
  }
  cat(basename(f), " run tau slowest 8:", round(sort(tau_run, TRUE)[1:8]), "\n"); print(out, digits = 3)
}
