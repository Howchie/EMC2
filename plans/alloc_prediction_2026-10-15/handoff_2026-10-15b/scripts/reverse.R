# How much of a shear displacement along the ridge do the following iterations' subject and group updates reverse?
# From the instrumented PM arm's per-family paths (state statistic: mu, lsd, abar per family per iteration).
A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/handoff_2026-10-15"
r <- readRDS(file.path(A, "arms/arm_path_2.rds")); u <- readRDS(file.path(A, "scripts/ridge_u.rds"))
FAM <- c("particle", "group", "shear", "collapsed", "hmc"); p <- length(u); d3 <- 3 * p; nch <- length(r$draws)
r$plans <- lapply(r$plans, function(h) { it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] })
paths <- lapply(1:nch, function(j) do.call(rbind, c(lapply(r$plans[[j]], `[[`, "path"), list(r$last_path[[j]]))))
wsd <- sqrt(mean(sapply(r$draws, function(g) var(drop(g[, 1:p] %*% u)))))   # within-chain sd along u (group mean)
for (coord in c("group mean", "collective")) {
  v <- if (coord == "group mean") c(u, numeric(2 * p)) else c(u / 2, numeric(p), u / 2)
  cat(sprintf("\n== ridge coordinate: %s (units: within-chain sd of the group mean along u)\n", coord))
  res <- NULL
  for (j in 1:nch) { P <- paths[[j]]; n <- nrow(P)
    D <- sapply(1:5, function(f) drop(P[, (f - 1) * d3 + 1:d3] %*% v)) / wsd; colnames(D) <- FAM
    sh <- D[, "shear"]; other <- D[, "particle"] + D[, "group"] + D[, "hmc"] + D[, "collapsed"]
    # regress the cumulative non-shear displacement over the next L iterations on the shear displacement at t
    for (L in c(1, 2, 5, 10, 25)) {
      cum <- sapply(1:(n - L), function(t) sum(other[(t + 1):(t + L)]))
      same <- sapply(1:(n - L), function(t) sum(other[t]))  # same-iteration (after the shear within the iteration: particle+group precede it, hmc follows)
      b <- coef(lm(cum ~ sh[1:(n - L)]))[2]
      bs <- coef(lm(sh[(1 + 1):(n - L + 1)] ~ sh[1:(n - L)]))[2]
      res <- rbind(res, data.frame(chain = j, L = L, reversed = -b, shear_next = bs, sh_var = var(sh), other_var = var(other), n = n)) } }
  agg <- aggregate(cbind(reversed, shear_next, sh_var, other_var) ~ L, res, mean)
  cat("fraction of a shear move undone by the other families over the next L iterations (mean over chains):\n"); print(round(agg, 3))
  # direct: persistent fraction of the shear move itself = cov(shear_t, total displacement over t..t+L) / var(shear_t)
  for (L in c(1, 5, 10, 25)) { fr <- sapply(1:nch, function(j) { P <- paths[[j]]; n <- nrow(P)
    D <- sapply(1:5, function(f) drop(P[, (f - 1) * d3 + 1:d3] %*% v)) / wsd; tot <- rowSums(D); sh <- D[, 3]
    cs <- sapply(1:(n - L), function(t) sum(tot[t:(t + L)])); cov(sh[1:(n - L)], cs) / var(sh[1:(n - L)]) })
    cat(sprintf("  persistent fraction of the shear move after %2d more iterations: %s (mean %.2f)\n", L, paste(round(fr, 2), collapse = " "), mean(fr))) }
  # the same for the particle family and the group update
  for (f in c("particle", "group", "hmc")) { fr <- sapply(1:nch, function(j) { P <- paths[[j]]; n <- nrow(P)
    D <- sapply(1:5, function(k) drop(P[, (k - 1) * d3 + 1:d3] %*% v)) / wsd; tot <- rowSums(D); x <- D[, match(f, FAM)]
    cs <- sapply(1:(n - 10), function(t) sum(tot[t:(t + 10)])); c(cov(x[1:(n - 10)], cs) / var(x[1:(n - 10)]), var(x)) })
    cat(sprintf("  %-8s: persistent fraction after 10 iterations %s | per-iteration variance %s\n", f, paste(round(fr[1, ], 2), collapse = " "), paste(signif(fr[2, ], 2), collapse = " "))) }
}
