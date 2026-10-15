# Attribution of persistent transport to move families, block by block, from an instrumented continuation.
r <- readRDS(commandArgs(TRUE)[1]); L <- 25; lagd <- 10; ndir <- 3
r$plans <- lapply(r$plans, function(h) { it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] })
ridge_mu <- readRDS(file.path(dirname(commandArgs(TRUE)[1]), "ridge_u.rds"))
FAM <- c("particle", "group", "shear", "collapsed", "hmc")
p <- ncol(r$alpha_mean[[1]]); nch <- length(r$draws); d3 <- 3 * p
Zst <- lapply(1:nch, function(j) cbind(r$draws[[j]], r$alpha_mean[[j]]))    # end-of-iteration state statistic
its <- vapply(r$plans[[1]], `[[`, numeric(1), "iteration"); n0 <- min(its) - 100
blocks <- lapply(seq_along(r$plans[[1]]), function(b) r$plans[[1]][[b]])
paths <- lapply(1:nch, function(j) { h <- r$plans[[j]]; c(lapply(h, `[[`, "path"), list(r$last_path[[j]])) })
ends <- c(its, n0 + nrow(Zst[[1]]))
starts <- c(n0 + 1, head(ends, -1) + 1)
famcost <- function(h) c(particle = sum(h$particle_rate * h$c_part, na.rm = TRUE),
  group = h$C_rest, shear = sum(h$rounds) * h$c_round + h$shear_fixed,
  collapsed = h$collapsed_rounds * h$collapsed_cost, hmc = sum(h$rate * h$c_hmc, na.rm = TRUE) + h$hmc_fixed)
for (b in seq_along(ends)) {
  rows <- (starts[b]:ends[b]) - n0
  P <- lapply(paths, `[[`, b); if (any(vapply(P, is.null, logical(1)))) next
  if (any(vapply(P, nrow, integer(1)) != length(rows))) { cat("block", b, "path rows", vapply(P, nrow, integer(1)), "vs", length(rows), "\n"); next }
  X <- lapply(Zst, function(z) z[rows, ])
  sdw <- sqrt(colMeans(t(sapply(X, function(x) apply(x, 2, var)))))
  Xs <- lapply(X, function(x) sweep(x, 2, sdw, "/")); all <- do.call(rbind, Xs); m <- colMeans(all)
  C0 <- crossprod(sweep(all, 2, m)) / nrow(all); C0 <- C0 + diag(1e-3 * mean(diag(C0)), d3)
  Cl <- Reduce(`+`, lapply(Xs, function(x) { n <- nrow(x); a <- sweep(x[1:(n - lagd), ], 2, m); bb <- sweep(x[-(1:lagd), ], 2, m); crossprod(a, bb) })) / (nch * (nrow(Xs[[1]]) - lagd))
  Cl <- (Cl + t(Cl)) / 2
  e0 <- eigen(C0, symmetric = TRUE); W <- e0$vectors %*% diag(1 / sqrt(e0$values)) %*% t(e0$vectors)
  ev <- eigen(W %*% Cl %*% W, symmetric = TRUE); V <- W %*% ev$vectors[, 1:ndir]   # v' C0 v = 1, scaled coords
  V <- V / sdw                                                                        # raw-unit directions
  rv <- c(ridge_mu, numeric(p), ridge_mu); rv <- rv / sqrt(drop(t(rv) %*% diag(sdw) %*% C0 %*% diag(sdw) %*% rv)); V <- cbind(rv, V)
  ev$values <- c(NA, ev$values)
  h <- blocks[[min(b, length(blocks))]]; hp <- if (b > 1) blocks[[b - 1]] else NULL
  cat(sprintf("\n== block %d (iters %d-%d) plan from rebuild %s: rounds %.2f hmc rate %.2f particle %.2f collapsed %.2f\n", b, starts[b], ends[b],
      if (is.null(hp)) "pre-continuation" else hp$iteration, if (is.null(hp)) NA else sum(hp$rounds), if (is.null(hp)) NA else mean(hp$rate),
      if (is.null(hp)) NA else mean(hp$particle_rate), if (is.null(hp)) NA else hp$collapsed_rounds))
  cost <- if (!is.null(hp)) famcost(hp) else rep(NA, 5)
  for (k in 1:(ndir + 1)) {
    v <- V[, k]; nm <- c(paste0("mu.", 1:p), paste0("lsd.", 1:p), paste0("abar.", 1:p))
    top <- order(-abs(v * sdw))[1:4]
    Tf <- E1 <- numeric(5)
    for (j in 1:nch) { Pj <- P[[j]]; n <- nrow(Pj)
      dfam <- sapply(1:5, function(f) Pj[, (f - 1) * d3 + 1:d3] %*% v)              # n x F
      Dl <- rowSums(dfam); wts <- 1 - abs(-(L - 1):(L - 1)) / L
      sm <- sapply(1:n, function(t) { l <- -(L - 1):(L - 1); ok <- t + l >= 1 & t + l <= n; sum(wts[ok] * Dl[t + l[ok]]) })
      Tf <- Tf + colMeans(dfam * sm) / nch; E1 <- E1 + colMeans(dfam^2) / nch }
    cat(sprintf("  dir %d (lag-%d acf %.2f; %s)\n", k, lagd, ev$values[k], paste(sprintf("%s %+.2f", nm[top], (v * sdw)[top] / sqrt(sum((v * sdw)^2))), collapse = " ")))
    cat(sprintf("    %-10s %s\n", "", paste(sprintf("%10s", FAM), collapse = "")))
    cat(sprintf("    %-10s %s   total %.4f\n", "one-step", paste(sprintf("%10.4f", E1), collapse = ""), sum(E1)))
    cat(sprintf("    %-10s %s   total %.4f\n", "persist", paste(sprintf("%10.4f", Tf), collapse = ""), sum(Tf)))
    cat(sprintf("    %-10s %s\n", "per cost", paste(sprintf("%10.4f", Tf / cost), collapse = "")))
  }
  cat("  cost/iter by family:", paste(sprintf("%s %.2f", FAM, cost), collapse = "  "), "\n")
}
