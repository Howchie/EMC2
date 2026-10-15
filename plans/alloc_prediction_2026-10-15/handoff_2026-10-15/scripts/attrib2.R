# Persistent-transport attribution per block of a full instrumented fit (saved emc with paths in the allocator history).
args <- commandArgs(TRUE); lib <- args[2]
suppressMessages(library(EMC2, lib.loc = lib)); ns <- asNamespace("EMC2")
emc <- ns$restore_duplicates(readRDS(args[1])); nch <- length(emc)
L <- 25; lagd <- 10; ndir <- 3; FAM <- ns$.EMC_ALLOC_PATH_FAMILIES; nf <- length(FAM)
p <- emc[[1]]$n_pars; d3 <- 3 * p
state <- function(ch, k) t(sapply(k, function(t) ns$.emc_alloc_state(ch$samples$theta_mu[, t], ch$samples$theta_var[, , t], ch$samples$alpha[, , t])))
hist <- lapply(emc, function(ch) { h <- ch$alloc$history; it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] })
its <- Reduce(intersect, lapply(hist, function(h) vapply(h, `[[`, numeric(1), "iteration")))
famcost <- function(h) c(particle = sum(h$particle_rate * h$c_part, na.rm = TRUE), group = h$C_rest,
  shear = sum(h$rounds) * h$c_round + h$shear_fixed, collapsed = h$collapsed_rounds * h$collapsed_cost,
  hmc = sum(h$rate * h$c_hmc, na.rm = TRUE) + h$hmc_fixed)
out <- NULL
for (R in its) {
  hs <- lapply(hist, function(h) Filter(function(x) x$iteration == R, h)[[1]])
  P <- lapply(hs, `[[`, "path"); if (any(vapply(P, is.null, logical(1)))) next
  n <- min(vapply(P, nrow, integer(1))); if (n < 50) next
  P <- lapply(P, function(x) tail(x, n)); rows <- (R - n + 1):R
  if (R > emc[[1]]$samples$idx || min(rows) < 2) next
  X <- lapply(emc, state, k = rows)
  X0 <- lapply(emc, state, k = rows - 1L)
  ok <- all(vapply(1:nch, function(j) max(abs(Reduce(`+`, lapply(1:nf, function(f) P[[j]][, (f - 1) * d3 + 1:d3])) - (X[[j]] - X0[[j]]))) < 1e-8, logical(1)))
  if (!ok) { cat("skip block", R, "(path does not match stored states)
"); next }
  sdw <- sqrt(colMeans(t(sapply(X, function(x) apply(x, 2, var))))); sdw[!is.finite(sdw) | sdw <= 0] <- 1
  Xs <- lapply(X, function(x) sweep(x, 2, sdw, "/")); all <- do.call(rbind, Xs); m <- colMeans(all)
  C0 <- crossprod(sweep(all, 2, m)) / nrow(all); C0 <- C0 + diag(1e-3 * mean(diag(C0)), d3)
  Cl <- Reduce(`+`, lapply(Xs, function(x) crossprod(sweep(x[1:(n - lagd), ], 2, m), sweep(x[-(1:lagd), ], 2, m)))) / (nch * (n - lagd)); Cl <- (Cl + t(Cl)) / 2
  e0 <- eigen(C0, symmetric = TRUE); W <- e0$vectors %*% diag(1 / sqrt(e0$values)) %*% t(e0$vectors)
  ev <- eigen(W %*% Cl %*% W, symmetric = TRUE); V <- (W %*% ev$vectors[, 1:ndir]) / sdw
  prev <- Filter(function(x) x$iteration < R, hist[[1]]); cost <- if (length(prev)) famcost(prev[[length(prev)]]) else rep(NA, nf)
  stage <- emc[[1]]$samples$stage[R]
  for (k in 1:ndir) { v <- V[, k]; Tf <- E1 <- numeric(nf); vg <- 0
    for (j in 1:nch) { dfam <- sapply(1:nf, function(f) P[[j]][, (f - 1) * d3 + 1:d3] %*% v); Dl <- rowSums(dfam)
      wts <- 1 - abs(-(L - 1):(L - 1)) / L
      sm <- sapply(1:n, function(t) { l <- -(L - 1):(L - 1); ok <- t + l >= 1 & t + l <= n; sum(wts[ok] * Dl[t + l[ok]]) })
      Tf <- Tf + colMeans(dfam * sm) / nch; E1 <- E1 + colMeans(dfam^2) / nch
      x <- drop(X[[j]] %*% v); vg <- vg + mean((x[-(1:L)] - x[1:(n - L)])^2) / L / nch }
    out <- rbind(out, data.frame(R = R, stage = stage, dir = k, acf10 = ev$values[k], vg = vg, tot = sum(Tf),
      t(setNames(Tf, paste0("T.", FAM))), t(setNames(E1, paste0("E.", FAM))), t(setNames(cost, paste0("c.", FAM))))) }
}
saveRDS(out, sub("\\.rds$", "_attrib.rds", args[1]))
cat("identity: sum of family transport vs lag-25 variogram / 25, correlation", round(cor(out$tot, out$vg), 3), " median ratio", round(median(out$tot / out$vg), 2), "\n")
for (st in unique(out$stage)) for (k in 1:ndir) { x <- out[out$stage == st & out$dir == k, ]; if (!nrow(x)) next
  sh <- colMeans(x[, paste0("T.", FAM)]); e <- colMeans(x[, paste0("E.", FAM)]); cc <- colMeans(x[, paste0("c.", FAM)], na.rm = TRUE)
  cat(sprintf("%-7s dir %d (%2d blocks, lag-10 acf %.2f) share of persistent transport: %s | persist/one-step: %s | per cost: %s\n", st, k, nrow(x), mean(x$acf10),
    paste(sprintf("%s %+.2f", FAM, sh / sum(sh)), collapse = " "), paste(sprintf("%.2f", sh / pmax(e, 1e-12)), collapse = " "),
    paste(sprintf("%.3f", sh / cc), collapse = " "))) }
