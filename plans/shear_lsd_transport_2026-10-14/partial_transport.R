# Partial transport: a shear step that carries only a fixed subset S of
# subjects by their Gaussian frame map and holds the rest. S is ranked by
# the prior's share of each subject's conditional precision on the moved
# parameter at the block reference (fixed per block, so the map depends on
# the group states only). Reports the Gaussian-predicted width of the orbit
# and the real-likelihood acceptance of jumps, by |S|.
#
# Rscript partial_transport.R <emc.rds> <chain> lsd:<par> <draws|auto> [sizes=0,10,20,40,80] [NC=8]
suppressPackageStartupMessages(library(EMC2, lib.loc = if (nzchar(Sys.getenv("EMC_LIB"))) Sys.getenv("EMC_LIB") else NULL))
ns <- asNamespace("EMC2")
a <- commandArgs(TRUE)
kv <- function(key, def) { m <- grep(paste0("^", key, "="), a, value = TRUE); if (length(m)) sub(".*=", "", m[1]) else def }
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(a[1])), quiet = TRUE)
j <- as.integer(a[2]); what <- sub(".*:", "", a[3]); NC <- as.integer(kv("NC", 8))
sizes <- as.integer(strsplit(kv("sizes", "0,10,20,40,80"), ",")[[1]])
ch0 <- emc[[j]]; p <- ch0$n_pars; n <- ch0$n_subjects; pn <- ch0$par_names; last <- ch0$samples$idx
draws <- if (a[4] == "auto") round(seq(last - 7 * 60, last, by = 60)) else as.integer(strsplit(a[4], ",")[[1]])
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample"); idx <- idx[idx <= last]
ctx <- ns$.emc_shear_context(ch0); H <- ctx$packed$cube; LIN <- ctx$packed$lin; spec <- ctx$spec; prior <- ch0$prior
l <- if (what %in% pn) match(what, pn) else as.integer(what)
wsd <- sd(0.5 * log(ch0$samples$theta_var[l, l, idx]))
state_at <- function(mu, V, u) { d <- rep(1, p); d[l] <- exp(u); list(mu = mu, V = V * outer(d, d)) }
jac_u <- function(u) (spec$size[l] + 1) * u
lp_group <- function(st) ns$.emc_shear_log_prior_mu(st$mu, prior) + ns$.emc_shear_log_prior_v(st$V, spec)
V_ref <- apply(ch0$samples$theta_var[, , idx, drop = FALSE], 1:2, mean); Vi_ref <- solve(V_ref)
share <- vapply(seq_len(n), function(s) solve(H[, , s] + Vi_ref, Vi_ref)[l, l], 1)
ord <- order(share, decreasing = TRUE)
cat(sprintf("%s chain %d log sd(%s), window sd %.3f; prior share quantiles %s\n", basename(a[1]), j, pn[l], wsd,
            paste(sprintf("%.2f", quantile(share, c(0, .25, .5, .75, 1))), collapse = " ")))
ug <- seq(-4, 4, by = 0.25) * wsd
sd_of <- function(lg) { w <- exp(lg - max(lg)); w <- w / sum(w); sqrt(sum(w * ug^2) - sum(w * ug)^2) / wsd }
frames <- function(st) { Vi <- solve(st$V)
  lapply(seq_len(n), function(s) { Q <- H[, , s] + Vi; R <- chol(Q); list(m = drop(solve(Q, LIN[, s] + Vi %*% st$mu)), R = R, ld = sum(log(diag(R)))) }) }
logN <- function(A, st) { R <- chol(st$V); z <- backsolve(R, A - st$mu, transpose = TRUE); -colSums(z^2) / 2 - sum(log(diag(R))) }
llb <- function(X, s) ns$.emc_burn_ll_batch(X, ch0$data[[s]], ch0$model, pn)
one_row <- function(A, ss) unlist(parallel::mclapply(ss, function(s) { X <- matrix(A[, s], 1); colnames(X) <- pn; llb(X, s) }, mc.cores = NC))
jumps <- c(-1, -0.5, 0.5, 1)
res <- array(NA, c(length(draws), length(sizes), 5), dimnames = list(draws, sizes, c("width", jumps)))
for (di in seq_along(draws)) { d <- draws[di]
  alpha <- matrix(ch0$samples$alpha[, , d], p); mu <- ch0$samples$theta_mu[, d]; V <- ch0$samples$theta_var[, , d]
  st0 <- state_at(mu, V, 0); f0 <- frames(st0); ll0 <- one_row(alpha, seq_len(n))
  sts <- lapply(jumps * wsd, function(u) state_at(mu, V, u)); fs <- lapply(sts, frames)
  for (si in seq_along(sizes)) { S <- ord[seq_len(sizes[si])]
    # predicted width: subjects in S integrated out under their frames, others held
    res[di, si, 1] <- sd_of(vapply(ug, function(u) { st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
      if (!is.finite(g)) return(-Inf); Vi <- solve(st$V); ldV <- determinant(st$V)$modulus; Vm <- drop(Vi %*% st$mu)
      held <- setdiff(seq_len(n), S); out <- g + jac_u(u) + if (length(held)) sum(logN(alpha[, held, drop = FALSE], st)) else 0
      out + sum(vapply(S, function(s) { Q <- H[, , s] + Vi; h <- LIN[, s] + Vm
        -0.5 * ldV - 0.5 * determinant(Q)$modulus + 0.5 * sum(h * solve(Q, h)) - 0.5 * sum(st$mu * Vm) }, 1)) }, 1))
    for (ji in seq_along(jumps)) {
      st1 <- sts[[ji]]; f1 <- fs[[ji]]; A1 <- alpha; lj <- 0
      for (s in S) { A1[, s] <- f1[[s]]$m + drop(backsolve(f1[[s]]$R, f0[[s]]$R %*% (alpha[, s] - f0[[s]]$m))); lj <- lj + f0[[s]]$ld - f1[[s]]$ld }
      dll <- if (length(S)) sum(one_row(A1, S) - ll0[S]) else 0
      lr <- lp_group(st1) + jac_u(jumps[ji] * wsd) - lp_group(st0) + sum(logN(A1, st1)) - sum(logN(alpha, st0)) + lj + dll
      res[di, si, 1 + ji] <- min(1, exp(lr))
    }
  }
  cat(sprintf("draw %4d: %s\n", d, paste(sprintf("|S|=%d width %.2f acc %s", sizes, res[di, , 1],
    apply(res[di, , -1, drop = FALSE], 2, function(x) paste(sprintf("%.2f", x), collapse = "/"))), collapse = " | ")))
}
cat("means over draws (width, then acceptance at -1/-0.5/+0.5/+1):\n")
for (si in seq_along(sizes)) cat(sprintf("  |S|=%2d: width %.2f, acc %s\n", sizes[si], mean(res[, si, 1]), paste(sprintf("%.2f", colMeans(res[, si, -1])), collapse = "/")))
