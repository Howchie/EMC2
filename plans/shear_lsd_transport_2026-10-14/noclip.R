# Laplace likelihood fits with negative curvature reflected (|lambda|) instead
# of clipped to zero, against clipped refits at the same state and the stored
# fits, along one log group sd (walls on, as production).
# Usage: Rscript noclip.R fixture.rds output.rds coord [chain] [draws]
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
j <- if (length(args) >= 4) as.integer(args[4]) else 1L
ch <- emc[[j]]; p <- ch$n_pars; n <- ch$n_subjects; pn <- ch$par_names
draws <- if (length(args) >= 5) as.integer(strsplit(args[5], ",")[[1]]) else c(1001, 1000, 980, 950, 900, 800, 700, 600)
k <- p + match(args[3], pn)
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
canon_clip <- ns$.emc_lik_canonical
canon_abs <- function(ll_hess, x, precision, mu) {
  e <- eigen((ll_hess + t(ll_hess)) / 2, symmetric = TRUE)
  prec <- e$vectors %*% (abs(e$values) * t(e$vectors))
  lin <- drop(prec %*% x + precision %*% (x - mu))
  if (!all(is.finite(prec)) || !all(is.finite(lin))) return(NULL)
  list(prec = prec, lin = lin, neg = sum(e$values < 0))
}
refit <- function(canon) {
  assignInNamespace(".emc_lik_canonical", canon, "EMC2")
  t0 <- proc.time()
  fits <- lapply(seq_len(n), function(s) {
    r <- tryCatch(ns$.emc_burn_subject_mode(ch, s), error = function(e) NULL)
    if (is.list(r) && is.list(r$lik)) r$lik else ch$lik_laplace[[s]]
  })
  assignInNamespace(".emc_lik_canonical", canon_clip, "EMC2")
  list(packed = ns$.emc_shear_pack(fits, p), fits = fits, cpu = sum((proc.time() - t0)[1:2]))
}
set.seed(1)
variants <- list(stored = list(packed = ns$.emc_shear_context(ch)$packed),
                 clip = refit(canon_clip), abs = refit(canon_abs))
cat("refit cpu: clip", round(variants$clip$cpu), "s, abs", round(variants$abs$cpu), "s\n")
neg <- vapply(variants$abs$fits, function(f) if (is.null(f$neg)) NA_integer_ else f$neg, 1L)
cat("subjects with negative curvature at their mode:", sum(neg > 0, na.rm = TRUE), "of", n,
    "; counts", paste(names(table(neg)), table(neg), sep = ":", collapse = " "), "\n")
grid <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2); g0 <- 5
out <- list()
for (d in draws) for (v in names(variants)) {
  c1 <- ch; c1$samples$idx <- d
  tm <- ns$.emc_shear_orbit_terms(c1, idx, k, grid, function(st) variants[[v]]$packed, 1L)
  per <- sweep(tm$per, 2, tm$per[g0, ])
  tot <- rowSums(per) + tm$lp - tm$lp[g0]
  acc <- pmin(exp(tot), 1)
  row <- data.frame(draw = d, frames = v, spread = sum(apply(per, 2, var)),
                    am1 = acc[3], amh = acc[4], aph = acc[6], ap1 = acc[7])
  print(row, digits = 3); out[[length(out) + 1]] <- row
}
res <- do.call(rbind, out)
saveRDS(list(results = res, neg = neg, head = system("git rev-parse HEAD", intern = TRUE)), args[2])
cat("\nmeans\n"); print(aggregate(cbind(am1, amh, aph, ap1) ~ frames, res, mean), digits = 3)
cat("median spread\n"); print(aggregate(spread ~ frames, res, median), digits = 3)
