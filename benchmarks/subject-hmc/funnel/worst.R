suppressPackageStartupMessages(library(posterior))
f <- commandArgs(TRUE)[1]; e <- readRDS(f); p <- e[[1]]$n_pars; n <- e[[1]]$n_subjects; pn <- e[[1]]$par_names
arr <- lapply(e, function(ch) { w <- which(ch$samples$stage == "sample")
  t(sapply(w, function(t) c(ch$samples$theta_mu[, t], 0.5 * log(diag(as.matrix(ch$samples$theta_var[, , t]))), as.vector(ch$samples$alpha[, , t])))) })
m <- min(sapply(arr, nrow)); a <- simplify2array(lapply(arr, function(x) x[1:m, ]))
nm <- c(paste0("mu_", pn), paste0("lsd_", pn), paste0("a", rep(1:n, each = p), "_", pn))
r <- apply(a, 2, rhat); eb <- apply(a, 2, ess_bulk)
o <- order(eb)[1:10]
blk <- cut(seq_len(m), 6, labels = FALSE)
for (k in o) { cat(sprintf("%-22s rhat %.3f ess %4.0f | block means by chain:\n", nm[k], r[k], eb[k]))
  print(round(sapply(1:dim(a)[3], function(c) tapply(a[, k, c], blk, mean)), 2)) }
