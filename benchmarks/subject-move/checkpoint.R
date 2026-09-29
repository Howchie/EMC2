# Inspect a completed stage without running another fit.
# Rscript checkpoint.R <library> <checkpoint.RData> [fixture.rds]
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = args[1]))
load(args[2])
emc <- EMC2:::restore_duplicates(emc)
cat("EMC2 from", find.package("EMC2"), "\n")
for (j in seq_along(emc)) {
  ch <- emc[[j]]; a <- ch$samples$alpha
  k <- seq.int(max(1L, dim(a)[3] - 249L), dim(a)[3])
  unique_states <- apply(a[1, , k, drop = FALSE], 2, function(x) length(unique(as.numeric(x))))
  moved <- apply(a[, , k, drop = FALSE], 2, function(x)
    mean(colSums(abs(x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE])) > 0))
  pm <- attr(ch$samples, "pm_settings")
  cat("chain", j, "stage", tail(ch$samples$stage, 1), "draws", dim(a)[3], "\n")
  cat(" unique states:", quantile(unique_states), "\n")
  cat(" movement:", signif(quantile(moved), 3), "\n")
  eps <- do.call(rbind, lapply(pm, function(s) s[[1]]$epsilon))
  cat(" epsilon quantiles by proposal:\n"); print(apply(eps, 2, quantile))
  block_scale <- unlist(lapply(pm, function(s) {
    bw <- s[[1]]$block_warmup
    if (is.null(bw)) numeric(0) else bw$scale
  }))
  if (length(block_scale)) {
    cat(" subject block proposal scale quantiles:\n"); print(quantile(block_scale))
  }
  mala_scale <- unlist(lapply(pm, function(s) {
    mt <- s[[1]]$mala
    if (is.null(mt)) numeric(0) else mt$scale
  }))
  if (length(mala_scale)) {
    cat(" subject MALA step scale quantiles:\n"); print(quantile(mala_scale))
  }
  if (!is.null(ch$chains_var)) {
    v <- sapply(ch$chains_var, function(v) eigen(v, symmetric = TRUE, only.values = TRUE)$values)
    cat(" covariance min eigenvalues:", signif(quantile(apply(v, 2, min)), 3), "\n")
    cat(" covariance condition numbers:", signif(quantile(apply(v, 2, function(x) max(x)/min(x))), 3), "\n")
  }
}
recent <- min(250L, min(vapply(emc, function(ch) dim(ch$samples$alpha)[3], 0L)))
rh <- sapply(seq_len(emc[[1]]$n_subjects), function(s)
  vapply(seq_len(emc[[1]]$n_pars), function(p) posterior::rhat(
    vapply(emc, function(ch) tail(ch$samples$alpha[p, s, ], recent), numeric(recent))), 0))
cat("Recent-window subject R-hat (warm-up values are descriptive):\n")
print(quantile(rh, na.rm = TRUE))
if (length(args) >= 3L) {
  fx <- readRDS(args[3]); truth <- t(as.matrix(fx$subject_effects))
  ch <- emc[[1]]
  if (identical(dim(truth), dim(ch$samples$alpha)[1:2])) {
    ll_truth <- vapply(seq_len(ch$n_subjects), function(s)
      as.numeric(EMC2:::calc_ll_manager(matrix(truth[, s], nrow = 1,
        dimnames = list(NULL, rownames(ch$samples$alpha))),
        dadm = ch$data[[s]], model = ch$model)), 0)
    cat("log likelihood at generating truth minus current state (positive = current worse):\n")
    print(sapply(emc, function(ch) quantile(ll_truth - ch$samples$subj_ll[, ch$samples$idx])))
  }
}
