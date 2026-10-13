# Single-core time of the pool likelihood for k fresh members on a few
# subjects of a fixture (rows = the subject's own recent posterior draws).
#   EMC_LIB=<lib> Rscript pool_cost.R <fit.rds>
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
fx <- commandArgs(TRUE)[1]
emc <- readRDS(fx)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
ch <- emc[[1]]
n <- ch$n_subjects; p <- ch$n_pars
idx <- seq_len(ch$samples$idx)
subs <- unique(round(seq(1, n, length.out = min(n, 6))))
ks <- c(1L, 7L, 31L, 105L)
tm <- sapply(ks, function(k) {
  mean(vapply(subs, function(s) {
    rows <- t(ch$samples$alpha[, s, tail(idx, k)])
    if (k == 1L) rows <- matrix(rows, 1)
    colnames(rows) <- ch$par_names
    reps <- max(3L, ceiling(200 / k))
    t0 <- proc.time()[["elapsed"]]
    for (r in seq_len(reps)) ll <- ns$.emc_ensemble_ll_candidate(rows, ch$data[[s]], ch$model)
    (proc.time()[["elapsed"]] - t0) / reps
  }, 0))
})
names(tm) <- ks
fixed <- tm[["1"]]; per_row <- (tm[["105"]] - tm[["7"]]) / 98
cat(sprintf("%s: n %d p %d | ms per subject call: %s | fixed %.3f ms + %.4f ms per row\n",
            basename(fx), n, p, paste(sprintf("k=%s %.3f", ks, 1000 * tm), collapse = ", "),
            1000 * fixed, 1000 * per_row))
cat(sprintf("  pool CPU per chain iteration: M=8 %.1f ms, M=32 %.1f ms, M=n=%d %.1f ms\n",
            1000 * n * (fixed + 6 * per_row), 1000 * n * (fixed + 30 * per_row), n,
            1000 * n * (fixed + (n - 2) * per_row)))
