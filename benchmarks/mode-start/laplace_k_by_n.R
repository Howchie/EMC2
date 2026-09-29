# Does the Laplace Pareto k depend on the number of importance draws?
# (It does not, from 256 up: it measures the approximation.)
#   Rscript benchmarks/mode-start/laplace_k_by_n.R
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1]))); library(parallel)
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("EMC2"); for (f in ls(ns, all.names = TRUE, pattern = "^\\.emc_")) assign(f, get(f, ns))
for (case in c("corr", "wide", "lba71")) {
  f <- file.path(out, paste0("preburn_", case, ".rds"))
  sm <- EMC2:::restore_duplicates(readRDS(f))[[1]]; pn <- rownames(sm$samples$alpha)
  res <- do.call(rbind, mclapply(seq_len(min(sm$n_subjects, 15)), function(s) {
    set.seed(s); r <- .emc_burn_subject_mode(sm, s); if (!is.matrix(r$covariance)) return(NULL)
    prior <- .emc_burn_subject_prior(sm, s)
    ll <- function(x) .emc_burn_ll_batch(x, sm$data[[s]], sm$model, pn)
    sapply(c(64, 256, 1024, 4096), function(n) { set.seed(n + s)
      .emc_laplace_is_start(r$alpha, r$covariance, prior, ll, n_draws = n)$pareto_k })
  }, mc.cores = 15))
  colnames(res) <- paste0("n=", c(64, 256, 1024, 4096))
  cat("==", case, "p =", sm$n_pars, " subjects:", nrow(res), "\n")
  print(round(rbind(median = apply(res, 2, median), q90 = apply(res, 2, quantile, .9),
    `k>0.7` = colSums(res > 0.7), `k>1` = colSums(res > 1)), 2))
  cat("per-subject spread of k across n (median sd):", round(median(apply(res, 1, sd)), 2), "\n")
}
