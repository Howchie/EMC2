# Compare the Step 2 fit-cache probe on the analytic 50/50 bimodal target.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: analyse_multimodal_cache.R results-dir")
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
summaries <- list()
for (arm in c("base", "new")) {
  x <- readRDS(file.path(root, paste0("multimodal_", arm, ".rds")))
  indicators <- lapply(x$fit, function(ch)
    as.numeric(ch$samples$alpha[1L, 1L, ] > 0))
  n <- min(vapply(indicators, length, integer(1)))
  draws <- array(NA_real_, c(n, length(indicators), 1L),
                 dimnames = list(NULL, NULL, "right_mode"))
  for (j in seq_along(indicators)) draws[, j, 1L] <- indicators[[j]][seq_len(n)]
  stats <- summarise_draws(as_draws_array(draws), "mean", "mcse_mean",
                           "ess_bulk", "ess_tail", "rhat")
  switches <- vapply(indicators, function(z) sum(diff(z) != 0), integer(1))
  summaries[[arm]] <- list(mean = stats$mean, mcse = stats$mcse_mean,
                           stats = stats, switches = switches,
                           elapsed = x$timing[["elapsed"]], n = n)
  cat(arm, "draws/chain", n, "right-mode fraction", stats$mean,
      "MCSE", stats$mcse_mean, "bulk ESS", stats$ess_bulk,
      "tail ESS", stats$ess_tail, "Rhat", stats$rhat,
      "switches/chain", switches, "elapsed", x$timing[["elapsed"]], "\n")
  cat("distance from exact 0.5 in MCSE", abs(stats$mean - .5) / stats$mcse_mean, "\n")
}
a <- summaries$base
b <- summaries$new
cat("candidate minus baseline right-mode fraction, combined MCSE units",
    abs(a$mean - b$mean) / sqrt(a$mcse^2 + b$mcse^2), "\n")
