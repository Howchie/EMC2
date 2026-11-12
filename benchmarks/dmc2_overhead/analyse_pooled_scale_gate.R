args <- commandArgs(TRUE)
if (length(args) < 2L) stop('usage: analyse_pooled_scale_gate.R results-dir seed [seed ...]')
root <- args[[1L]]
seeds <- as.integer(args[-1L])
rows <- list()
for (seed in seeds) {
  runs <- lapply(c('base', 'new'), function(arm) {
    readRDS(file.path(root, sprintf('pooled_%d_%s_stop.rds', seed, arm)))
  })
  names(runs) <- c('base', 'new')
  stopifnot(identical(runs$base$data_md5, runs$new$data_md5))
  for (arm in names(runs)) {
    x <- runs[[arm]]
    within_limit <- isTRUE(x$reached) &&
      !is.null(x$target) && x$target$sample_draws <= 4000L
    target <- if (within_limit) x$target else NULL
    eligible_checks <- Filter(function(m) m$sample_draws <= 4000L, x$checks)
    last_check <- if (length(eligible_checks)) eligible_checks[[length(eligible_checks)]] else NULL
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, arm = arm, reached_within_default = within_limit,
      sample_draws = if (is.null(target)) NA_integer_ else target$sample_draws,
      censored_at_draws = if (within_limit || is.null(last_check)) NA_integer_ else last_check$sample_draws,
      elapsed_seconds = if (within_limit) target$elapsed_seconds else
        if (is.null(last_check)) NA_real_ else last_check$elapsed_seconds,
      alpha_rhat_q99 = if (is.null(last_check)) NA_real_ else last_check$alpha_rhat_q99,
      mu_rhat_max = if (is.null(last_check)) NA_real_ else last_check$mu_rhat_max,
      sigma2_rhat_max = if (is.null(last_check)) NA_real_ else last_check$sigma2_rhat_max,
      group_log_sd_ess_bulk_min = if (is.null(last_check)) NA_real_ else last_check$group_log_sd_ess_bulk_min
    )
  }
}
out <- do.call(rbind, rows)
cat('Shared scale-aware target patterned after DMC2 hardfit (full prefixes; no discard, gap, spread, or culling; three chains):\n')
print(out, row.names = FALSE)
paired <- split(out, out$seed)
paired <- paired[vapply(paired, function(x) all(x$reached_within_default), logical(1))]
if (!length(paired)) {
  cat('No seed pair reached the target in both arms within 4,000 draws; time-to-target speedup is censored.\n')
} else {
  speedup <- vapply(paired, function(x)
    x$elapsed_seconds[x$arm == 'base'] / x$elapsed_seconds[x$arm == 'new'], 0)
  cat('paired time-to-target speedups (base/new):\n')
  print(speedup)
  if (length(speedup) >= 2L) {
    z <- log(speedup)
    ci <- exp(mean(z) + c(-1, 1) * qt(.975, df = length(z) - 1) * sd(z) / sqrt(length(z)))
    cat('geometric mean speedup', exp(mean(z)), '95% t interval', ci, '\n')
  } else cat('geometric mean speedup', exp(mean(log(speedup))), '(one paired seed)\n')
}
