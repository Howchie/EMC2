# Posterior-level comparison for ordinary adaptive allocator runs.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_pooled_posterior.R results-dir')
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
seeds <- as.integer(strsplit(Sys.getenv(
  "POOLED_SEEDS", "101 202 303 404 505 606 707"), "\\s+")[[1L]])
quantiles <- c(.05, .5, .95)
collect <- function(x) {
  chains <- lapply(x$fit, function(ch) {
    a <- aperm(ch$samples$alpha, c(3L, 1L, 2L))
    dim(a) <- c(dim(a)[1L], dim(a)[2L] * dim(a)[3L])
    colnames(a) <- unlist(lapply(dimnames(ch$samples$alpha)[[2L]], function(s)
      paste(dimnames(ch$samples$alpha)[[1L]], s, sep = '__')), use.names = FALSE)
    pooled <- t(ch$samples$pooled)
    colnames(pooled) <- rownames(ch$samples$pooled)
    cbind(a, pooled)
  })
  vars <- colnames(chains[[1L]])
  draws <- array(NA_real_, c(nrow(chains[[1L]]), length(chains), length(vars)),
                 dimnames = list(NULL, NULL, vars))
  for (j in seq_along(chains)) draws[, j, ] <- chains[[j]][, vars, drop = FALSE]
  summarise_draws(as_draws_array(draws), 'mean', 'mcse_mean', 'ess_bulk', 'ess_tail', 'rhat',
    ~quantile2(.x, probs = quantiles),
    ~setNames(vapply(quantiles, function(p) mcse_quantile(.x, probs = p), numeric(1)),
              c('mcse_q05', 'mcse_q50', 'mcse_q95')))
}
ratios <- matrix(NA_real_, length(seeds), 2L,
                 dimnames = list(seeds, c('bulk', 'tail')))
time_speedups <- numeric(length(seeds))
for (i in seq_along(seeds)) {
  seed <- seeds[i]
  a <- readRDS(file.path(root, if (seed == 101) 'pooled_base.rds' else paste0('pooled_', seed, '_base.rds')))
  b <- readRDS(file.path(root, if (seed == 101) 'pooled_new.rds' else paste0('pooled_', seed, '_new.rds')))
  stopifnot(identical(a$data_md5, b$data_md5))
  time_speedups[i] <- a$timing[['elapsed']] / b$timing[['elapsed']]
  sa <- collect(a); sb <- collect(b)
  stopifnot(identical(sa$variable, sb$variable))
  draws_per_chain <- dim(a$fit[[1L]]$samples$alpha)[3L]
  total_draws <- draws_per_chain * length(a$fit)
  for (kind in c('bulk', 'tail')) {
    col <- paste0('ess_', kind)
    ratios[i, kind] <- min(sb[[col]]) / min(sa[[col]])
  }
  diffs <- cbind(abs(sa$q5 - sb$q5) / sqrt(sa$mcse_q05^2 + sb$mcse_q05^2),
                 abs(sa$q50 - sb$q50) / sqrt(sa$mcse_q50^2 + sb$mcse_q50^2),
                 abs(sa$q95 - sb$q95) / sqrt(sa$mcse_q95^2 + sb$mcse_q95^2))
  colnames(diffs) <- c('q05', 'q50', 'q95')
  cat('\nseed', seed, 'draws per chain', draws_per_chain,
      'elapsed base/candidate', a$timing[['elapsed']], b$timing[['elapsed']],
      'data hash', a$data_md5, '\n')
  cat('minimum bulk ESS/draw base/candidate', min(sa$ess_bulk) / total_draws,
      min(sb$ess_bulk) / total_draws, 'minimum tail ESS/draw', min(sa$ess_tail) / total_draws,
      min(sb$ess_tail) / total_draws, '\n')
  cat('minimum bulk ESS/s base/candidate', min(sa$ess_bulk) / a$timing[['elapsed']],
      min(sb$ess_bulk) / b$timing[['elapsed']], 'minimum tail ESS/s',
      min(sa$ess_tail) / a$timing[['elapsed']], min(sb$ess_tail) / b$timing[['elapsed']], '\n')
  cat('quantile differences beyond 4 combined MCSE:', sum(diffs > 4), '/', length(diffs),
      'largest standardized difference', max(diffs), '\n')
  cat('largest quantile standardized differences by variable:\n')
  print(data.frame(variable = sa$variable, max_mcse = apply(diffs, 1L, max)), row.names = FALSE)
  cat('candidate diagnostics:\n'); print(sb, width = Inf)
}
cat('\nPer-seed minimum ESS ratios (candidate/base):\n'); print(ratios)
for (kind in colnames(ratios)) {
  z <- log(ratios[, kind]); ci <- exp(mean(z) + c(-1, 1) * qt(.975, df = length(z) - 1) * sd(z) / sqrt(length(z)))
  cat(kind, 'geometric mean ratio', exp(mean(z)), '95% t interval', ci, '\n')
}
z <- log(time_speedups)
time_ci <- exp(mean(z) + c(-1, 1) * qt(.975, df = length(z) - 1) * sd(z) / sqrt(length(z)))
cat('elapsed-time geometric speedup', exp(mean(z)), '95% t interval', time_ci,
    'mean paired elapsed reduction percent',
    100 * mean(1 - 1 / time_speedups), '\n')
# Same-code repeat is a useful control for timing-driven allocator variation.
repeat_files <- file.path(root, c('pooled_base.rds', 'pooled_audit_base_repeat.rds'))
if (all(file.exists(repeat_files))) {
  a <- readRDS(repeat_files[[1L]]); r <- readRDS(repeat_files[[2L]])
  sa <- collect(a); sr <- collect(r)
  stopifnot(identical(a$data_md5, r$data_md5))
  repeat_diff <- cbind(abs(sa$q5 - sr$q5) / sqrt(sa$mcse_q05^2 + sr$mcse_q05^2),
                       abs(sa$q50 - sr$q50) / sqrt(sa$mcse_q50^2 + sr$mcse_q50^2),
                       abs(sa$q95 - sr$q95) / sqrt(sa$mcse_q95^2 + sr$mcse_q95^2))
  cat('\nSame-baseline seed-101 repeat: data hash', a$data_md5,
      'draws identical', identical(a$fit[[1]]$samples$alpha, r$fit[[1]]$samples$alpha),
      'max quantile difference in combined MCSE', max(repeat_diff),
      'quantile checks beyond 4 combined MCSE', sum(repeat_diff > 4), '/', length(repeat_diff), '\n')
}
