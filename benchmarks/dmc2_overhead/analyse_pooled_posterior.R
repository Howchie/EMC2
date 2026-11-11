# Posterior-level comparison for ordinary adaptive allocator runs.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_pooled_posterior.R results-dir')
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
seeds <- c(101L, 202L, 303L, 404L, 505L, 606L, 707L)
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
for (i in seq_along(seeds)) {
  seed <- seeds[i]
  a <- readRDS(file.path(root, if (seed == 101) 'pooled_base.rds' else paste0('pooled_', seed, '_base.rds')))
  b <- readRDS(file.path(root, if (seed == 101) 'pooled_new.rds' else paste0('pooled_', seed, '_new.rds')))
  stopifnot(identical(a$data_md5, b$data_md5))
  sa <- collect(a); sb <- collect(b)
  stopifnot(identical(sa$variable, sb$variable))
  for (kind in c('bulk', 'tail')) {
    col <- paste0('ess_', kind)
    ratios[i, kind] <- min(sb[[col]]) / min(sa[[col]])
  }
  diffs <- cbind(abs(sa$q5 - sb$q5) / sqrt(sa$mcse_q05^2 + sb$mcse_q05^2),
                 abs(sa$q50 - sb$q50) / sqrt(sa$mcse_q50^2 + sb$mcse_q50^2),
                 abs(sa$q95 - sb$q95) / sqrt(sa$mcse_q95^2 + sb$mcse_q95^2))
  colnames(diffs) <- c('q05', 'q50', 'q95')
  cat('\nseed', seed, 'elapsed base/candidate', a$timing[['elapsed']], b$timing[['elapsed']],
      'data hash', a$data_md5, '\n')
  cat('minimum bulk ESS/draw base/candidate', min(sa$ess_bulk) / (400 * 3),
      min(sb$ess_bulk) / (400 * 3), 'minimum tail ESS/draw', min(sa$ess_tail) / 1200,
      min(sb$ess_tail) / 1200, '\n')
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
# Same-code repeat is a useful control for timing-driven allocator variation.
a <- readRDS(file.path(root, 'pooled_base.rds'))
r <- readRDS(file.path(root, 'pooled_audit_base_repeat.rds'))
sa <- collect(a); sr <- collect(r)
stopifnot(identical(a$data_md5, r$data_md5))
repeat_diff <- cbind(abs(sa$q5 - sr$q5) / sqrt(sa$mcse_q05^2 + sr$mcse_q05^2),
                     abs(sa$q50 - sr$q50) / sqrt(sa$mcse_q50^2 + sr$mcse_q50^2),
                     abs(sa$q95 - sr$q95) / sqrt(sa$mcse_q95^2 + sr$mcse_q95^2))
cat('\nSame-baseline seed-101 repeat: data hash', a$data_md5,
    'draws identical', identical(a$fit[[1]]$samples$alpha, r$fit[[1]]$samples$alpha),
    'max quantile difference in combined MCSE', max(repeat_diff),
    'quantile checks beyond 4 combined MCSE', sum(repeat_diff > 4), '/', length(repeat_diff), '\n')
