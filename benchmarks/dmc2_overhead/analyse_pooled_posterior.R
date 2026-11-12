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
      paste('alpha', dimnames(ch$samples$alpha)[[1L]], s, sep = '__')), use.names = FALSE)
    mu <- t(ch$samples$theta_mu)
    colnames(mu) <- paste('mu', rownames(ch$samples$theta_mu), sep = '__')
    V <- ch$samples$theta_var
    sigma2 <- t(vapply(seq_len(dim(V)[3L]), function(i) diag(V[, , i]),
                       numeric(dim(V)[1L])))
    colnames(sigma2) <- paste('sigma2', rownames(V), sep = '__')
    log_sd <- 0.5 * log(sigma2)
    colnames(log_sd) <- paste('log_sd', rownames(V), sep = '__')
    pooled <- t(ch$samples$pooled)
    colnames(pooled) <- paste('pooled', rownames(ch$samples$pooled), sep = '__')
    cbind(a, mu, sigma2, log_sd, pooled)
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
hardfit_vars <- NULL
hardfit_ratios <- list(bulk = NULL, tail = NULL)
time_speedups <- numeric(length(seeds))
hardfit_pass <- matrix(FALSE, length(seeds), 2L,
                       dimnames = list(seeds, c('base', 'candidate')))
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
  hardfit_ess_vars <- grepl('^(mu|log_sd)__', sa$variable)
  if (is.null(hardfit_vars)) {
    hardfit_vars <- sa$variable[hardfit_ess_vars]
    hardfit_ratios <- lapply(c('bulk', 'tail'), function(x)
      matrix(NA_real_, length(seeds), length(hardfit_vars),
             dimnames = list(seeds, hardfit_vars)))
    names(hardfit_ratios) <- c('bulk', 'tail')
  }
  for (kind in c('bulk', 'tail')) {
    col <- paste0('ess_', kind)
    ratios[i, kind] <- min(sb[[col]][hardfit_ess_vars]) /
      min(sa[[col]][hardfit_ess_vars])
    hardfit_ratios[[kind]][i, ] <- sb[[col]][hardfit_ess_vars] /
      sa[[col]][hardfit_ess_vars]
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
  get_rhat <- function(fit, selection) unlist(EMC2::gd_summary(
    fit, selection = selection, stage = 'sample', stat = NULL, version = 'old'))
  ar <- lapply(list(a$fit, b$fit), get_rhat, selection = 'alpha')
  mr <- lapply(list(a$fit, b$fit), get_rhat, selection = 'mu')
  vr <- lapply(list(a$fit, b$fit), get_rhat, selection = 'sigma2')
  hess <- c(min(sa$ess_bulk[hardfit_ess_vars]), min(sb$ess_bulk[hardfit_ess_vars]))
  rule_max <- c(max(quantile(ar[[1L]], .99), max(mr[[1L]]), max(vr[[1L]])),
                max(quantile(ar[[2L]], .99), max(mr[[2L]]), max(vr[[2L]])))
  hardfit_pass[i, ] <- c(rule_max[1L] < 1.1 && hess[1L] >= 100,
                         rule_max[2L] < 1.1 && hess[2L] >= 100)
  cat('DMC2 hardfit-style R-hat q99(alpha), max(mu), max(sigma2), and min bulk ESS; base/candidate:',
      sprintf('%.3f/%.3f, %.3f/%.3f, %.3f/%.3f, %.1f/%.1f',
              quantile(ar[[1L]], .99), quantile(ar[[2L]], .99),
              max(mr[[1L]]), max(mr[[2L]]), max(vr[[1L]]), max(vr[[2L]]),
              hess[1L], hess[2L]), '\n')
  cat('DMC2 hardfit R-hat+ESS subset pass (max R-hat < 1.1; min group/log-SD bulk ESS >= 100), base/candidate:',
      rule_max[1L] < 1.1 && hess[1L] >= 100,
      rule_max[2L] < 1.1 && hess[2L] >= 100,
      '(omits likelihood-gap and subject-spread checks)\n')
  cat('quantile differences beyond 4 combined MCSE:', sum(diffs > 4), '/', length(diffs),
      'largest standardized difference', max(diffs), '\n')
  cat('largest quantile standardized differences by variable:\n')
  print(data.frame(variable = sa$variable, max_mcse = apply(diffs, 1L, max)), row.names = FALSE)
  cat('candidate diagnostics:\n'); print(sb, width = Inf)
}
cat('\nPer-seed minimum group-mean/log-SD ESS ratios (candidate/base):\n'); print(ratios)
cat('Partial DMC2 hardfit R-hat+ESS subset passes (omits likelihood-gap and subject-spread): base',
    sum(hardfit_pass[, 'base']), '/', nrow(hardfit_pass), '; candidate',
    sum(hardfit_pass[, 'candidate']), '/', nrow(hardfit_pass), '\n')
for (kind in colnames(ratios)) {
  z <- log(ratios[, kind]); ci <- exp(mean(z) + c(-1, 1) * qt(.975, df = length(z) - 1) * sd(z) / sqrt(length(z)))
  cat(kind, 'geometric mean ratio', exp(mean(z)), '95% t interval', ci, '\n')
  v <- as.vector(hardfit_ratios[[kind]])
  cat(kind, 'DMC2 group/log-SD ESS ratios per variable/run: median', median(v),
      '10th percentile', unname(quantile(v, .1)), 'minimum', min(v),
      'below 0.90', sum(v < .90), '/', length(v), '\n')
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
  hvars <- grepl('^(mu|log_sd)__', sa$variable)
  for (kind in c('bulk', 'tail')) {
    col <- paste0('ess_', kind)
    x <- sr[[col]][hvars] / sa[[col]][hvars]
    cat('Same-baseline repeat group-mean/log-SD', kind,
        'ESS ratio median / geometric mean / minimum:', median(x),
        exp(mean(log(x))), min(x), '\n')
  }
}
