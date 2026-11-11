args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_reference.R results-dir')
suppressPackageStartupMessages(library(posterior))
for (case in c('EXG', 'LNR')) for (seed in c(101, 202, 303)) {
  path <- file.path(args[[1]], paste(case, seed, sep = '_'))
  a <- readRDS(paste0(path, '_base_shared.rds'))
  b <- readRDS(paste0(path, '_archive_shared.rds'))
  cat('\n', case, seed, 'data hash equal', identical(a$data_md5, b$data_md5), '\n')
  for (arm in c('base', 'archive')) {
    x <- if (arm == 'base') a else b
    arr <- aperm(simplify2array(lapply(x$fit, function(ch) ch$samples$alpha[, 1, ])), c(2,3,1))
    dimnames(arr) <- list(NULL, NULL, x$fit[[1]]$par_names)
    stats <- summarise_draws(as_draws_array(arr), 'ess_bulk', 'ess_tail', 'rhat', 'mcse_mean',
      ~quantile2(.x, probs = c(.05, .5, .95)),
      ~setNames(vapply(c(.05,.5,.95), function(p) mcse_quantile(.x, probs = p), numeric(1)),
                c('mcse_q05','mcse_q50','mcse_q95')))
    cat(arm, 'elapsed', x$timing[['elapsed']], 'min bulk/draw', min(stats$ess_bulk) / prod(dim(arr)[1:2]),
        'min tail/draw', min(stats$ess_tail) / prod(dim(arr)[1:2]),
        'min bulk/s', min(stats$ess_bulk) / x$timing[['elapsed']],
        'min tail/s', min(stats$ess_tail) / x$timing[['elapsed']], '\n')
    print(stats, width = Inf)
  }
}
