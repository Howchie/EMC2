# Compare paired fixed-transition runs produced by run_paper.sh.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_paper.R results-dir')
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
for (case in c('EXG', 'LNR')) for (seed in c(101, 202, 303)) {
  prefix <- file.path(root, paste(case, seed, sep = '_'))
  a <- readRDS(paste0(prefix, '_base.rds'))
  b <- readRDS(paste0(prefix, '_new.rds'))
  stopifnot(identical(a$data_md5, b$data_md5), identical(a$initial_seed, b$initial_seed))
  sample_fields <- c('alpha', 'stage', 'subj_ll', 'idx')
  chain_equal <- vapply(seq_along(a$fit), function(j) {
    all(vapply(sample_fields, function(f) identical(a$fit[[j]]$samples[[f]], b$fit[[j]]$samples[[f]]), logical(1)))
  }, logical(1))
  rng_equal <- vapply(seq_along(a$fit), function(j) identical(a$fit[[j]]$rng, b$fit[[j]]$rng), logical(1))
  continuation <- setdiff(union(names(a$fit[[1]]), names(b$fit[[1]])),
                          c('data', 'prior', 'model'))
  same_fields <- vapply(continuation, function(f) all(vapply(seq_along(a$fit),
    function(j) identical(a$fit[[j]][[f]], b$fit[[j]][[f]]), logical(1))), logical(1))
  arr <- simplify2array(lapply(b$fit, function(ch) ch$samples$alpha[, 1, ]))
  arr <- aperm(arr, c(2, 3, 1))
  dimnames(arr) <- list(NULL, NULL, b$fit[[1]]$par_names)
  stats <- summarise_draws(as_draws_array(arr), 'mean', 'mcse_mean', 'ess_bulk', 'ess_tail',
                           ~quantile2(.x, probs = c(.05, .5, .95)),
                           ~setNames(vapply(c(.05, .5, .95), function(prob) mcse_quantile(.x, probs = prob), numeric(1)),
                                     c('mcse_q05', 'mcse_q50', 'mcse_q95')))
  cat('\n', case, seed, 'base elapsed', a$timing[['elapsed']], 'new elapsed', b$timing[['elapsed']],
      'data md5', a$data_md5, '\ndraws:', chain_equal, 'RNG:', rng_equal,
      'caller RNG:', identical(a$final_seed, b$final_seed),
      '\ncontinuation fields differing:', paste(continuation[!same_fields], collapse = ', '), '\n')
  print(stats, width = Inf)
  cat('minimum bulk ESS per draw', min(stats$ess_bulk) / (dim(arr)[1] * dim(arr)[2]),
      'minimum tail ESS per draw', min(stats$ess_tail) / (dim(arr)[1] * dim(arr)[2]),
      'minimum bulk ESS/s', min(stats$ess_bulk) / b$timing[['elapsed']],
      'minimum tail ESS/s', min(stats$ess_tail) / b$timing[['elapsed']], '\n')
}
