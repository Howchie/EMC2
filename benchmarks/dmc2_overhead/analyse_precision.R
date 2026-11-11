args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_precision.R results-dir')
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
for (case in c('EXG', 'LNR')) for (seed in c(101, 202, 303)) {
  path <- file.path(args[[1]], paste(case, seed, sep = '_'))
  for (arm in c('base', 'new')) {
    x <- readRDS(paste0(path, '_', arm, '_precision.rds'))
    ess <- EMC2::ess_summary(x$fit, selection = 'alpha', stage = 'sample', stat_only = TRUE, n_cores = 1)
    n <- sum(x$fit[[1]]$samples$stage == 'sample', na.rm = TRUE)
    cat(case, seed, arm, 'elapsed', x$timing[['elapsed']], 'sample draws/chain', n,
        'min ESS', min(ess), 'target reached', min(ess) > 800, '\n')
  }
}
