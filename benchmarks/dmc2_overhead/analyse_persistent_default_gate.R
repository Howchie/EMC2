args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_persistent_default_gate.R results-dir')
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
root <- args[[1L]]
fields <- c('alpha', 'stage', 'subj_ll', 'idx')

for (case in c('EXG', 'LNR')) for (seed in c(101, 202, 303)) {
  prefix <- file.path(root, paste(case, seed, sep = '_'))
  baseline <- readRDS(paste0(prefix, '_fork_per_block_default.rds'))
  candidate <- readRDS(paste0(prefix, '_persistent_default.rds'))
  stopifnot(identical(baseline$data_md5, candidate$data_md5),
            identical(baseline$initial_seed, candidate$initial_seed),
            isFALSE(baseline$options$emc2.persistent_chains),
            isTRUE(candidate$options$emc2.persistent_chains))

  sample_equal <- vapply(seq_along(baseline$fit), function(i)
    all(vapply(fields, function(f)
      identical(baseline$fit[[i]]$samples[[f]], candidate$fit[[i]]$samples[[f]]),
      logical(1))), logical(1))
  rng_equal <- vapply(seq_along(baseline$fit), function(i)
    identical(baseline$fit[[i]]$rng, candidate$fit[[i]]$rng), logical(1))
  continuation <- setdiff(union(names(baseline$fit[[1L]]),
                                names(candidate$fit[[1L]])),
                          c('data', 'prior', 'model', 'samples', 'rng'))
  state_equal <- vapply(continuation, function(f)
    all(vapply(seq_along(baseline$fit), function(i)
      identical(baseline$fit[[i]][[f]], candidate$fit[[i]][[f]]), logical(1))),
    logical(1))
  counts_equal <- identical(EMC2::chain_n(baseline$fit),
                            EMC2::chain_n(candidate$fit))
  caller_rng_equal <- identical(baseline$final_seed, candidate$final_seed)
  timing <- c(fork_per_block = baseline$timing[['elapsed']],
              persistent = candidate$timing[['elapsed']])
  cat('\n', case, seed,
      'data', baseline$data_md5,
      'elapsed', paste(names(timing), signif(timing, 5), collapse = ' / '),
      'speedup', signif(timing[['fork_per_block']] / timing[['persistent']], 5),
      '\n  sample fields identical:', sample_equal,
      'chain RNG identical:', rng_equal,
      'caller RNG identical:', caller_rng_equal,
      'draw counts identical:', counts_equal,
      '\n  continuation fields differing:',
      paste(continuation[!state_equal], collapse = ', '), '\n')
  if (!all(sample_equal) || !all(rng_equal) || !caller_rng_equal ||
      !counts_equal || !all(state_equal))
    stop('persistent-chain parity gate failed for ', case, ' seed ', seed)
}
