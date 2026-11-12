args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_persistent_default_gate.R results-dir')
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
root <- args[[1L]]
fields <- c('alpha', 'stage', 'subj_ll', 'idx')
rows <- list()

time_metrics <- function(path) {
  lines <- readLines(path, warn = FALSE)
  value <- function(label) {
    hit <- lines[startsWith(trimws(lines), label)]
    if (length(hit) != 1L) return(NA_real_)
    as.numeric(gsub(',', '', sub('.*: *', '', hit)))
  }
  wall_text <- lines[startsWith(trimws(lines), 'Elapsed (wall clock) time')]
  wall <- NA_real_
  if (length(wall_text) == 1L) {
    stamp <- sub('^.*\\): *', '', wall_text)
    parts <- as.numeric(strsplit(stamp, ':', fixed = TRUE)[[1L]])
    wall <- if (length(parts) == 2L) 60 * parts[1L] + parts[2L]
      else if (length(parts) == 3L) 3600 * parts[1L] + 60 * parts[2L] + parts[3L]
      else NA_real_
  }
  c(user = value('User time (seconds):'),
    system = value('System time (seconds):'),
    wall = wall,
    max_rss_kb = value('Maximum resident set size (kbytes):'),
    minor_faults = value('Minor (reclaiming a frame) page faults:'))
}

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
  counts <- EMC2::chain_n(baseline$fit)
  counts_equal <- identical(counts, EMC2::chain_n(candidate$fit))
  caller_rng_equal <- identical(baseline$final_seed, candidate$final_seed)
  timing <- c(fork_per_block = baseline$timing[['elapsed']],
              persistent = candidate$timing[['elapsed']])
  resources <- rbind(
    fork_per_block = time_metrics(paste0(prefix, '_fork_per_block_default.time')),
    persistent = time_metrics(paste0(prefix, '_persistent_default.time')))
  cat('\n', case, seed,
      'data', baseline$data_md5,
      'fit elapsed', paste(names(timing), signif(timing, 5), collapse = ' / '),
      'fit speedup', signif(timing[['fork_per_block']] / timing[['persistent']], 5),
      'wall speedup',
      signif(resources['fork_per_block', 'wall'] / resources['persistent', 'wall'], 5),
      '\n  stage draws (preburn/burn/adapt/sample):',
      paste(counts[1L, c('preburn', 'burn', 'adapt', 'sample')], collapse = '/'),
      '\n  sample fields identical:', sample_equal,
      'chain RNG identical:', rng_equal,
      'caller RNG identical:', caller_rng_equal,
      'draw counts identical:', counts_equal,
      '\n  process resources (user/system/wall seconds, max RSS MiB, minor faults):\n',
      paste(capture.output(print(round(cbind(
        resources[, c('user', 'system', 'wall')],
        max_rss_mib = resources[, 'max_rss_kb'] / 1024,
        minor_faults = resources[, 'minor_faults']), 2))), collapse = '\n'),
      '\n  continuation fields differing:',
      paste(continuation[!state_equal], collapse = ', '), '\n')
  if (!all(sample_equal) || !all(rng_equal) || !caller_rng_equal ||
      !counts_equal || !all(state_equal))
    stop('persistent-chain parity gate failed for ', case, ' seed ', seed)
  rows[[paste(case, seed, sep = '_')]] <- data.frame(
    case = case, seed = seed,
    fit_speedup = timing[['fork_per_block']] / timing[['persistent']],
    wall_speedup = resources['fork_per_block', 'wall'] / resources['persistent', 'wall'],
    user_ratio = resources['persistent', 'user'] / resources['fork_per_block', 'user'],
    system_ratio = resources['persistent', 'system'] / resources['fork_per_block', 'system'],
    max_rss_ratio = resources['persistent', 'max_rss_kb'] / resources['fork_per_block', 'max_rss_kb'],
    minor_fault_ratio = resources['persistent', 'minor_faults'] / resources['fork_per_block', 'minor_faults'])
}

all_rows <- do.call(rbind, rows)
for (case in unique(all_rows$case)) {
  x <- all_rows[all_rows$case == case, , drop = FALSE]
  z <- log(x$fit_speedup)
  ci <- exp(mean(z) + c(-1, 1) * stats::qt(.975, df = length(z) - 1L) *
              stats::sd(z) / sqrt(length(z)))
  cat('\n', case, 'geometric fit speedup', signif(exp(mean(z)), 5),
      'paired 95% t interval', paste(signif(ci, 5), collapse = '–'),
      'geometric wall speedup', signif(exp(mean(log(x$wall_speedup))), 5),
      '\n  mean persistent/fork ratios user/system/max-RSS/minor-faults:',
      paste(signif(colMeans(x[, c('user_ratio', 'system_ratio',
                                  'max_rss_ratio', 'minor_fault_ratio')]), 5),
            collapse = '/'), '\n')
}
