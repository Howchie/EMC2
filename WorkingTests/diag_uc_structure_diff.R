#!/usr/bin/env Rscript

.libPaths(c('/data/ubuntu-relocated/R/x86_64-pc-linux-gnu-library/4.5', .libPaths()))
suppressPackageStartupMessages(library(EMC2))

set.seed(123)

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

build_design <- function(posdrift = TRUE) {
  design(
    factors = list(subjects = 1, S = c('left', 'right')),
    Rlevels = c('left', 'right'),
    matchfun = matchfun,
    model = LBA(posdrift = posdrift),
    formula = list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ 0 + lM),
    constants = c(A = log(0.4))
  )
}

build_true_pars <- function(designLBA) {
  p <- sampled_pars(designLBA, doMap = FALSE)
  p[['B']] <- log(1.2)
  p[['A']] <- log(0.4)
  p[['t0']] <- log(0.15)
  p[['v_lMFALSE']] <- 0.4
  p[['v_lMTRUE']] <- 1.2
  p[['sv_lMFALSE']] <- log(1.2)
  p[['sv_lMTRUE']] <- log(1)
  p
}

summarise_dadm <- function(dadm, tag) {
  cat('\n===== ', tag, ' =====\n', sep = '')
  cat('nrow:', nrow(dadm), ' ncol:', ncol(dadm), '\n')
  cat('colnames:', paste(colnames(dadm), collapse = ', '), '\n')

  for (nm in c('LT', 'LC', 'UC', 'UT')) {
    if (nm %in% names(dadm)) {
      x <- dadm[[nm]]
      cat(sprintf('%s: class=%s unique_n=%d na=%d inf=%d min=%s max=%s\n',
                  nm,
                  paste(class(x), collapse='|'),
                  length(unique(x)),
                  sum(is.na(x)),
                  sum(is.infinite(x)),
                  suppressWarnings(as.character(min(x, na.rm = TRUE))),
                  suppressWarnings(as.character(max(x, na.rm = TRUE)))))
    } else {
      cat(nm, ': <missing column>\n', sep = '')
    }
  }

  key_attrs <- c('expand', 'unique_nort', 'expand_nort', 'unique_nortR', 'expand_nortR',
                 'emc2_ll_cache_version', 'emc2_all_finite_trials',
                 'finite_rt_unique_trial_indices', 'other_unique_trial_indices',
                 'RACE_nacc_by_row', 'RACE_mask')
  for (nm in key_attrs) {
    v <- attr(dadm, nm)
    if (is.null(v)) {
      cat('attr(', nm, '): <NULL>\n', sep = '')
    } else {
      cat('attr(', nm, '): len=', length(v), ' class=', paste(class(v), collapse='|'), sep = '')
      if (is.atomic(v) && length(v) <= 8) {
        cat(' values=', paste(v, collapse=','))
      }
      cat('\n')
    }
  }

  ds <- attr(dadm, 'designs')
  if (is.null(ds)) {
    cat('attr(designs): <NULL>\n')
  } else {
    cat('attr(designs): names=', paste(names(ds), collapse=', '), '\n')
    for (nm in names(ds)) {
      dm <- ds[[nm]]
      cat(sprintf('  design[%s]: dim=%s expand_len=%d\n',
                  nm,
                  paste(dim(dm), collapse='x'),
                  length(attr(dm, 'expand'))))
    }
  }
}

compare_vectors <- function(a, b, name) {
  if (is.null(a) && is.null(b)) {
    cat(name, ': both NULL\n')
    return(invisible(NULL))
  }
  if (is.null(a) || is.null(b)) {
    cat(name, ': NULL mismatch\n')
    return(invisible(NULL))
  }
  if (length(a) != length(b)) {
    cat(name, ': length mismatch ', length(a), ' vs ', length(b), '\n', sep = '')
    return(invisible(NULL))
  }
  same <- isTRUE(all.equal(a, b, check.attributes = FALSE))
  cat(name, ': ', if (same) 'identical' else 'DIFFERENT', '\n', sep = '')
  if (!same && is.atomic(a) && is.atomic(b)) {
    idx <- which(!(a == b) & !(is.na(a) & is.na(b)))
    idx <- idx[seq_len(min(length(idx), 10))]
    if (length(idx) > 0) {
      cat('  first mismatches idx=', paste(idx, collapse=','), '\n', sep='')
      cat('  A=', paste(a[idx], collapse=','), '\n', sep='')
      cat('  B=', paste(b[idx], collapse=','), '\n', sep='')
    }
  }
}

compare_dadm <- function(a, b) {
  cat('\n===== A vs B structural diff =====\n')
  cat('same colnames:', identical(colnames(a), colnames(b)), '\n')
  cat('same nrow:', nrow(a) == nrow(b), '\n')

  common_cols <- intersect(names(a), names(b))
  for (nm in c('rt', 'R', 'lR', 'LT', 'LC', 'UC', 'UT')) {
    if (nm %in% common_cols) {
      compare_vectors(a[[nm]], b[[nm]], paste0('col:', nm))
    }
  }

  for (nm in c('expand', 'unique_nort', 'expand_nort', 'unique_nortR', 'expand_nortR',
               'emc2_ll_cache_version', 'emc2_all_finite_trials',
               'finite_rt_unique_trial_indices', 'other_unique_trial_indices')) {
    compare_vectors(attr(a, nm), attr(b, nm), paste0('attr:', nm))
  }

  dsa <- attr(a, 'designs')
  dsb <- attr(b, 'designs')
  if (!is.null(dsa) && !is.null(dsb)) {
    cat('design names identical:', identical(names(dsa), names(dsb)), '\n')
    for (nm in intersect(names(dsa), names(dsb))) {
      dm_a <- dsa[[nm]]
      dm_b <- dsb[[nm]]
      cat(sprintf('design[%s] same dim: %s\n', nm, identical(dim(dm_a), dim(dm_b))))
      compare_vectors(attr(dm_a, 'expand'), attr(dm_b, 'expand'), paste0('design:', nm, ':expand'))
      if (identical(dim(dm_a), dim(dm_b))) {
        same_dm <- isTRUE(all.equal(dm_a, dm_b, check.attributes = FALSE))
        cat(sprintf('design[%s] values: %s\n', nm, if (same_dm) 'identical' else 'DIFFERENT'))
      }
    }
  }
}

designLBA <- build_design()
p <- build_true_pars(designLBA)

# A: UC unspecified (TC omitted)
set.seed(9991)
dat_A <- make_data(p, designLBA, n_trials = 10000)
# B: UC explicitly finite but practically non-censoring
set.seed(9991)
dat_B <- make_data(p, designLBA, n_trials = 10000, TC = list(UC = 100))

emc_A <- make_emc(dat_A, designLBA, type = 'single')
emc_B <- make_emc(dat_B, designLBA, type = 'single')

dadm_A <- emc_A[[1]]$data[[1]]
dadm_B <- emc_B[[1]]$data[[1]]

summarise_dadm(dadm_A, 'A: UC unspecified')
summarise_dadm(dadm_B, 'B: UC = 100')
compare_dadm(dadm_A, dadm_B)

cat('\nDone.\n')
