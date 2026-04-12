#!/usr/bin/env Rscript
.libPaths(c('/data/ubuntu-relocated/R/x86_64-pc-linux-gnu-library/4.5', .libPaths()))
suppressPackageStartupMessages(library(EMC2))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260410)

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designLBA <- design(
  factors = list(subjects = 1, S = c('left','right')),
  Rlevels = c('left','right'),
  matchfun = matchfun,
  model = LBA(posdrift = TRUE),
  formula = list(B~1, v~0+lM, A~1, t0~1, sv~0+lM),
  constants = c(A = log(0.4))
)

p <- sampled_pars(designLBA, doMap = FALSE)
p[['B']] <- log(1.2); p[['A']] <- log(0.4); p[['t0']] <- log(0.15)
p[['v_lMFALSE']] <- 0.4; p[['v_lMTRUE']] <- 1.2
p[['sv_lMFALSE']] <- log(1.2); p[['sv_lMTRUE']] <- log(1)

# Base data once (no censoring in simulated RT/R)
set.seed(99)
dat_base <- make_data(p, designLBA, n_trials = 4000)

# Two copies differing only in UC metadata
A <- dat_base
B <- dat_base
A$UC <- Inf
B$UC <- 100
# Keep the rest explicitly aligned
A$LC <- 0; A$LT <- 0; A$UT <- Inf
B$LC <- 0; B$LT <- 0; B$UT <- Inf

run_fit <- function(dat, label) {
  emc <- make_emc(dat, designLBA, type = 'single')
  set.seed(777)
  emc <- fit(
    emc,
    fileName = paste0('tmp_', label, '.RData'),
    stop_criteria = list(
      sample = list(iter = 250, max_sample_iter = 350, max_gd = 1.20,
                    max_flat_loc = 0.8, flat_selection = c('alpha','subj_ll')),
      cores_per_chain = 1,
      cores_for_chains = 1
    ),
    max_tries = 1
  )
  alpha <- get_pars(emc, selection = 'alpha', stage = 'sample', return_mcmc = FALSE, merge_chains = TRUE, remove_constants = FALSE)
  ll <- get_pars(emc, selection = 'LL', stage = 'sample', return_mcmc = FALSE, merge_chains = TRUE)
  list(emc = emc, alpha = alpha, ll = ll)
}

cat('Running fit A (UC=Inf) ...\n')
fitA <- run_fit(A, 'uc_inf')
cat('Running fit B (UC=100) ...\n')
fitB <- run_fit(B, 'uc_100')

# Summaries
ma <- apply(fitA$alpha, 1:2, mean)
mb <- apply(fitB$alpha, 1:2, mean)

cat('\n=== Posterior mean alpha diff (A-B) ===\n')
print(round(ma - mb, 6))

cat('\n=== LL summary ===\n')
lla <- fitA$ll[[1]][[1]]
llb <- fitB$ll[[1]][[1]]
cat('A mean LL:', mean(lla), ' max LL:', max(lla), '\n')
cat('B mean LL:', mean(llb), ' max LL:', max(llb), '\n')
cat('Delta mean:', mean(lla)-mean(llb), ' Delta max:', max(lla)-max(llb), '\n')
