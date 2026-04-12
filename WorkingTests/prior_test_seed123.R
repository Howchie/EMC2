## Test whether centering v_lMTRUE prior on 2 fixes seed=123 recovery failure
## Also notes: v is identity-scale (N(0,1) = drift in [-3,3], fine)
##             t0 is log-scale, true=-1.90 is 1.9sd from default prior mean

rm(list = ls())
library(EMC2)

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
designLBA <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = LBA(posdrift = TRUE),
  formula = list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM),
  constants = c(A = log(0.4))
)

p_true <- c(B=log(1.2), v_lMFALSE=0.4, v_lMTRUE=1.2, t0=log(0.15),
            sv_lMFALSE=log(1.2), sv_lMTRUE=log(1))
free_pars <- p_true  # all are free (A is constant)

# Generate data with seed=123 (the failing seed)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
pv <- sampled_pars(designLBA, doMap=FALSE)
pv[names(p_true)] <- p_true
dat <- make_data(pv, designLBA, n_trials=10000,
  TC=list(UC=Inf, UT=Inf, LC=0, LT=0,
          LCresponse=FALSE, UCresponse=FALSE,
          LCdirection=TRUE, UCdirection=TRUE))

cat("Data generated with seed=123\n")
cat(sprintf("RT range: [%.3f, %.3f], mean=%.3f\n",
            min(dat$rt[is.finite(dat$rt)]),
            max(dat$rt[is.finite(dat$rt)]),
            mean(dat$rt[is.finite(dat$rt)])))

# Build prior: v_lMTRUE centred on 2, everything else default N(0,1)
# Default: theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(1, n_pars)
# Only change: v_lMTRUE mean from 0 -> 2
pmean_mod <- c(B=0, v_lMFALSE=0, v_lMTRUE=2, t0=0, sv_lMFALSE=0, sv_lMTRUE=0)
psd_mod   <- c(B=1, v_lMFALSE=1, v_lMTRUE=1, t0=1, sv_lMFALSE=1, sv_lMTRUE=1)

prior_mod <- prior(designLBA, type="single",
                   pmean = pmean_mod,
                   psd   = psd_mod)
cat("\nModified prior summary:\n")
print(summary(prior_mod))

emc <- make_emc(dat, designLBA, type="single", rt_resolution=NULL,
                prior = prior_mod)

emc <- fit(emc,
  fileName = "samples_ss_seed123_prior_vtrue2.RData",
  stop_criteria = list(
    sample = list(
      iter = 1000,
      max_gd = 1.10,
      max_flat_loc = 0.5,
      flat_selection = c("alpha", "subj_ll"),
      flat_p1 = 1/3,
      flat_p2 = 1/3,
      max_sample_iter = 5000
    )
  ),
  cores_per_chain = 3, cores_for_chains = 3, max_tries = 30
)

cat("\n=== RECOVERY: seed=123, v_lMTRUE prior centred on 2 ===\n")
rec <- recovery(emc, true_pars = free_pars)
print(rec)
