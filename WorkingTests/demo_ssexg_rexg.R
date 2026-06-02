library(EMC2)
rm(list = ls())
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

N_TRIALS <- 1e4
RT_RESOLUTION <- 1e-3

# Collaborator-style 1-choice detection setup for SSEXG:
# force all trials to be "go" by setting lI = 2 and SSD = Inf.
lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)
SSD <- function(d) rep(Inf, nrow(d))

design_EXG <- design(
  factors = list(subjects = 1, S = 1),
  Rlevels = 1,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, exg_lb ~ 1),
  constants = c(exg_lb = -Inf),
  functions = list(lI = lIfun_go, SSD = SSD),
  model = SSEXG
)

pexg1 <- log(c(mu = .5, sigma = .1, tau = .2))

dat <- make_data(
  pexg1,
  design_EXG,
  n_trials = N_TRIALS,
  rt_resolution = RT_RESOLUTION
)


exg <- make_emc(dat, design_EXG, type = "single")
exg <- fit(exg)
cat("\nUncensored SSEXG fit:\n")
print(credint(exg, map = TRUE))
recovery(exg,true_pars = pexg1)


datC <- make_missing(
  dat,
  LC = .47,
  UC = .95,
  verbose = TRUE,
  rt_resolution = RT_RESOLUTION
)

plot_density(datC)


exgC <- make_emc(datC, design_EXG, type = "single")
exgC <- fit(exgC)
cat("\nCensored SSEXG fit:\n")
print(credint(exgC, map = TRUE))
recovery(exgC,true_pars = pexg1)

## New - Racing Ex-Gaussian model to keep it simple
# Equivalent 1-choice ex-Gaussian race model.
design_REXG <- design(
  factors = list(subjects = 1, S = 1),
  Rlevels = 1,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1),
  model = REXG
)

dat_race <- make_data(
  pexg1,
  design_REXG,
  n_trials = N_TRIALS,
  rt_resolution = RT_RESOLUTION
)

# Fit against new model as generator

rexg <- make_emc(dat_race, design_REXG, type = "single")
rexg <- fit(rexg)
cat("\nUncensored REXG fit:\n")
print(credint(rexg, map = TRUE))
recovery(rexg,true_pars = pexg1)

# Fit against the SSEXG-generated data to show equivalence

rexg2 <- make_emc(dat, design_REXG, type = "single")
rexg2 <- fit(rexg2)
cat("\nUncensored REXG fit:\n")
print(credint(rexg2, map = TRUE))
recovery(rexg2,true_pars = pexg1)

datC_race <- make_missing(
  dat_race,
  LC = .47,
  UC = .95,
  verbose = TRUE,
  rt_resolution = RT_RESOLUTION
)

plot_density(datC_race)

# Fit against new model as generator with UC/LC
rexgC <- make_emc(datC_race, design_REXG, type = "single")
rexgC <- fit(rexgC)
cat("\nCensored REXG fit:\n")
print(credint(rexgC, map = TRUE))
recovery(rexgC,true_pars = pexg1)

# Fit against the SSEXG-generated data with UC/LC to show equivalence

rexgC2 <- make_emc(datC, design_REXG, type = "single")
rexgC2 <- fit(rexgC2)
cat("\nCensored REXG fit:\n")
print(credint(rexgC2, map = TRUE))
recovery(rexgC2,true_pars = pexg1)
