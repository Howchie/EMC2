rm(list = ls())
library(EMC2)
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
designROU <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = ROU(parameterization="equilibrium"),
  formula = list(theta ~ lM, B ~ 1, A ~ 1, t0 ~ 1, chi ~ 1, tk ~ 1),
  constants = c(B = log(1),A=log(0))
)
p_vector <- sampled_pars(designROU, doMap = FALSE)
p_vector["t0"] <- log(0.15)
p_vector["chi"] <- log(1)
p_vector["theta"] <- log(.5)
p_vector["theta_lMTRUE"] <- log(2.5)
p_vector["tk"] <- log(.5)

set.seed(42)
dat <- make_data(p_vector, designROU, n_trials = 1000)
dat$Correct <- as.numeric(dat$S==dat$R)

designRDM <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RDMSWTN(),
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
  constants = c(s = log(1),A=log(0))
)

emc2 <- make_emc(dat, designRDM, type = "single", rt_resolution = 1/60)
fit_res2 <- fit(
  emc2,
  cores_for_chains = 3,
  cores_per_chain = 2,
  stop_criteria = list(
    sample = list(
      iter = 200,
      max_gd = 1.10,
      max_flat_loc = 2,
      flat_selection = c("alpha", "subj_ll"),
      flat_p1 = 1/3,
      flat_p2 = 1/3,
      max_sample_iter = 500
    )
  ),
  max_tries = 1
)

print(summary(fit_res2))
