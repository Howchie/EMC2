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
dat <- make_data(p_vector, designROU, n_trials = 100000)
dat$Correct <- as.numeric(dat$S==dat$R)
cat("Means:\n")
print(tapply(dat$rt, dat$Correct, mean))
cat("Quantiles correct:\n")
print(quantile(dat$rt[dat$Correct==1], probs=c(.1,.3,.5,.7,.9)))
cat("Quantiles error:\n")
print(quantile(dat$rt[dat$Correct==0], probs=c(.1,.3,.5,.7,.9)))
