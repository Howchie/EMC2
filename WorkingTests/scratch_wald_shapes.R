rm(list = ls())
library(EMC2)
designRDM <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
  model = RDM(),
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
  constants = c(s = log(1),A=log(0))
)
p_vector <- sampled_pars(designRDM, doMap = FALSE)
p_vector["t0"] <- log(0.15)
p_vector["B"] <- log(1)
p_vector["v"] <- log(1)
p_vector["v_lMTRUE"] <- log(3)
set.seed(42)
dat <- make_data(p_vector, designRDM, n_trials = 100000)
dat$Correct <- as.numeric(dat$S==dat$R)
cat("Wald Means:\n")
print(tapply(dat$rt, dat$Correct, mean))
cat("Wald Quantiles correct:\n")
print(quantile(dat$rt[dat$Correct==1], probs=c(.1,.3,.5,.7,.9)))
cat("Wald Quantiles error:\n")
print(quantile(dat$rt[dat$Correct==0], probs=c(.1,.3,.5,.7,.9)))
