#remotes::install_github("https://github.com/Howchie/EMC2/",ref="censoring_truncation_ah")
# Notes:
# Go/NoGo designs are detected via the use of a "nogo" level in the response.
# Key likelihood function is in c_log_likelihood_race in particle_ll.cpp (ignore any R likelihood). THere's a gng branch for missing RTs and otherwise finite RTs use the normal branch (treats go as having beaten nogo)
# model_LBA.R, design.R (has branches when any Rlevel=="nogo"), make_missing sets nogo wins to rt=Inf, R=NA and design sets winner=lR: nogo for those combinations. 
rm(list=ls())
library(EMC2)
set.seed(123)

## Test a small easy model
# NB setting a different lMFALSE for S failed to converge
# With separate drifts it fails to converge at normal levels
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 0.5, -.5),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  model=LBA,UC=3,
  formula=list(v~match:S,A~1,B~1,t0~1,sv~1),
  constants = c(A=log(0.5),sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1,1, 1.5, log(1.2), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)
# All responses that would have been nogo are coded as NA with rt=Inf
unique(dat$R)
# Response probabilities as expected
tapply(is.na(dat$R),dat$S,mean)
print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designLBA,type="single")
system.time({emc1 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3
))})
print(recovery(emcc,p_vector,selection="alpha"))

# Test more realistic structure
# Go/NoGo designs are detected via the use of a "nogo" level in the response.
n_subj=15
designLBA <- design(
  factors=list(subjects=1:n_subj,S=c("go","nogo","nogo","nogo"), L=c("Low","Med","High")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,UC=1.5,
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 0.5, -0.5),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  formula=list(v~S:match:L, B~lR, t0~1, A~1,sv~1),
  constants = c(sv=log(1),A=log(.5))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1,1.5,1.8,1,1.5,1,1, log(1.15),log(1), log(0.15))

# IMPORTANT:
# profile_plot assumes one shared parameter vector for the full dataset.
# So to validate the likelihood shape, use fixed-effect simulated data.
dat_fixed <- make_data(p_vector,designLBA,n_trials=200)

# Hierarchical simulation for recovery should use an explicit covariance matrix.
# The default make_random_effects() scales variances by abs(mean), which can
# create very large variances and exact zero variances (if a mean is 0).
re_sd <- c(rep(0.15, 7), 0.10, 0.10, 0.08)
re_cov <- diag(re_sd^2)
tmp <- make_random_effects(designLBA,p_vector,n_subj=n_subj,covariances = re_cov)
dat <- make_data(tmp,designLBA,n_trials=200)

# All responses that would have been nogo are coded as NA with rt=Inf
unique(dat$R)

# Response probabilities as expected
tapply(is.na(dat$R),interaction(dat$S,dat$L),mean)

# Likelihood profile sanity check (fixed-effect data only)
print(profile_plot(dat_fixed,designLBA,p_vector,n_cores=1,layout=c(2,2)))

# Fitting
emc <- make_emc(dat,designLBA,type="standard", chains=3)
system.time({emc2 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), fileName = "gng_fit.RData")})
save.image("gng_test.RData")
# Both good
print(recovery(emcc,p_vector,selection="mu"))

# Recovers data trends
pred=predict(emcc)
tapply(is.na(dat$R),dat$S,mean)
tapply(is.na(pred$R),pred$S,mean)
tapply(dat$rt,dat$S,function(x) {mean(x[is.finite(x)])})
tapply(pred$rt,pred$S,function(x) {mean(x[is.finite(x)])})
