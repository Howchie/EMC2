
rm(list=ls())
library(EMC2)
set.seed(123)

## Test 1 - Single error drift, match x S, threshold by S. -- CONVERGES
# First use separate drifts (not base + diff)

designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  model=LBA,UC=3,
  formula=list(v~0+mismatch + match:S,A~1,B~0+lR,t0~1,sv~1),
  constants = c(sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.5, 2.2, 2.5, log(0.4), log(.75), log(0.6), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
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
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc1,p_vector,selection="alpha"))

## Test 2, identical to above but different mismatch accumulators -- DOESN'T RELIABLY CONVERGE
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  model=LBA,UC=3,
  formula=list(v~0+mismatch:S + match:S,A~1,B~0+lR,t0~1,sv~1),
  constants = c(A=log(0.4),sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.6,.4, 2.2, 2.5, log(.75), log(.6), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designLBA,type="single")
system.time({emc2 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc2,p_vector,selection="alpha"))

## Test 3, identical to above but single threshold  -- CONVERGES
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  model=LBA,UC=3,
  formula=list(v~0+mismatch:S + match:S,A~1,B~1,t0~1,sv~1),
  constants = c(A=log(0.4),sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.5,.35, 2.2, 2.5, log(.3), log(.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designLBA,type="single")
system.time({emc3 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc3,p_vector,selection="alpha"))

## Test 4, individual match/mismatch x S, B mean vs difference form -- DOESN'T CONVERGE
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  model=LBA,UC=3,
  formula=list(v~0+mismatch:S + match:S,A~1,B~lR,t0~1,sv~1),
  constants = c(A=log(0.4),sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.5,.35, 2.2, 2.5, log(.4), log(0.9), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designLBA,type="single")
system.time({emc4 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc4,p_vector,selection="alpha"))

## Test 5 - Threshold + drifts difference form -- DOESN'T CONVERGE

designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 0.5, -.5)
  ),
  model=LBA,UC=3,
  formula=list(v~match:S,A~1,B~lR,t0~1,sv~1),
  constants = c(A=log(0.4),sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1,2.2, 2.5, log(.5), log(0.85), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)
tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designLBA,type="single")
system.time({emc5 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc5,p_vector,selection="alpha"))

# Test more realistic structure with drift (correct) by load + threshold differs between the accumulators. Also unfix A.
# hierarchical as well
n_subj=15

designLBA <- design(
  factors=list(subjects=1:n_subj,S=c("go","nogo","nogo","nogo"), L=c("Low","Med","High")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,UC=3,
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  formula=list(v~0+mismatch + S:match:L, B~0+lR, t0~1, A~1,sv~1),
  constants = c(sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.2,2.1,2.3,1.65,1.95,1.1,1.65, log(.4),log(.25), log(0.15), log(0.5))

# IMPORTANT:
# profile_plot assumes one shared parameter vector for the full dataset.
# So to validate the likelihood shape, use fixed-effect simulated data.
dat_fixed <- make_data(p_vector,designLBA,n_trials=200)

# Hierarchical simulation for recovery should use an explicit covariance matrix.
# The default make_random_effects() scales variances by abs(mean), which can
# create very large variances and exact zero variances (if a mean is 0).
# I modified it to baseline at natural mean on the transformed scale zero (e.g. for the default log(1)=0, we'd get a variance of variance_prop * exp(0))

tmp <- make_random_effects(designLBA,p_vector,n_subj=n_subj, variance_proportion = .1)
dat <- make_data(tmp,designLBA,n_trials=200)

# Response probabilities as expected
tapply(is.na(dat$R),interaction(dat$S,dat$L),mean)
tapply(dat$rt,interaction(dat$S,dat$L),function(x){mean(x[is.finite(x)],na.rm=TRUE)})
# Likelihood profile sanity check (fixed-effect data only)
#print(profile_plot(dat_fixed,designLBA,p_vector,n_cores=1,layout=c(2,2)))

# Fitting
emc <- make_emc(dat,designLBA,type="standard", chains=3)
system.time({emc_big <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries= 30, fileName = "gng_fit.RData")})
# Both good
print(recovery(emc_big,p_vector,selection="mu"))
save.image("gng_test.RData")


# Add threshold load effect as well (this didn't converge by max_tries = 20 but the output of recovery shows it was getting there)
designLBA <- design(
  factors=list(subjects=1:n_subj,S=c("go","nogo","nogo","nogo"), L=c("Low","Med","High")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,UC=3,
  functions = list(
    match = function(d) ifelse(d$lM==TRUE, 0.5, -0.5),
    mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
  ),
  formula=list(v~0+mismatch + S:match:L, B~0+lR:L, t0~1, A~1,sv~1),
  constants = c(sv=log(1))
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.2,2.1,2.3,1.65,1.95,1.1,1.65, log(.3),log(.2), log(.45),log(.35), log(.7),log(.6),log(0.15), log(0.5))


tmp2 <- make_random_effects(designLBA,p_vector,n_subj=n_subj,variance_proportion=.1)
dat <- make_data(tmp,designLBA,n_trials=200)

tapply(is.na(dat$R),interaction(dat$S,dat$L),mean)
tapply(dat$rt,interaction(dat$S,dat$L),function(x){mean(x[is.finite(x)],na.rm=TRUE)})
# Fitting
emc <- make_emc(dat,designLBA,type="standard", chains=3)
system.time({emc_full2 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries = 30, fileName = "gng_fit.RData")})

print(recovery(emc_full2,p_vector,selection="mu"))
save.image("WorkingTests/gng_unidentifiabilities.RData")


# Recovers data trends
pred=predict(emc_full2)
tapply(is.na(dat$R),dat$S,mean)
tapply(is.na(pred$R),pred$S,mean)
tapply(dat$rt,dat$S,function(x) {mean(x[is.finite(x)])})
tapply(pred$rt,pred$S,function(x) {mean(x[is.finite(x)])})
