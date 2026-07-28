## Script to simulate and test recovery of Go/No-go Race Models
# Using the RDM for now as it is less unidentified than LBA (it seems to be the combined ridge of A + B + t0 that hurts most)
# Also exploring the "marginalize" option which integrates across the uncertainty of t0.
rm(list=ls())
library(EMC2)
set.seed(123)

## Test 1 - Shared base drift, match x S, threshold by S, single t0. 
# Easy converge, marginalise gets 'closer' to t0 truth but is much slower (2.4 mins vs 21 seconds)
marginalize="t0" # NULL or "t0"
designRDM <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)),
  model=RDM,UC=3,
  formula=list(v~match:S,A~1,B~lR,t0~1,s~lM),
  constants = c(s=log(1),A=log(0)),
  marginalise = marginalize
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(log(2.5), log(2.5), log(3), log(1.4), log(0.85), log(0.2), log(0.8))
prior_pvec = sampled_pars(designRDM);
prior_pvec[]=c(log(2),log(2),log(2),log(1),log(0.8),log(.3),log(1))
prior_svec = sampled_pars(designRDM)+1;
prior_svec["t0"]=.5
dat <- make_data(p_vector,designRDM, n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designRDM,type="single",
                prior = prior(designRDM,theta_mu_mean=prior_pvec,
                               theta_mu_sd = prior_svec,type="single"))
system.time({emc1 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
  )})
print(recovery(emc1,p_vector,selection="alpha"))

## Test 2, identical to above but different base drifts as well
# both converge, marginalise jumped around a bit and took more iterations to get a worse recovery
marginalize=NULL # or "t0"
designRDM <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)),
  model=RDM,UC=3,
  formula=list(v~S+match:S,A~1,B~lR,t0~1,s~lM),
  constants = c(s=log(1),A=log(0)),
  marginalise = marginalize
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(log(2.5), log(.7), log(2.5), log(3), log(1.4), log(0.85), log(0.2), log(0.8))
prior_pvec = sampled_pars(designRDM);
prior_pvec[]=c(log(2),log(2),log(2),log(2),log(1),log(0.8),log(.3),log(1))
prior_svec = sampled_pars(designRDM)+1;
prior_svec["t0"]=.5
dat <- make_data(p_vector,designRDM, n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designRDM,type="single",
                prior = prior(designRDM,theta_mu_mean=prior_pvec,
                              theta_mu_sd = prior_svec,type="single"))
system.time({emc1 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
)})
print(recovery(emc1,p_vector,selection="alpha"))

## Test 3, LBA version of the above model (all the same except with added drift)
# both converge fine
marginalize= NULL # NULL or "t0"
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)),
  model=LBA,UC=3,
  formula=list(v~S+match:S,A~1,B~lR,t0~1,sv~lM),
  constants = c(s=log(1)),
  marginalise = marginalize
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(2, -.5, 2.5, 3, log(1.4), log(0.85), log(0.2), log(0.8))
prior_pvec = sampled_pars(designRDM);
prior_pvec[]=c(2,0,2,2,log(1),log(0.8),log(.3),log(1))
prior_svec = sampled_pars(designRDM)+1;
prior_svec["t0"]=.5
dat <- make_data(p_vector,designRDM, n_trials=10000)

tapply(is.na(dat$R),dat$S,mean)
tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
#print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
emc <- make_emc(dat,designRDM,type="single",
                prior = prior(designRDM,theta_mu_mean=prior_pvec,
                              theta_mu_sd = prior_svec,type="single"))
system.time({emc1 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
)})
print(recovery(emc1,p_vector,selection="alpha"))
