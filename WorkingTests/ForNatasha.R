## Script to simulate and test recovery of Go/No-go Race Models
# Using the RDM for now as it is less unidentified than LBA
library(EMC2)
library(dplyr)
rm(list=ls())
set.seed(123)
## Note - the issue when we chatted was setting TC=list(UC=3) globally.
# It seems like for the current build there's  acouple of unique quirks.
# For "nogo" outcomes, UC needs to be set to 0, so that it then integrates across the entire time range
# (previously it was trying to integrate from [UC, Inf])
# The complication there means we should *not* use any UC for the other trials because it won't do anything anyway
# i.e. for now, every "nogo" is considered a "nogo accumulator won".
# In future I assume this will change
# NB for using real data, drop the factors argument and just add something like data = data
# ensuring your data has columns with *factors* S, subjects, R (levels includes "nogo" but responses coded as NA)
designRDM <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)), # match here is identical to the ADMat lM from the help file, either way of specifying it works
  model=RDM,TC=list(UC=999), # NB: For *simulation* only, we need a finite UC specified here
  formula=list(v~match:S,A~1,B~lR,t0~1,s~lM),
  constants = c(s=log(1),A=log(0)) 
  # I set no start-point variability for RDM because it often doesn't need it, 
  # but if you want it just remove the A=log(0) from constants and it will estimate
  # a value instead
)

## This code is only used for the simulation, it sets the parameter values (in log-scale for the RDM)
p_vector <- sampled_pars(designRDM)
p_vector["v"] = log(2) # intercept or base speed
p_vector["v_match:Sgo"] = log(2.5) # +- 0.5 match effect for GO stimuli
p_vector["v_match:Snogo"] = log(3) # +- 0.5 match effect for NOGO stimuli, here NOGO will have better discriminability
p_vector["B"] = log(1.4) # base threshold, here applies to the GO accumulator (because it's the first level of R)
p_vector["B_lRnogo"] = log(0.85) # *proportional* multiplier for NOGO threshold, i.e. <1 will be a lower threshold because of the log scale
p_vector["t0"] = log(.2) # shared t0
p_vector["s_lMTRUE"] = log(.8) # shared *multiplier* on within-trial diffusion variability for whichever accumulator matches the stimulus, here it will have lower variance as it is <1
# Call make_data to simulate, noting we then correct a couple of things until the bugs are patched
# That is, we set UC = 0 for nogo outcomes and Inf for go outcomes, and we ensure R=="nogo" for nogo outcomes (not NA as I previosuly indicated)
dat <- make_data(p_vector,designRDM, n_trials=1000) %>%
  mutate(R = factor(ifelse(is.na(rt),"nogo","go"),levels=c("go","nogo")),
         UC = ifelse(is.na(rt),0,Inf))
## Summarise our data to make sure it looks like go/nogo
tmp = tapply(is.na(dat$rt),dat$S,mean)
cat("Number of unobserved RTs per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
tmp = tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
cat("Median observed RT per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))

## This code sets slightly more informative priors
# pvec for the group mean value
# svec for the standard deviation of that (leaviung most at the default 1)
prior_pvec = sampled_pars(designRDM);
prior_svec = sampled_pars(designRDM)+1; # +1 so it doesn't accidentally set 0 variance on the prior and crash
prior_pvec[grepl("^v",names(prior_pvec))] = log(2) # set the intercept term at a non-zero value
prior_pvec[grepl("^v_",names(prior_pvec))] = 0 # terms other than the intercept are *differences* so set them at 0 by default
prior_pvec[grepl("^B",names(prior_pvec))] = log(1) # set the intercept term for threshold nonzero
prior_pvec[grepl("^B_",names(prior_pvec))] = 0 # and keep the difference terms centered on zero
prior_pvec[grepl("^t0",names(prior_pvec))] = log(.3) # set the prior scale for t0 much closer to what we actually expect
prior_svec[grepl("^t0",names(prior_pvec))] = .5 # set the piror width for t0 much smaller than default, especially for nogo data where it is less informed

## Make the emc object
# by default it will use rt_resolution = 1/60 and compress the data, which gives a speed boost and doesn't seem to reduce the quality of estimation
emc <- make_emc(dat,designRDM,type="single", # NB "single" here is because this is a single-subject simulation, use "standard" (which is the default) for normal hierarchical fits
                prior = prior(designRDM,theta_mu_mean=prior_pvec, # we could make a prior object but here I just define the prior within the make_emc call using the above vectors
                              theta_mu_sd = prior_svec,type="single"))

# Fit emc object (note, make sure cores_per_chain * cores_for_chains doesn't exceed what your computer has!)
# and cores_per_chain can't be used if you're on windows so just remove it if that's the case
emc <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_sample_iter = 1000
  ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
)

# This is purely for simulation, it plots the recovered parameters against the true generating ones - useful for checking designs actually recover
# Ideally the points roughly fall on the AB line
print(recovery(emc,p_vector,selection="alpha"))

## Generate some posterior predictive data
# here we've got to hack the emc object a little until the bugs get fixed
# because otherwise UC=0 for all the nogo rows and it won't be able to predict anything different
for (i in 1:length(emc)) {
  for (k in names(emc[[i]]$data)) {
    emc[[i]]$data[[k]]$UC = 999
  }
}
synth.data = predict(emc,n_post = 50) %>%
  mutate(R = factor(ifelse(is.na(rt),"nogo","go"),levels=c("go","nogo")),
         UC = ifelse(is.na(rt),0,Inf))

## The current plot_cdf function doesn't handle the missing data properly so we can reconstruct it here
plot_cdf(dat,synth.data,defective_factor = "R", factors = "S", remove_na = FALSE)

## Summarise our data to make sure it looks like go/nogo
tmp = tapply(is.na(synth.data$rt),synth.data$S,mean)
cat("Number of unobserved RTs per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
tmp = tapply(synth.data$rt,synth.data$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
cat("Median observed RT per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
