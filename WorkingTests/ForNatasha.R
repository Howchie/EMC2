## Script to simulate and test recovery of Go/No-go Race Models
# Install EMC2 from this branch until bugs are fixed
#remotes::install_github("https://github.com/ampl-psych",ref="gng")

library(EMC2)
library(dplyr)
rm(list=ls())
set.seed(123)
## Note - the issue when we chatted was setting TC=list(UC=3) globally.
# It seems like for the current build there's a couple of unique quirks.
# For "nogo" outcomes, UC needs to be set to 0, so that it then integrates across the entire time range
# (previously it was trying to integrate from [UC, Inf])
# The complication there means we should *not* use any UC for the other trials because it won't do anything anyway
# i.e. for now, every "nogo" is considered a "nogo accumulator won". 
# The fixed branch handles this aytomatically for you.
# In future I assume this will change
# NB for using real data, drop the factors argument and just add something like data = data
# ensuring your data has columns with *factors* S, subjects, R (levels includes "nogo" but responses coded as NA)
designRDM <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)), # match here is identical to the ADMat lM from the help file, either way of specifying it works
  model=RDM,
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
dat <- make_data(p_vector,designRDM, n_trials=1000)
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
synth.data = predict(emc,n_post = 50)

## Plot the defective CDFs, which should show how well the data predicts both the observed RT shape and the overall "nogo" rates
plot_cdf(dat,synth.data,defective_factor = "R", factors = "S", remove_na = FALSE)

## Summarise our data to make sure it looks like go/nogo
tmp = tapply(is.na(synth.data$rt),synth.data$S,mean)
cat("Number of unobserved RTs per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
tmp = tapply(synth.data$rt,synth.data$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
cat("Median observed RT per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))

## Small example of fitting a mis-specified model and doing the model comparison
# We will fit a model with no stimulus drift difference (although we know we had one in the generating data)
designRDM_threshold <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)), # match here is identical to the ADMat lM from the help file, either way of specifying it works
  model=RDM,
  formula=list(v~match,A~1,B~lR,t0~1,s~lM),
  constants = c(s=log(1),A=log(0)) 
  # I set no start-point variability for RDM because it often doesn't need it, 
  # but if you want it just remove the A=log(0) from constants and it will estimate
  # a value instead
)
## Make the emc object
# by default it will use rt_resolution = 1/60 and compress the data, which gives a speed boost and doesn't seem to reduce the quality of estimation
emc_nodrift <- make_emc(dat,designRDM_threshold,type="single", # NB "single" here is because this is a single-subject simulation, use "standard" (which is the default) for normal hierarchical fits
                prior = prior(designRDM_threshold,theta_mu_mean=prior_pvec, # we could make a prior object but here I just define the prior within the make_emc call using the above vectors
                              theta_mu_sd = prior_svec,type="single"))

# Fit emc object (note, make sure cores_per_chain * cores_for_chains doesn't exceed what your computer has!)
# and cores_per_chain can't be used if you're on windows so just remove it if that's the case
emc_nodrift <- fit(emc_nodrift,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_sample_iter = 1000
  ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
)

## Run the compare() function and hopefully the right model wins! 
sList = list("drift" = emc, "nodrift" = emc_nodrift)
comparison = compare(sList,BayesFactor = TRUE) # BF can be expensive
# My results:
# MD   wMD DIC  wDIC BPIC wBPIC EffectiveN meanD Dmean minD
# drift   692 0.852 661 0.947  667 0.918          7   654   648  647
# nodrift 695 0.148 666 0.053  672 0.082          6   661   655  655
# so the 'right' model was selected by all metrics, but in this toy example there wouldn't be a huge amount to distinguish the two models


## DDM version
# For the DDM the process is basically identical. It's currently going to assume an implicit lower boundary that is hit for all R="nogo" trials.
designDDM <- design(
  factors=list(subjects=1,S=c("go","nogo")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  functions=list(match=function(d) ifelse(d$S=="go",-1,1)), # match here is used to keep drift on a positive scale. The DDM model uses the lower boundary for the *first* Rlevel, so the first drift gets multiplied by -1
  model=DDM,
  formula=list(v~0 + match:S,Z~1,a~1,t0~1,s~1,sv~1), # 0 + creates an absolute term per S rather than a contrast for easier interpretation in the DDM
  constants = c(s=log(1)) 
)

## This code is only used for the simulation, it sets the parameter values (drift on natural scale, others in log-scale or qnorm() for Z)
p_vector <- sampled_pars(designDDM)
p_vector["v_match:Sgo"] = .9 # match drift for GO stimuli
p_vector["v_match:Snogo"] = 1.25 # match drift for NOGO stimuli, here NOGO will have better discriminability *and* faster speed
p_vector["a"] = log(1.2) # boundary separation
p_vector["Z"] = qnorm(0.4) # *proportional* start point bias toward the upper boundary, which is the *second* level of R (at 0.4 we're slightly biasing towards "go" responses here)
p_vector["t0"] = log(.2) # shared t0
p_vector["sv"] = log(.2) # shared drift variability

# That is, we set UC = 0 for nogo outcomes and Inf for go outcomes, and we ensure R=="nogo" for nogo outcomes (not NA as I previosuly indicated)
dat <- make_data(p_vector,designDDM, n_trials=1000)
## Summarise our data to make sure it looks like go/nogo
tmp = tapply(is.na(dat$rt),dat$S,mean)
cat("Number of unobserved RTs per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
tmp = tapply(dat$rt,dat$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
cat("Median observed RT per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))

## This code sets slightly more informative priors
# pvec for the group mean value
# svec for the standard deviation of that (leaving most at the default 1)
prior_pvec = sampled_pars(designDDM);
prior_svec = sampled_pars(designDDM)+1; # +1 so it doesn't accidentally set 0 variance on the prior and crash
prior_pvec[grepl("^v_",names(prior_pvec))] = log(2) 
prior_pvec[grepl("^a",names(prior_pvec))] = log(1) # set the intercept term for boundary separation
prior_pvec[grepl("^Z",names(prior_pvec))] = qnorm(0.5) # and keep the start point bias centered on zero
prior_pvec[grepl("^t0",names(prior_pvec))] = log(.3) # set the prior scale for t0 much closer to what we actually expect
prior_svec[grepl("^t0",names(prior_pvec))] = .5 # set the piror width for t0 smaller than default, especially for nogo data where it is less informed
prior_svec[grepl("^Z",names(prior_pvec))] = .5 # set the piror width for Z smaller than default
prior_svec[grepl("^sv",names(prior_pvec))] = .5 # set the piror width for sv smaller than default

## Make the emc object
# by default it will use rt_resolution = 1/60 and compress the data, which gives a speed boost and doesn't seem to reduce the quality of estimation
emc <- make_emc(dat,designDDM,type="single", # NB "single" here is because this is a single-subject simulation, use "standard" (which is the default) for normal hierarchical fits
                prior = prior(designDDM,theta_mu_mean=prior_pvec, # we could make a prior object but here I just define the prior within the make_emc call using the above vectors
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

synth.data = predict(emc,n_post = 50)

## Plot the defective CDFs, which should show how well the data predicts both the observed RT shape and the overall "nogo" rates
plot_cdf(dat,synth.data,defective_factor = "R", factors = "S")


## Summarise our data to make sure it looks like go/nogo
tmp = tapply(is.na(synth.data$rt),synth.data$S,mean)
cat("Number of unobserved RTs per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
tmp = tapply(synth.data$rt,synth.data$S,function(x){mean(x[is.finite(x)],na.rm=TRUE)})
cat("Median observed RT per stimulus type:\n",paste(names(tmp),round(tmp,2),sep=":"))
