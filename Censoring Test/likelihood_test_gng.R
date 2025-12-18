## This script is written by Dr Zachary Howard (UWA)
# It demonstrates the Go/No-Go implementation in EMC2

library(EMC2)
# This package constructs and draws profile plots on the natural scale, centered on the true value
source("Censoring Test/test_likelihood_plotfuns.R")

## Construct a simulation using the new "LBAGNG" model
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designLBA <- design(
  factors=list(subjects=1,S=c("go","go","go","nogo"), Difficulty=c("Low","High")),
  Rlevels=c("go","nogo"),
  matchfun=matchfun,
  model=LBAGNG,
  formula = list(B ~ 0+S, v ~ 0+(lM:Difficulty), A ~ 1, sv ~ 0+lM, t0 ~ 1),
  constants = c(sv_lMFALSE = log(1))
)

## Set some snesible parameters for simulating
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2),log(3), 1, 3, 1.5, 2.5, log(1), log(.75),log(.2))
dat <- make_data(p_vector,designLBA,n_trials=1000, UC=3)

## Check trial numbers
table(dat$S, dat$Difficulty)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
Rtfun <- function(d){ d$rt[d$rt==Inf]=NA; d$rt}

## Print some sanity checks to show accuracy and RT by cell (showing there are errors)
print("Accuracy by Cell: ")
tapply(Cfun(dat),dat[,c("S","Difficulty")],function(x){mean(x,na.rm=TRUE)})
print("RT by Cell: ")
tapply(Rtfun(dat),dat[,c("S","Difficulty")],function(x){mean(x,na.rm=TRUE)})

## Sanity check that there are no "no-go" responses with an RT, and no "go" responses with an undefined RT
print("RT by Response: ")
tapply(Rtfun(dat),dat[,c("R")],function(x){mean(x,na.rm=TRUE)})
print("Proportion of Undefined RTs by Response: ")
tapply(Rtfun(dat),dat[,c("R")],function(x){mean(is.na(x))})

## Construct profile plots showing relative recovery
# R and C functions give identical output, C is much faster
library(parallel)
profile_plot_test(dat,designLBA,p_vector,n_cores=1,layout=c(3,3)) # R likelihood
profile_plot_test(dat,designLBA,p_vector,n_cores=1,layout=c(3,3),use_c=TRUE) # C likelihood

## Perform a parameter recovery fit
emc <- make_emc(dat,designLBA,type="single")
emc <- fit(emc,cores_for_chains = 3,fileName = 'samples.RData')

## Inspect output
post_predict = predict(emc, n_post=50)
plot_pars(emc,post_predict=post_predict, true_pars = p_vector)