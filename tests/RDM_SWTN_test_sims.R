## Script to simulate from, and test recovery of, RDMSWTN
library(EMC2)
devtools::load_all(reset = TRUE)
RNGkind("L'Ecuyer-CMRG")
source("tests/test_likelihood_plotfuns.R")
set.seed(123)
nsubs=1; ntrials=2000
## For maximum testing purposes we'll simulate from a model with varying n_acc
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designRDMSWTN <- design(
  factors=list(subjects=1:nsubs,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=RDMSWTN,constants=c(s=log(1)),
  formula=list(v~0+lM,B~1,t0~1,zA~1,s~1,cv~1),
)

p_vector <- sampled_pars(designRDMSWTN,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(1), log(2),log(2),log(0.2),qnorm(0.6),qnorm(0.3))
mapped_pars(designRDMSWTN,p_vector)

dat <- make_data(p_vector,designRDMSWTN,n_trials=ntrials)

profile_plot_test(dat,designRDMSWTN,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dat,designRDMSWTN,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
emc <- make_emc(dat,designRDMSWTN,type="single")
emc <- fit(emc,cores_for_chains = 1,cores_per_chain=1,fileName = 'samples.RData')
mapped_pars(emc)
tmp=predict(emc); tmp=tmp[tmp$rt<5,]; dat=dat[dat$rt<5,]
q1 = quantile(tmp$rt); q2=quantile(dat$rt)
recovery(emc,p_vector)
plot_pars(emc,layout = c(1,1))

