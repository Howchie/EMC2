## Script to simulate and test recovery of Go/No-go Race Models
# Using the RDM for now as it is less unidentified than LBA (it seems to be the combined ridge of A + B + t0 that hurts most)
# Also exploring the "marginalize" option which integrates across the uncertainty of t0.
rm(list=ls())

for (lib in c("EMC2","EMC2ct")) {
  library(lib,character.only = TRUE)
  set.seed(123)
  ## Test 1 - Shared base drift, match x S, threshold by S, single t0. 
  designRDM <- design(
    factors=list(subjects=1,S=c("go","no-go")),Rlevels=c("go","no-go"),
    matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
    functions=list(match=function(d) ifelse(d$lM==TRUE,.5,-.5)),
    model=RDM,TC = list(UC=3),
    formula=list(v~match:S,A~1,B~lR,t0~1,s~lM),
    constants = c(s=log(1),A=log(0))
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
  emc <- make_emc(dat,designRDM,type="single",
                  prior = prior(designRDM,theta_mu_mean=prior_pvec,
                                 theta_mu_sd = prior_svec,type="single"))
  system.time({emc1 <- fit(emc,stop_criteria = list(
    sample = list(
      iter = 1000,
      max_gd = 1.10,
      max_sample_iter = 1000
    ),cores_per_chain=6, cores_for_chains = 3), max_tries=30,
    )})
  print(recovery(emc1,p_vector,selection="alpha"))
  unloadNamespace(lib)
}