
rm(list=ls())
library(EMC2)

# Here is a simple 5 parameter RDM model with contamination
designRDM <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=RDM,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(log(2), log(1.5), log(2), log(0.2), qnorm(.05))


### TRUNCATION + CENSORING + CONTAMINATION ----

TC=list(LT=.5,UT=1.3,LC=.55,UC=1.2,verbose=TRUE)
datCT <- make_data(p_vector,designRDM,n_trials=10000,TC=TC)

# Fitting
emc <- make_emc(datCT,designRDM,type="single")
emc <- fit(
  emc,
  stop_criteria = list(
    sample = list(
      iter = 2000,
      max_gd = 1.10,
      max_flat_loc = 0.5,
      flat_selection = c("alpha", "subj_ll"),
      flat_p1 = 1/3,
      flat_p2 = 1/3,
      max_sample_iter = 5000
    )
  )
)

print(recovery(emc,p_vector,selection="alpha"))
