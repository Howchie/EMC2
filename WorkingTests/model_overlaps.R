rm(list=ls())
library(EMC2)

matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designRDM <- design(
  factors=list(subjects=1,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=RDM,constants=c(s=log(1)),
  formula=list(v~lM,B~1,t0~1,A~1)
)
designRDMSWTN <- design(
  factors=list(subjects=1,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=RDMSWTN,constants=c(s=log(1),A=log(0)),
  formula=list(v~lM,B~1,t0~1,s~1,sv~1,A~1)
)
designRDMGBM <- design(
  factors=list(subjects=1,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=RDMGBM,constants=c(s=log(1)),
  formula=list(v~lM,B~1,t0~1,A~1)
)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=BAwL,constants=c(sv=log(1)),
  formula=list(v~lM,B~1,t0~1,A~1)
)
designFRQ <- design(
  factors=list(subjects=1,S=c("left","right")),
  Rlevels=c("left","right"),
  matchfun=matchfun,
  model=FRQ,constants=c(h=qnorm(1),delta=log(0)),
  formula=list(lambda~lM,alpha~1,beta~1,h~1,delta~1,cv_u~1,t0~1)
)

p_vec <- sampled_pars(designRDMSWTN,doMap = FALSE)
p_vec[1:length(p_vec)] <- c(log(1.5), log(2), log(1), log(0.5),log(.3))
dat = make_data(p_vec,designRDMSWTN, n_trials = 5000)

emc <- make_emc(dat, designRDM, type = "single", compress = T)
emcRDM <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)

#print(recovery(emcRDM,p_vec,selection="alpha"))

emc <- make_emc(dat, designRDMSWTN, type = "single", compress = T)
emcRDMSWTN <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)
print(recovery(emcRDMSWTN,p_vec,selection="alpha"))

# emc <- make_emc(dat, designRDMGBM, type = "single", compress = T)
# emcRDMGBM <- fit(emc,stop_criteria = list(
#   sample = list(
#     iter = 1000,
#     max_gd = 1.10,
#     max_flat_loc = 0.5,
#     flat_selection = c("alpha", "subj_ll"),
#     flat_p1 = 1/3,
#     flat_p2 = 1/3,
#     max_sample_iter = 5000
#   ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)

emc <- make_emc(dat, designLBA, type = "single", compress = T)
emcLBA <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)
# print(recovery(emcRDMGBM,p_vec,selection="alpha"))
emc <- make_emc(dat, designFRQ, type = "single", compress = T)
emcFRQ <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)
comp = compare(sList=list(emcRDM,emcRDMSWTN,emcLBA,emcFRQ),BayesFactor=FALSE,cores_for_loo = 8)
pred=predict(emcFRQ,n_cores = 24)

plot_cdf(dat,pred,defective_factor = "R", factors = "S")
gbm_pars = get_pars(emcRDMGBM,merge_chains = TRUE, return_mcmc = FALSE)[,,,drop=TRUE]
swtn_pars = get_pars(emcRDMSWTN,merge_chains = TRUE, return_mcmc = FALSE)[,,,drop=TRUE]
