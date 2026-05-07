## Production EMC2 script for the Forstmann et al. (2008) perceptual data.
## This mirrors the TRDM specification used in the paper:
## - race model with an additional timing accumulator
## - evidence drift and threshold vary by speed/neutral/accuracy
## - timer threshold and timer onset are fixed
## - across-trial drift variability is fixed to zero

rm(list = ls())
library(EMC2)
library(dplyr)
source("Functions.R")

# Prevent nested threading inside forked workers (BLAS/OpenMP/etc.).
# This keeps effective parallelism at cores_for_chains * cores_per_chain.
limit_nested_threads <- function(n = 1L) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 1L) n <- 1L
  Sys.setenv(
    OMP_NUM_THREADS = n,
    OPENBLAS_NUM_THREADS = n,
    MKL_NUM_THREADS = n,
    BLIS_NUM_THREADS = n,
    VECLIB_MAXIMUM_THREADS = n,
    NUMEXPR_NUM_THREADS = n,
    GOTO_NUM_THREADS = n
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n)
    RhpcBLASctl::omp_set_num_threads(n)
  }
}
limit_nested_threads(1L)

Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
# load FDBNCRW2008 data for fitting
load("data/FDBNCRW2008.RData")
data_trdm = data %>%
  mutate(sa = factor(instruction, levels = c("speed", "neutral", "accuracy")),
         S  = factor(S, levels = 1:2, labels = c("left", "right")),
         R  = factor(data$R, levels = 1:3, labels = c("left", "right", "time")), # adding time instantiates the lR for time
         subjects = factor(subj),UT=1.5,LT=.25) %>% #
  rename(rt=RT) %>%
  select(-c(correct,subj,instruction))

matchfun <- function(d) as.numeric(d$lR) == as.numeric(d$S)

design_TRDM <- design(
  data = data_trdm,
  factors = list(Rlevels=c("left","right","time")),
  model = RDMSWTN(),
  matchfun = matchfun,
  functions = list(match = function(d){dplyr::case_when(d$lM==TRUE & !(d$lR=="time") ~ 0.5,
                                                        d$lM==FALSE & !(d$lR=="time") ~-.5,
                                                        d$lR=="time" ~ 0)},
                   E = function(d){ifelse(d$lR=="time",0,1)},
                   Time = function(d){ifelse(d$lR=="time",1,0)},
                   Resp = function(d){dplyr::case_when(d$lR=="left" ~ 0.5,
                                                       d$lR=="right" ~ -0.5,
                                                       d$lR=="time" ~ 0)}),
  formula = list(
    v ~ 0 + sa:E + match + sa:Time,
    B ~ 0 + sa:E + Resp + Time,
    t0 ~ 0 + (E:sa + Time),
    s ~ 0 + (E + Time),
    A ~ 1,
    sv ~ 1
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    lambda_g = log(0),
    lambda_k = log(0),
    "B_saspeed:E" = log(1),
    "B_Time" = log(1),
    t0_Time = log(0.05)
  ),
  pre_transform_terms = list(B=c("B_saspeed:E","B_saneutral:E","B_saaccuracy:E","B_Time"),
                             v=c("v_saspeed:E","v_saneutral:E","v_saaccuracy:E",
                                 "v_saspeed:Time","v_saneutral:Time","v_saaccuracy:Time")),
  #transform = list(func=c(B="identity",v="identity")),
  bound = list(minmax = cbind(v = c(-Inf,Inf))),
  report_p_vector = TRUE
)
## Set priors
pvec = sampled_pars(design_TRDM)

## Test using TRDM paper
pvec_paper = c(c(3.99),log(c(2.2, 1.89, 1.59, 1.66, 1.12, 0.86)), -.11, log(c(1.05, 1.08, .24, .29, .29, 1.33, .51)))
names(pvec_paper) = names(pvec)
#pvec_paper[grepl("^B_",names(pvec_paper))] = exp(pvec_paper[grepl("^B_",names(pvec_paper))])
#pvec_paper["B_Resp"]=-.11
mapped_pars(design_TRDM,pvec_paper)
pred = make_data(pvec_paper,design_TRDM,n_trials=500)
plot_cdf(data_trdm, post_predict=pred, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")

pvec[grepl("^v",names(pvec))] = log(2)
pvec["v_match"] = 3
pvec[grepl("^B",names(pvec))] = log(1)
pvec[grepl("^B_Resp",names(pvec))] = 0
pvec[grepl("^t0",names(pvec))] = log(.2)
pvec[grepl("^s",names(pvec))] = log(1)
svec = sampled_pars(design_TRDM)
svec[] = 1
svec["B_Resp"] = .25
svec[grepl("^v_match",names(svec))] = 1
svec[grepl("^t0",names(svec))] = .3
prior_TRDM = prior(design_TRDM, theta_mu_mean = pvec, theta_mu_sd = svec, type="standard")
plot(prior_TRDM, map=FALSE)
emc <- make_emc(data_trdm, design_TRDM, n_chains = 3, compress = TRUE, prior=prior_TRDM)

## This is the production fit used for comparison to the paper.
## The saved object is the same one consumed by the downstream comparison scripts.
fit_TRDM <- fit(emc, cores_for_chains=3,cores_per_chain = 3,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 1,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max.tries = 30, fileName  = "samples.RData")
pp_TRDM = predict(fit_TRDM)

pp_sum = pp_TRDM %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(Acc = mean(Correct), rt=mean(rt))


plot_cdf(data_trdm, post_predict=pp_TRDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data_trdm,pp_TRDM,factor="sa")

post_means = rowMeans(get_pars(fit_TRDM,return_mcmc = FALSE, merge_chains = TRUE)[,,drop=TRUE])
mapped_pars(design_TRDM,post_means)
mapped_pars(design_TRDM,pvec_paper)
save.image(file = "modelFits/fit_TRDM.RData")
pvec_fit = rowMeans(get_pars(fit_TRDM,merge_chains = TRUE, return_mcmc = FALSE)[,,drop=TRUE])

par(mfrow=c(2,2))
tmp1 = rwald(10000,a=(exp(0)+0.5*pvec_fit["B_Resp"])/exp(pvec_fit["s_E"]),m=(exp(pvec_fit["v_saspeed:E"]) + 
                                                                               0.5*pvec_fit["v_match"])/exp(pvec_fit["s_E"]),exp(pvec_fit["t0_E:saspeed"]))
tmp2 = rwald(10000,a=(exp(0)+0.5*pvec_paper["B_Resp"])/exp(pvec_paper["s_E"]),m=(exp(pvec_paper["v_saspeed:E"]) + 
                                                                                   0.5*pvec_paper["v_match"])/exp(pvec_fit["s_E"]),exp(pvec_paper["t0_E:saspeed"]))
plot(density(tmp2), col='blue')
lines(density(tmp1), col='red', lty='dashed')

tmp1 = rwald(10000,a=1/exp(pvec_fit["s_Time"]),m=exp(pvec_fit["v_saspeed:Time"])/exp(pvec_fit["s_Time"]),0.05)
tmp2 = rwald(10000,a=1/exp(pvec_paper["s_Time"]),m=exp(pvec_paper["v_saspeed:Time"])/exp(pvec_paper["s_Time"]),0.05)
plot(density(tmp2), col='blue', xlim=c(0,3))
lines(density(tmp1), col='red', lty='dashed')

tmp1 = rwald(10000,a=exp(pvec_fit["B_saneutral:E"]+0.5*pvec_fit["B_Resp"])/exp(pvec_fit["s_E"]),m=exp((pvec_fit["v_saspeed:E"]+
                                                                                                         0.5*pvec_fit["v_E:match"]))/exp(pvec_fit["s_E"]),0)
tmp2 = rwald(10000,a=(1+0.5*pvec_paper["B_Resp"])/pvec_paper["s_E"],m=(pvec_paper["v_saspeed:E"]+
                                                                         0.5*pvec_paper["v_E:match"])/pvec_paper["s_E"])
plot(density(tmp2), col='blue')
lines(density(tmp1), col='red', lty='dashed')
#
## Fit the Erlang Process model with exponential timing (constant hazard)
load("data/FDBNCRW2008.RData")
data = data %>%
  mutate(sa = factor(instruction, levels = c("speed", "neutral", "accuracy")),
         S  = factor(S, levels = 1:2, labels = c("left", "right")),
         R  = factor(data$R, levels = 1:2, labels = c("left", "right")), # adding time instantiates the lR for time
         subjects = factor(subj)) %>%
  rename(rt=RT)

design_Erlang1 <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(erlang_shape = 1, erlang_type = "local_guess"),
  matchfun = matchfun,
  functions = list(match = function(d){dplyr::case_when(d$lM==TRUE ~ .5,
                                                        d$lM==FALSE ~-.5)},
                   Resp = function(d){dplyr::case_when(d$lR=="left" ~ 0.5,
                                                       d$lR=="right" ~ -.5)}),
  Rlevels = c("left", "right"),
  formula = list(
    v ~ 0 + sa + match,
    B ~ 0 + sa + Resp,
    t0 ~ 0 + sa,
    s ~ 1,
    A ~ 1,
    sv ~ 1,
    lambda_g ~ 0 + sa
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    lambda_k = log(0),
    "B_saspeed" = log(1)
  ),
  report_p_vector = TRUE
)
## Set priors
pvec = sampled_pars(design_Erlang1)
pvec[grepl("^v",names(pvec))] = log(2)
pvec["v_match"] = log(1)
pvec[grepl("^B",names(pvec))] = log(2)
pvec[grepl("^t0",names(pvec))] = log(.2)
pvec[grepl("^s",names(pvec))] = log(1)
pvec[grepl("^lambda",names(pvec))] = log(1)
svec = sampled_pars(design_Erlang1)
svec[] = 1
svec[grepl("^t0",names(svec))] = .3
prior_Erlang1 = prior(design_Erlang1, theta_mu_mean = pvec, theta_mu_sd = svec, type="standard")
plot(prior_Erlang1, map=FALSE)
emc <- make_emc(data, design_Erlang1, n_chains = 3, compress = TRUE, prior=prior_Erlang1)

## This is the production fit used for comparison to the paper.
## The saved object is the same one consumed by the downstream comparison scripts.

fit_Erlang1 <- fit(emc, cores_for_chains=3,cores_per_chain = 3,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max.tries = 30, fileName  = "samples2.RData")

pp_Erlang1 = predict(fit_Erlang1)
plot_cdf(data, post_predict=pp_Erlang1, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_Erlang1,factor="sa")
save.image(fit_Erlang1, file = "modelFits/fit_Erlang1.RData")

## Erlang2
design_Erlang2 <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(erlang_shape = 2, erlang_type = "local_guess"),
  matchfun = matchfun,
  functions = list(match = function(d){dplyr::case_when(d$lM==TRUE ~ .5,
                                                        d$lM==FALSE ~-.5)},
                   Resp = function(d){dplyr::case_when(d$lR=="left" ~ 0.5,
                                                       d$lR=="right" ~ -.5)}),
  Rlevels = c("left", "right"),
  formula = list(
    v ~ 0 + sa + match,
    B ~ 0 + sa + Resp,
    t0 ~ 0 + sa,
    s ~ 1,
    A ~ 1,
    sv ~ 1,
    lambda_g ~ 0 + sa
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    lambda_k = log(0),
    "B_saspeed" = log(1)
  ),
  report_p_vector = TRUE
)
## Set priors
pvec = sampled_pars(design_Erlang2)
pvec[grepl("^v",names(pvec))] = log(2)
pvec["v_match"] = log(1)
pvec[grepl("^B",names(pvec))] = log(1)
pvec[grepl("^t0",names(pvec))] = log(.2)
pvec[grepl("^s",names(pvec))] = log(1)
pvec[grepl("^lambda",names(pvec))] = log(1)
svec = sampled_pars(design_Erlang2)
svec[] = 1
svec[grepl("^t0",names(svec))] = .3
prior_Erlang2 = prior(design_Erlang2, theta_mu_mean = pvec, theta_mu_sd = svec, type="standard")
plot(prior_Erlang1, map=FALSE)
emc <- make_emc(data, design_Erlang2, n_chains = 3, compress = TRUE, prior=prior_Erlang2)

## This is the production fit used for comparison to the paper.
## The saved object is the same one consumed by the downstream comparison scripts.

fit_Erlang2 <- fit(emc, cores_for_chains=3,cores_per_chain = 3,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max.tries = 30, fileName  = "samples2.RData")

pp_Erlang2 = predict(fit_Erlang2)
plot_cdf(data, post_predict=pp_Erlang2, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_Erlang2,factor="sa")

save.image(file = "modelFits/fit_Erlang2.RData")

## RDM with no extra accumulator for comparison
design_RDM <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(),
  matchfun = matchfun,
  functions = list(match = function(d){dplyr::case_when(d$lM==TRUE ~ .5,
                                                        d$lM==FALSE ~-.5)},
                   Resp = function(d){dplyr::case_when(d$lR=="left" ~ 0.5,
                                                       d$lR=="right" ~ -.5)}),
  formula = list(
    v ~ 0 + sa + match,
    B ~ 0 + sa + Resp,
    t0 ~ 0 + sa,
    s ~ 1,
    A ~ 1,
    sv ~ 1
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    lambda_g = log(0),
    lambda_k = log(0),
    "B_saspeed" = log(1)
  ),
  report_p_vector = TRUE
)
## Set priors
pvec = sampled_pars(design_RDM)
pvec[grepl("^v",names(pvec))] = log(2)
pvec["v_match"] = log(1)
pvec[grepl("^B",names(pvec))] = log(1)
pvec[grepl("^t0",names(pvec))] = log(.2)
pvec[grepl("^s",names(pvec))] = log(1)
svec = sampled_pars(design_RDM)
svec[] = 1
svec[grepl("^t0",names(svec))] = .3
prior_RDM = prior(design_RDM, theta_mu_mean = pvec, theta_mu_sd = svec, type="standard")
plot(prior_RDM, map=FALSE)
emc <- make_emc(data, design_RDM, n_chains = 3, compress = TRUE, prior=prior_RDM)

## This is the production fit used for comparison to the paper.
## The saved object is the same one consumed by the downstream comparison scripts.
fit_RDM <- fit(emc, cores_for_chains=3,cores_per_chain = 3,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 1,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max.tries = 30, fileName  = "samples.RData")
pp_RDM = predict(fit_RDM)

plot_cdf(data, post_predict=pp_RDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_RDM,factor="sa")

save.image(file = "modelFits/fit_RDM.RData")

## Plot all fits
pdf("FitPlots.pdf")
plot_cdf(data_trdm, post_predict=pp_TRDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("TRDM",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_Erlang1, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("Erlang 1",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_Erlang2, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("Erlang 2",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_RDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("RDM (no T)",outer=TRUE,padj=2,side=3)
