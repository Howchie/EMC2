rm(list = ls())
library(EMC2)
library(dplyr)
source("Functions.R")
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)

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
  model = RDMSWTN(posdrift=FALSE),
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
    "B_saspeed:E" = log(1),
    "B_Time" = log(1),
    t0_Time = log(0.05)
  )
)

## Test using TRDM paper
#pvec_paper = c(c(3.99),log(c(2.2, 1.89, 1.59, 1.66, 1.12, 0.86)), -.11, log(c(1.05, 1.08, .24, .29, .29, 1.33, .51)))

prior_TRDM = make_priors(design_TRDM)
emc <- make_emc(data_trdm, design_TRDM, n_chains = 3, compress = TRUE, prior=prior_TRDM)

fit_TRDM <- fit(emc, cores_for_chains=3,cores_per_chain = 3,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = .5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max_tries = 50, fileName  = "samples.RData",rhat_version = "new", save_pw_ll = TRUE)
pp_TRDM = predict(fit_TRDM,cores=4)

summary_TRDM = pp_TRDM %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(pred_Acc = mean(Correct),
            pred_Rt_C = mean(rt[Correct]),
            pred_RT_E = mean(rt[!Correct]),
            TimedRespErrors = mean(isTime[!Correct]),
            TimedRespCorrect = mean(isTime[Correct])) %>%
  left_join(data_trdm %>%
              mutate(Correct=Cfun(.)) %>%
              group_by(sa) %>%
              summarise(dat_Acc = mean(Correct),
                        dat_Rt_C = mean(rt[Correct]),
                        dat_RT_E = mean(rt[!Correct])))

plot_cdf(data_trdm, post_predict=pp_TRDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data_trdm,pp_TRDM,factor="sa")

save.image(file = "modelFits/fit.RData")

## Fit the Erlang Process model with exponential timing (constant hazard)
load("data/FDBNCRW2008.RData")
data = data %>%
  mutate(sa = factor(instruction, levels = c("speed", "neutral", "accuracy")),
         S  = factor(S, levels = 1:2, labels = c("left", "right")),
         R  = factor(data$R, levels = 1:2, labels = c("left", "right")), # adding time instantiates the lR for time
         subjects = factor(subj), LT=.25, UT=1.5) %>%
  rename(rt=RT) %>%
  select(-c(correct,subj,instruction))

design_Erlang1 <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(erlang_shape = 1, erlang_type = "local_guess", posdrift=FALSE),
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
    mG ~ 0 + sa
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    "B_saspeed" = log(1)
  )
)
## Set priors
prior_Erlang1 = make_priors(design_Erlang1)
emc <- make_emc(data, design_Erlang1, n_chains = 3, compress = TRUE, prior=prior_Erlang1)

fit_Erlang1 <- fit(emc, cores_for_chains=3,cores_per_chain = 4,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = .5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max_tries = 50, fileName  = "samples.RData",rhat_version = "new", save_pw_ll = TRUE)

pp_Erlang1 = predict(fit_Erlang1)

summary_Erlang1 = pp_Erlang1 %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(pred_Acc = mean(Correct),
            pred_Rt_C = mean(rt[Correct]),
            pred_RT_E = mean(rt[!Correct]),
            TimedRespErrors = mean(isTime[!Correct]),
            TimedRespCorrect = mean(isTime[Correct])) %>%
  left_join(data %>%
              mutate(Correct=Cfun(.)) %>%
              group_by(sa) %>%
              summarise(dat_Acc = mean(Correct),
                        dat_Rt_C = mean(rt[Correct]),
                        dat_RT_E = mean(rt[!Correct])))

plot_cdf(data, post_predict=pp_Erlang1, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_Erlang1,factor="sa")
save.image(file = "modelFits/fit.RData")

## Erlang2
design_Erlang2 <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(erlang_shape = 2, erlang_type = "local_guess", posdrift=FALSE),
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
    mG ~ 0 + sa
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    "B_saspeed" = log(1)
  )
)
## Set priors
prior_Erlang2 = make_priors(design_Erlang2)
emc <- make_emc(data, design_Erlang2, n_chains = 3, compress = TRUE, prior=prior_Erlang2)

fit_Erlang2 <- fit(emc, cores_for_chains=3,cores_per_chain = 4,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = .5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max_tries = 50, fileName  = "samples.RData",rhat_version = "new", save_pw_ll = TRUE)

pp_Erlang2 = predict(fit_Erlang2)

summary_Erlang2 = pp_Erlang2 %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(pred_Acc = mean(Correct),
            pred_Rt_C = mean(rt[Correct]),
            pred_RT_E = mean(rt[!Correct]),
            TimedRespErrors = mean(isTime[!Correct]),
            TimedRespCorrect = mean(isTime[Correct])) %>%
  left_join(data %>%
              mutate(Correct=Cfun(.)) %>%
              group_by(sa) %>%
              summarise(dat_Acc = mean(Correct),
                        dat_Rt_C = mean(rt[Correct]),
                        dat_RT_E = mean(rt[!Correct])))

plot_cdf(data, post_predict=pp_Erlang2, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_Erlang2,factor="sa")

save.image(file = "modelFits/fit.RData")

## Erlang Mixed
design_Erlang3 <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(erlang_shape = "mixed", erlang_type = "local_guess", posdrift=FALSE),
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
    mG ~ 0 + sa,
    omega ~ 0 + sa
  ),
  constants = c(
    A = log(0),
    sv = log(0),
    "B_saspeed" = log(1)
  )
)
## Set priors
prior_Erlang3 = make_priors(design_Erlang3)
emc <- make_emc(data, design_Erlang3, n_chains = 3, compress = TRUE, prior=prior_Erlang3)

fit_Erlang3 <- fit(emc, cores_for_chains=3,cores_per_chain = 4,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = .5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max_tries = 50, fileName  = "samples.RData",rhat_version = "new", save_pw_ll = TRUE)

pp_Erlang3 = predict(fit_Erlang3)

summary_Erlang3 = pp_Erlang3 %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(pred_Acc = mean(Correct),
            pred_Rt_C = mean(rt[Correct]),
            pred_RT_E = mean(rt[!Correct]),
            TimedRespErrors = mean(isTime[!Correct]),
            TimedRespCorrect = mean(isTime[Correct])) %>%
  left_join(data %>%
              mutate(Correct=Cfun(.)) %>%
              group_by(sa) %>%
              summarise(dat_Acc = mean(Correct),
                        dat_Rt_C = mean(rt[Correct]),
                        dat_RT_E = mean(rt[!Correct])))

plot_cdf(data, post_predict=pp_Erlang3, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_Erlang3,factor="sa")

save.image(file = "modelFits/fit.RData")

## RDM with no extra accumulator for comparison
design_RDM <- design(
  data = data,
  factors = list(Rlevels=c("left","right")),
  model = RDMSWTN(posdrift=FALSE),
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
    "B_saspeed" = log(1)
  )
)
## Set priors
prior_RDM = make_priors(design_RDM)
emc <- make_emc(data, design_RDM, n_chains = 3, compress = TRUE, prior=prior_RDM)

fit_RDM <- fit(emc, cores_for_chains=3,cores_per_chain = 4,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = .5,
    flat_selection = c("alpha", "subj_ll","theta_mu"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  )), max_tries = 50, fileName  = "samples.RData",rhat_version = "new", save_pw_ll = TRUE)
pp_RDM = predict(fit_RDM)

summary_RDM = pp_RDM %>%
  mutate(Correct=Cfun(.)) %>%
  group_by(sa) %>%
  summarise(pred_Acc = mean(Correct),
            pred_Rt_C = mean(rt[Correct]),
            pred_RT_E = mean(rt[!Correct]),
            TimedRespErrors = 0, TimedRespCorrect=0) %>%
  left_join(data %>%
              mutate(Correct=Cfun(.)) %>%
              group_by(sa) %>%
              summarise(dat_Acc = mean(Correct),
                        dat_Rt_C = mean(rt[Correct]),
                        dat_RT_E = mean(rt[!Correct])))

plot_cdf(data, post_predict=pp_RDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
plot_subj_estimates(data,pp_RDM,factor="sa")

save.image(file = "modelFits/fit.RData")

models_list = list(fit_TRDM,fit_Erlang1,fit_Erlang2, fit_Erlang3, fit_RDM)
model_comparison = compare(models_list, WAIC = TRUE, LOO = TRUE, cores=4,cores_per_prop = 4)
model_comparison_subj = compare(models_list, WAIC = TRUE, LOO = TRUE, pointwise = "subject", cores=4,cores_per_prop = 4)

## Plot all fits
pdf("FitPlots.pdf")
plot_cdf(data_trdm, post_predict=pp_TRDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("TRDM",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_Erlang1, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("Erlang 1",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_Erlang2, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("Erlang 2",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_Erlang3, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("Erlang 3",outer=TRUE,padj=2,side=3)
plot_cdf(data, post_predict=pp_RDM, functions=list(Correct=Cfun), defective_factor = "Correct", factors="sa")
mtext("RDM (no T)",outer=TRUE,padj=2,side=3)
dev.off()
save.image(file = "modelFits/fit.RData")