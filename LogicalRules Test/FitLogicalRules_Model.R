## Script to simulate from, and test recovery of, RDMSWTN
rm(list=ls())
Packages <- c("EMC2","dplyr","sft")
lapply(Packages, library, character.only = TRUE) # to install change 'library' to 'install.packages'

## Read in Eidels 2015 data
load("LogicalRules Test/Eidels2015.RData")
data = all_data_acc
data = all_data_acc[all_data_acc$subject=="AW", ]

source("LogicalRules Test/Functions.R")
source("LogicalRules Test/Contrasts.R")
designBushmakin <- design(
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=correctfun,
  model=function(){LogicalRulesLBA(posdrift=TRUE)},constants=c('sv_lMFALSE:LogicalRuleAND'=log(1),'sv_lMFALSE:LogicalRuleOR'=log(1),p_LogicalRuleAND=qnorm(1),q_LogicalRuleAND=qnorm(0.5)),
  formula=list(v~0+Stim,B~(0+Resp:LogicalRule),t0~0+LogicalRule,A~0+LogicalRule,sv~0+lM:LogicalRule,p~0+LogicalRule,q~0+LogicalRule),
  functions=list(Stim=RateYNMatchMismatch_Contrast,Resp=YN_Contrast,nDots=nTargets_Contrast),
  factors=list(subjects=factor(rep(1,8)),S=factor(rep(c("DT","UT","LT","NT"),2),levels=c("DT","UT","LT","NT")),LogicalRule=factor(rep(c("OR","AND"),each=4),levels=c("AND","OR"))),
  Rlevels=c("yes","no")
)

BushmakinModel = make_emc(data,designBushmakin,type="single",fileName = 'samples.RData')
BushmakinModel=fit(BushmakinModel,max_tries = 50, cores_per_chain=1, cores_for_chains=3, 
    stop_criteria = list(
    preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
    sample = list(iter = 25, max_gd = 2), verbose = FALSE, particle_factor = 30, step_size = 25))
save.image("BushmakinModel.RData")
tmp=predict(BushmakinModel); pars = get_pars(BushmakinModel,merge_chains=TRUE,return_mcmc = FALSE)[,1000]; acc = mean(correctfun(BushmakinModel))