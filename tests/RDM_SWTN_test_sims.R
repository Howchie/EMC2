### Test implementation of RDMSWTN
devtools::load_all(reset = TRUE)
Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
source("forceSerial_TEST.R")         # your breakpoints now trigger
#Sys.unsetenv("PAR_DEBUG")     # remove when finished
## Script to simulate from, and test recovery of, RDMSWTN
library(EMC2)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

## Set simulation parameters:
nsubs=10; ntrials=500
## For maximum testing purposes we'll simulate from a model with varying n_acc
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designRDMSWTN <- design(
  factors=list(subjects=1:nsubs,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=RDMSWTN,constants=c(v_RACE3=0,s=log(1)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,s~1,cv~1),
)

p_vector <- sampled_pars(designRDMSWTN,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(1), log(1.5), log(2),log(2),log(0.2),log(0.5),log(0.14))
mapped_pars(designRDMSWTN,p_vector)
reRDMSWTN = make_random_effects(design = designRDMSWTN, group_means = p_vector,n_subj=nsubs, variance_proportion = 0.05)
# Make square data so can remove pm in RACE = 2
template <- make_data(reRDMSWTN,designRDMSWTN,n_trials=ntrials)
#attr(template,"UC")=Inf
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(reRDMSWTN,designRDMSWTN,data=template)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R) | (d$R=="pm" & as.numeric(d$S)>2)
print("Simulated Accuray Per Cell:")
tapply(Cfun(dat),dat[,c("S","RACE")],mean)
# dadm <- EMC2:::design_model(dat,designRDMSWTN,compress=FALSE)

p_names = names(sampled_pars(designRDMSWTN))
## Test setting priors using new reverse_transform logic
# Here we set our prior means
p_vector <- sampled_pars(designRDMSWTN,doMap=FALSE)
p_vector[grepl("B", names(p_vector))] <- 1
p_vector[grepl("t0", names(p_vector))] <- 0.3
p_vector[grepl("A", names(p_vector))] <- 1
p_vector[grepl("cv", names(p_vector))] <- 0.5
p_vector[regexpr("^v.*TRUE",names(p_vector))==1] <- 2
p_vector[regexpr("^v.*FALSE",names(p_vector))==1] <- 1
p_vector[names(p_vector)=="v"]<-1
# variances
s_vector <- sampled_pars(designRDMSWTN,doMap=FALSE)
s_vector[grepl("B", names(p_vector))] <- 1
s_vector[grepl("t0", names(p_vector))] <- 0.2
s_vector[grepl("A", names(p_vector))] <- 0.5
s_vector[grepl("cv", names(p_vector))] <- 0.25
s_vector[regexpr("^v.*TRUE",names(p_vector))==1] <- 1
s_vector[regexpr("^v.*FALSE",names(p_vector))==1] <- 1
s_vector[which(names(s_vector)=="v")]<-1
prior_pars=do_reverse_transform_variance(p_vector,s_vector,designRDMSWTN$model(),FALSE)
priorRDMSWTN = prior(designRDMSWTN,mu_mean=prior_pars$pars,mu_sd=sqrt(prior_pars$var),type="standard")

emc <- make_emc(dat,designRDMSWTN,type="standard",prior=priorRDMSWTN)
emc <- fit(emc,cores_for_chains = 1,cores_per_chain=1,fileName = 'samples.RData')
#recovery(emc,p_vector)
plot_pars(emc)

