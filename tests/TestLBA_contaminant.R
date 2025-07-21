#### RACE LBA ----
devtools::load_all(reset = TRUE)
## Script to simulate from, and test recovery of, LBASWTN
library(EMC2)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
# if (identical(Sys.getenv("PAR_DEBUG"), "1")) {      # opt-in via env-var
#   env <- asNamespace("parallel")
#   unlockBinding("parLapply", env)                   # temporarily unlock
#   assign("parLapply", function(cl, X, FUN, ...) {
#     lapply(X, FUN, ...)                             # <- runs in master R
#   }, envir = env)
#   lockBinding("parLapply", env)
# }
## Set simulation parameters:
nsubs=1; ntrials=2000
## For maximum testing purposes we'll simulate from a model with varying n_acc
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designLBA <- design(
  factors=list(subjects=1:nsubs,S=c("yes"), block=c("low","medium","high")),
  Rlevels=c("yes"),
  matchfun=matchfun,
  model=LBAIO,constants=c(B=log(2),sv=log(1)),
  formula=list(v~0+block,B~1,t0~1,A~1,sv~1,pContaminant~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(3,2,1,log(0.2),log(0.5),qnorm(0.2))
mapped_pars(designLBA,p_vector)
# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designLBA,n_trials=ntrials,UC=4)
dat <- make_data(p_vector,designLBA,data=template)
p_names = names(sampled_pars(designLBA))

emc <- make_emc(dat,designLBA,type="single")
emc <- fit(emc,cores_for_chains = 3,fileName = 'samples.RData')
get_pars(emc,merge_chains=TRUE,return_mcmc = FALSE)[,,1000]
pred = predict(emc)
mean(is.na(dat$R[dat$block=="low"]));mean(is.na(pred$R[pred$block=="low"]))
mean(is.na(dat$R[dat$block=="high"]));mean(is.na(pred$R[pred$block=="high"]))
quantile(dat$rt[is.finite(dat$rt)&dat$block=="low"],seq(0.05,0.95,by=0.05), na.rm=TRUE)
quantile(pred$rt[is.finite(pred$rt)&pred$block=="low"],seq(0.05,0.95,by=0.05),na.rm=TRUE)
quantile(dat$rt[is.finite(dat$rt)&dat$block=="high"],seq(0.05,0.95,by=0.05),na.rm=TRUE)
quantile(pred$rt[is.finite(pred$rt)&pred$block=="high"],seq(0.05,0.95,by=0.05),na.rm=TRUE)
recovery(emc,p_vector)


