#### RACE LBA ----
devtools::load_all(reset = TRUE)
## Script to simulate from, and test recovery of, LBASWTN
library(EMC2)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
if (identical(Sys.getenv("PAR_DEBUG"), "1")) {      # opt-in via env-var
  env <- asNamespace("parallel")
  unlockBinding("parLapply", env)                   # temporarily unlock
  assign("parLapply", function(cl, X, FUN, ...) {
    lapply(X, FUN, ...)                             # <- runs in master R
  }, envir = env)
  lockBinding("parLapply", env)
}
## Set simulation parameters:
nsubs=10; ntrials=200
## For maximum testing purposes we'll simulate from a model with varying n_acc
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designLBA <- design(
  factors=list(subjects=1:nsubs,S=c("yes"), block=c("low","high")),
  Rlevels=c("yes"),
  matchfun=matchfun,
  model=LBAIO,constants=c(sv=log(1)),
  formula=list(v~0+block,B~1,t0~1,A~1,sv~1,pContaminant~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(2,1,log(2),log(0.2),log(0.5),qnorm(0.2))
mapped_pars(designLBA,p_vector)
reLBA = make_random_effects(design = designLBA, group_means = p_vector,n_subj=nsubs, variance_proportion = 0.05)
# Make square data so can remove pm in RACE = 2
template <- make_data(reLBA,designLBA,n_trials=ntrials)
dat <- make_data(reLBA,designLBA,data=template)
p_names = names(sampled_pars(designLBA))
## Test setting priors using new reverse_transform logic
# Here we set our prior means
# designLBA <- design(
#   factors=list(subjects=1:nsubs,S=c("yes")),
#   Rlevels=c("yes"),
#   matchfun=matchfun,
#   model=LBA,constants=c(s=log(1)),
#   formula=list(v~1,B~1,t0~1,A~1,s~1),
# )
p_vector <- sampled_pars(designLBA,doMap=FALSE)
p_vector[grepl("B", names(p_vector))] <- 1
p_vector[grepl("t0", names(p_vector))] <- 0.3
p_vector[grepl("A", names(p_vector))] <- 1
p_vector[grepl("cv", names(p_vector))] <- 0.5
p_vector[regexpr("^v.*",names(p_vector))==1] <- 2
p_vector[regexpr("^v.*FALSE",names(p_vector))==1] <- 1
p_vector[names(p_vector)=="v"]<-1
p_vector[regexpr("pContaminant",names(p_vector))==1] <- 0.1
# variances
s_vector <- sampled_pars(designLBA,doMap=FALSE)
s_vector[grepl("B", names(s_vector))] <- 1
s_vector[grepl("t0", names(s_vector))] <- 1
s_vector[grepl("A", names(s_vector))] <- 1
s_vector[grepl("cv", names(s_vector))] <- 0.25
s_vector[regexpr("^v.*",names(s_vector))==1] <- 1
s_vector[regexpr("^v.*FALSE",names(s_vector))==1] <- 1
s_vector[regexpr("^v.*TRUE",names(s_vector))==1] <- 1
s_vector[regexpr("^v",names(s_vector))==1] <- 1
s_vector[regexpr("pContaminant",names(s_vector))==1] <- 0.1
prior_pars=do_reverse_transform_variance(p_vector,s_vector,designLBA$model(),FALSE)
priorLBA = prior(designLBA,mu_mean=prior_pars$pars,mu_sd=sqrt(prior_pars$var),type="standard")

emc <- make_emc(dat,designLBA,type="standard",prior=priorLBA)
emc <- fit(emc,cores_for_chains = 3,fileName = 'samples.RData')
#recovery(emc,p_vector)
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(2, 1,log(2),log(0.2),log(0.5))
plot_pars(emc,true_pars = p_vector)

