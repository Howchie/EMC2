#### RACE LBA ----
devtools::load_all()
devtools::document()
roxygen2::roxygenise()
devtools::install(upgrade = "never")  # rebuild & install
## then, in every run
library(EMC2)          # now workers load the same code automatically
#library(tictoc)
#source("tests/test_likelihood_plotfuns.R")
Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
source("forceSerial_TEST.R")         # your breakpoints now trigger
#RNGkind("L'Ecuyer-CMRG")
#set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBA,constants=c(v_RACE3=0,sv=log(1)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,s~1,cv~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(1), log(2), log(1),log(2),log(0.2),log(0.5))

# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designLBA,n_trials=1000)
#attr(template,"UC")=Inf
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(p_vector,designLBA,data=template)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R) | (d$R=="pm" & as.numeric(d$S)>2)
tapply(Cfun(dat),dat[,c("S","RACE")],mean)
# dadm <- EMC2:::design_model(dat,designLBA,compress=FALSE)

# Check likelihood
dadmLBA <- EMC2:::design_model(dat,designLBA)
dadmRDM <- EMC2:::design_model(dat,designRDM)
pars <- EMC2:::get_pars_matrix(p_vector, dadmLBA, model = attr(dadmLBA, "model")())

library(parallel)
# tic()
#profile_plot_test(dat,designRDM,p_vector,n_cores=1,layout=c(2,3)) # good
#profile_plot_test(dat,designRDM,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
#profile_plot_test(dat,designLBA,c(p_vector),n_cores=1,layout=c(2,3)) # good
#profile_plot_test(dat,designLBA,c(p_vector),n_cores=1,layout=c(2,3),use_c=TRUE) # ?
# toc()
## Set reduced model to test
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBA,constants=c(v_RACE3=0,s=log(1),A=log(0.5)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,s~1,cv~1),
)
p_names = names(sampled_pars(designLBA))
## Test setting priors using new reverse_transform logic
# Here we set our prior means
p_vector <- sampled_pars(designLBA,doMap=FALSE)
p_vector[grepl("B", names(p_vector))] <- 1
p_vector[grepl("t0", names(p_vector))] <- 0.3
p_vector[grepl("A", names(p_vector))] <- 1
p_vector[grepl("cv", names(p_vector))] <- 0.5
p_vector[regexpr("^v.*TRUE",names(p_vector))==1] <- 2
p_vector[regexpr("^v.*FALSE",names(p_vector))==1] <- 1
p_vector[names(p_vector)=="v"]<-1
# variances
s_vector <- sampled_pars(designLBA,doMap=FALSE)
s_vector[grepl("B", names(p_vector))] <- 1
s_vector[grepl("t0", names(p_vector))] <- 0.2
s_vector[grepl("A", names(p_vector))] <- 0.5
s_vector[grepl("cv", names(p_vector))] <- 0.25
s_vector[regexpr("^v.*TRUE",names(p_vector))==1] <- 1
s_vector[regexpr("^v.*FALSE",names(p_vector))==1] <- 1
s_vector[which(names(s_vector)=="v")]<-1
prior_pars=do_reverse_transform_variance(p_vector,s_vector,designLBA$model(),FALSE)
priorLBA = prior(designLBA,mu_mean=prior_pars$pars,mu_sd=sqrt(prior_pars$var),type="single")

emc <- make_emc(dat,designLBA,type="single",prior=priorLBA)
emc <- fit(emc,cores_for_chains = 3,fileName = 'samples.RData')
#recovery(emc,p_vector)
plot_pars(emc)



t_adj=1.328; alpha=2; mu_drift=1.5; sigma_drift=9.582; v=sigma_drift^2
denom = sqrt(t_adj*(1+t_adj*v))
prob_d_gt_0 = pnorm(mu_drift/sigma_drift)
rho = (sigma_drift*t_adj)/denom
h1 =  (mu_drift*t_adj - alpha)/denom
k1 =  mu_drift / sigma_drift
term1 = pbvn_drezner(h1, k1, rho)
mu_p = mu_drift + 2*alpha*v
h2 = (-mu_p*t_adj - alpha)/denom
k2 = mu_p/sigma_drift
term2 = pbvn_drezner(h2, k2, -rho)
cdf_val  = (term1 + exp(2*alpha*mu_drift + 2*(alpha^2)*v)*term2) / prob_d_gt_0
