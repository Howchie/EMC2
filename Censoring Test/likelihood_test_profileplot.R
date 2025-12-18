#### RACE LBA ----
devtools::load_all(reset = TRUE)
## then, in every run
#library(EMC2)          # now workers load the same code automatically
source("Censoring Test/test_likelihood_plotfuns.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
designLBA <- design(
  factors=list(subjects=1,S=c("go","nogo")),
  Rlevels=c("go","nogo"),
  matchfun=matchfun,
  model=LBAGNG,
  formula = list(B ~ 1, v ~ lM, A ~ 1, sv ~ lM, t0 ~ 1),
  constants = c(sv = 0)
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2),3,1,log(2),log(.75),log(.2))

# Make square data so can remove pm in RACE = 2
dat <- make_data(p_vector,designLBA,n_trials=1000, UC=3)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
Rtfun <- function(d){ d$rt[d$rt==Inf]=NA; d$rt}
tapply(Cfun(dat),dat[,c("S")],function(x){mean(x,na.rm=TRUE)})
tapply(Rtfun(dat),dat[,c("S")],function(x){mean(x,na.rm=TRUE)})

# Check likelihood
dadmLBA <- EMC2:::design_model(dat,designLBA)
pars <- EMC2:::get_pars_matrix(p_vector, dadmLBA, model = attr(dadmLBA, "model")())

#library(parallel)
profile_plot_test(dat,designLBA,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dat,designLBA,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
emc <- make_emc(dat,designLBA,type="single")
emc <- fit(emc,cores_for_chains = 3,fileName = 'samples.RData')
#data=dat;design=designLBA;n_point=100;range=0.5;dadm=design_model(data, design, verbose = FALSE);min_ll=log(1e-10)
#model=attr(dadm, "model")();pars <- get_pars_matrix(p_vector, dadm, model)
