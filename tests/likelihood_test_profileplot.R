#### RACE LBA ----
devtools::load_all(reset = TRUE)
devtools::document()
roxygen2::roxygenise()
devtools::install(upgrade = "never")  # rebuild & install
## then, in every run
#library(EMC2)          # now workers load the same code automatically
source("tests/test_likelihood_plotfuns.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBAIO,constants=c(v_RACE3=0,sv=log(1)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,sv~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(0.5), log(4), log(1),log(1.5),log(0.2),log(0.5))

# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designLBA,n_trials=1000,UT=3)
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(p_vector,designLBA,data=template)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R) | (d$R=="pm" & as.numeric(d$S)>2)
tapply(Cfun(dat),dat[,c("S","RACE")],mean)

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
