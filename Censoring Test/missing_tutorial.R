#remotes::install_github("https://github.com/Howchie/EMC2/",ref="censoring_truncation_ah")
unloadNamespace("EMC2")
devtools::install_local(force=TRUE)
rm(list=ls())
library(EMC2)
## This is a helper function to directly check likelihoods
lfun <- function(p_vector, dadm, use_c=T) {
  if (use_c) {
    p_matrix <- matrix(p_vector,nrow=1)
    colnames(p_matrix) <- names(p_vector)
    model <- attr(dadm, "model")()
    p_types=names(model$p_types)
    designs <- list()
    for (p in p_types) {
      designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
    }
    constants <- attr(dadm,"constants")
    if (is.null(constants)) constants <- NA
    EMC2:::calc_ll(p_matrix, dadm, constants,designs,model$c_name,
                   model$bound,model$transform,model$pre_transform,p_types,log(1e-10),model$trend)
  } else {
    EMC2:::calc_ll_R(p_vector, attr(dadm, "model")(), dadm)
  }
}
pdf("Censoring Test/Test_Plots.pdf")

#### Using make_missing to censor and truncate data ----

# Make a design with two subjects (for showcasing subjectwise missing functionality)
designRDM <- design(model = RDM,
  factors = list(subjects = 1:2, S = c("left", "right")),Rlevels = c("left", "right"),
  matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
  formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM, pContaminant~1),
  constants = c(s = log(1)))

# First address a case with no contamination.
p_vector <-c(log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8)),pContaminant=qnorm(0))

# Default data generation, LT/UT/LC/UC columns at default values
dat <- make_data(p_vector, designRDM,n_trials = 50)

# Make this into a typical data file with no LT/UT/LC/UC columns
data <- dat[,1:5]

# Apply make_missing, showing different ways the truncate and censor can be specified

# Single value lower truncation
LT = .7
# Subject-specific upper truncation
UT = c('1'=1.8,'2'=1.85)
# Lower truncation equal to the number of data rows
LC = rep(.75,nrow(data))
# Functional subject-specific upper truncation at 90th percentile
UC = \(d) {
  ok <- is.finite(d$rt)
  tapply(d$rt[ok],d$subjects[ok],quantile,probs=.9)
}

# Running this command with verbose TRUE reports the effects of filtering
dat1 <- make_missing(data,LT,UT,LC,UC,verbose=TRUE)

# NB: If data that has LT/UT/LC/UC columns make_missing does not change them
dat2 <- make_missing(dat1)
all(c(dat1$LC==dat2$LC,dat1$UC==dat2$UC,dat1$LC==dat2$LC,dat1$UC==dat2$UC))

# Now make a data file with a 3 second response window, so trials with model
# generated rt>3 are given rt = Inf and R = NA. This can be achieved in make_data
# by passing a list named TC with the arguments required my make_missing.
# We use a large number of trials so that there are reliably a few rts outside the
# window, and we up the default digits to get a precise report.
dat <- make_data(p_vector, designRDM,n_trials = 5000,
                 TC=list(UC=3,verbose=TRUE,digits=3))

# Lets look at the result
hist(dat$rt,breaks=seq(0,3,.025))

# Suppose we decide that trials with rt > 2.5 are caused by distraction, so we
# want to censor these rts, but we decide to keep the response. However, we
# don't want to affect the timeouts, so we use the no_censor argument, identifying
# the relevant trials with an argument.
dat1 <- make_missing(dat,UC=2.5,verbose=TRUE,UCresponse=TRUE,
                     no_censor = \(d) is.infinite(d$rt))

# We see there are now two types of UC
head(dat1[dat1$UC==3,],2)
head(dat1[dat1$UC==2.5,],2)

# Another way to make trials with R = NA and rt=Inf is via contamination (e.g.,
# as might be due to mind wandering), where their occurrence is random with
# probability. We will first do this through make_missing called by make_data
dat <- make_data(p_vector, designRDM,n_trials = 5000,
                 TC=list(pContaminant=.1,verbose=TRUE))

# We can achieve the same result through the model
p_vector["pContaminant"] <- qnorm(.1)
dat <- make_data(p_vector, designRDM,n_trials = 5000,
                 TC=list(verbose=TRUE))

# Note that the parameter setting is overridden if pContaminant is supplied
dat <- make_data(p_vector, designRDM,n_trials = 5000,
                 TC=list(pContaminant=.2,verbose=TRUE))

# NB: It is important that censoring and truncation match the rt_resolution
#     assumed when fitting (or more generally in evaluating likelihoods, as in
#     the next section). Previously rt_resolution could only be specified
#     through make_emc when making an emc object just prior to fitting. However,
#     to accommodate truncation and censoring, both make_missing and make data
#     have rt_resolution argument. For make_missing the default is the same as
#     for make_emc (1/60), whereas for make_data the default is no binning
#     (i.e., the same was was the case previously).

# Default no binning
dat <- make_data(p_vector, designRDM,n_trials = 2,TC=list(LT=.69,LC=.76,UC=1.59,UT=1.71))
dat

# Applying the default version of make_missing bins at the make_emc default of 1/60
make_missing(dat,LT=.69,LC=.76,UC=1.59,UT=1.71)

# This aligns the output of make_data with the make_emc default
make_data(p_vector, designRDM,n_trials = 2,
          TC=list(LT=.69,LC=.76,UC=1.59,UT=1.71,rt_resolution=1/60))


#### Check profiles ----

# Look at two cases, no contamination
pContaminant <- NULL; const <- c(pContaminant=qnorm(0))
# OR
# contamination
pContaminant <- qnorm(.1); const <- NULL

# Move back to a single-subject case
designRDM <- design(model = RDM,constants = c(s = log(1),const),
  factors = list(subjects = 1, S = c("left", "right")),Rlevels = c("left", "right"),
  matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
  formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM, pContaminant~1))

p_vector <- c(log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8)),pContaminant=pContaminant)

# Also define a model using Andrew's R
RDMR <- function() {
  m <- RDM()
  m$c_name=NULL
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m
}
designRDM1 <- design(model = RDMR,constants = c(s = log(1),const),
  factors = list(subjects = 1, S = c("left", "right")),Rlevels = c("left", "right"),
  matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
  formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM, pContaminant~1))

# Now lets look at speed and profile recovery, with three cases

# Either with only one value of bounds. This will be fastest, very little
# difference between C and R in speed and identical estimates.
LT=.75;UT=1.8;LC=.85;UC=1.6
# OR
# with 11 values each of truncate, now C 3x
LT <- \(d,lo=.7,hi=.8) round(runif(nrow(d),lo,hi),2)
UT <- \(d,lo=1.75,hi=1.85) round(runif(nrow(d),lo,hi),2)
LC=.85;UC=1.6
# OR
# with 11 values each of censor, not much slower
LT=.75;UT=1.8
UC <- \(d,lo=1.55,hi=1.65) round(runif(nrow(d),lo,hi),2)
LC <- \(d,lo=.8,hi=.9) round(runif(nrow(d),lo,hi),2)

# Make data
TC <- list(LT=LT,UT=UT,LC=LC,UC=UC,verbose=TRUE)
dat <- make_data(p_vector,designRDM,n_trials=5000,TC=TC)
c(nLT=length(unique(dat$LT)),nUT=length(unique(dat$UT)),
  nLC=length(unique(dat$LC)),nUC=length(unique(dat$UC)))

# Note that I reinstated a default use of C by plot_profile

# C profile
system.time({
  print(profile_plot(dat,designRDM,p_vector,n_cores=1,layout=c(2,4)))
})

# R profile, need to set use_c to False
system.time({
 print(profile_plot(dat,designRDM1,p_vector,n_cores=1,layout=c(2,4),use_c=FALSE))
})



#### RDM Parameter recovery ----

# Start with no contamination, truncation or censoring as a baseline comparison

# Here is a simple 4 parameter RDM model that uses Zach's C
designRDM <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=RDM,
  formula=list(v~lM,B~1,t0~1)
)

# This parameter vector produces reasonable RTs and accuracies with a small
# level of contamination.
p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(log(2), log(1.5), log(2), log(0.2))

## This is a hack to get access to my R
RDMR <- function() {
  m <- RDM()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}

# and the design that uses it.
designRDM1 <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=RDMR,
  formula=list(v~lM,B~1,t0~1)
)

# Here we set rt_resolution to match the default for make_emc
TC=list(LT=0,UT=Inf,LC=0,UC=Inf,rt_resolution=1/60)
dat <- make_data(p_vector,designRDM,n_trials=2000,TC=TC)
emc <- make_emc(dat,designRDM,type="single")
# # NB: we can see that rt_resolution matches up in the two cases
# head(dat)
# head(emc[[1]]$data[[1]][emc[[1]]$data[[1]]$lR=="left",])

# Fit C
system.time({emcc <- fit(emc)})
# Time difference of 17.11186 secs

# Fit R
emc <- make_emc(dat,designRDM1,type="single")
system.time({emcr <- fit(emc)})
# Time difference of 41.82084 secs

# Check recovery
print(recovery(emcc,p_vector,selection="alpha"))
print(recovery(emcr,p_vector,selection="alpha"))


### CONTAMINATION only ----

# I saved off generated data and fits
load("Censoring Test/RDM.RData")

# Here is a simple 5 parameter RDM model with contamination that uses Zach's C
designRDM <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=RDM,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[] <- c(log(2), log(1.5), log(2), log(0.2),qnorm(.05))


RDMR <- function() {
  m <- RDM()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}
designRDM1 <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=RDMR,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

# Make this if you want to repeat sampling
TC=list(LT=0,UT=Inf,LC=0,UC=Inf,rt_resolution=1/60)
dat <- make_data(p_vector,designRDM,n_trials=2000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(dat,designRDM)
dadmR <- EMC2:::design_model(dat,designRDM1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Fit C
emc <- make_emc(dat,designRDM,type="single")
system.time({emcc <- fit(emc)})
# Time difference of 1.29533 mins

# Fit R
emc <- make_emc(dat,designRDM1,type="single")
system.time({emcr <- fit(emc)})
# Time difference of 1.10279 mins

# Check recovery
print(recovery(emcc,p_vector,selection="alpha"))
print(recovery(emcr,p_vector,selection="alpha"))


### CENSORING + CONTAMINATION ----

# Use a bit simulation to get good results, but not slow due to compression
# I also did (slow) checks without compression with similar results.

TC=list(LC=.5,UC=1.2,verbose=TRUE,rt_resolution=1/60)
datC <- make_data(p_vector,designRDM,n_trials=10000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(datC,designRDM)
dadmR <- EMC2:::design_model(datC,designRDM1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Fitting
emc <- make_emc(datC,designRDM,type="single")
system.time({emccC <- fit(emc)})
# Time difference of 1.930384 mins
emc <- make_emc(datC,designRDM1,type="single")
system.time({emcrC <- fit(emc)})
# Time difference of 3.344983 mins

# Check recovery
print(recovery(emccC,p_vector,selection="alpha"))
print(recovery(emcrC,p_vector,selection="alpha"))


### TRUNCATION + CONTAMINATION ----

TC=list(LC=.5,UC=1.2,verbose=TRUE)
datT <- make_data(p_vector,designRDM,n_trials=10000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(datT,designRDM)
dadmR <- EMC2:::design_model(datT,designRDM1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Fitting
emc <- make_emc(datT,designRDM,type="single")
system.time({emccT <- fit(emc)})
# Time difference of 1.79773 mins
emc <- make_emc(datT,designRDM1,type="single")
system.time({emcrT <- fit(emc)})
# Time difference of 3.282718 mins

# Recovery
print(recovery(emccT,p_vector,selection="alpha"))
print(recovery(emcrT,p_vector,selection="alpha"))


### TRUNCATION + CENSORING + CONTAMINATION ----

TC=list(LT=.5,UT=1.3,LC=.55,UC=1.2,verbose=TRUE)
datCT <- make_data(p_vector,designRDM,n_trials=10000,TC=TC)
# Compare likelihoods
dadmC <- EMC2:::design_model(datCT,designRDM)
dadmR <- EMC2:::design_model(datCT,designRDM1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Fitting
emc <- make_emc(datCT,designRDM,type="single")
system.time({emccCT <- fit(emc)})
# Time difference of 5.866442 mins
emc <- make_emc(datCT,designRDM1,type="single")
system.time({emcrCT <- fit(emc)})
# Time difference of 3.282171 mins


# C fails, R fine (same pattern with no contamination)
print(recovery(emccCT,p_vector,selection="alpha"))
print(recovery(emcrCT,p_vector,selection="alpha"))

# As likelihoods match at the generating value it seems likely
# that there is a failure in the measures taken to make the likelihood
# robust in the C ... will have to leave that to you Zach, sorry.

# save(dat,emcc,emcr,datC,emccC,emcrC,datT,emccT,emcrT,datCT,emccCT,emcrCT,
#      file="Censoring TEST/RDM.RData")

#### LBA, Parameter recovery, test memory consumption of old vs new C ----

load("Censoring Test/LBA.RData")

# Here is a simple 5 parameter RDM model with contamination that uses Zach's C
designLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

# Old C likelihood implementation (for memory/speed comparisons)
designOLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=OLBA,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

# This parameter vector produces reasonable RTs and accuracies with a small
# level of contamination.
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1, 1.25, log(2), log(0.2),qnorm(.05))


## This is a hack to get access to my R
LBAR <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}

# and the design that uses it.
designLBA1 <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBAR,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

TC=list(LT=0,UT=Inf,LC=0,UC=Inf,verbose=TRUE)
datCT <- make_data(p_vector,designLBA,n_trials=10000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(datCT,designLBA)
dadmO <- EMC2:::design_model(datCT,designOLBA)
lfun(p_vector,dadmC)
lfun(p_vector,dadmO)


# Fitting
emc <- make_emc(datCT,designLBA,type="single")
emccCT <- profile_fit_rss_live(fit(emc), label = "fit: LBA (C)")
# Time difference of 8.741704 mins
emc <- make_emc(datCT,designLBA1,type="single")
emcrCT <- profile_fit_rss_live(fit(emc), label = "fit: LBA (R)")
# Time difference of 2.771561 mins

# Optional: fit old C (OLBA) with the same data, to compare memory behaviour.
RUN_OLD_C_FIT <- FALSE
if (isTRUE(RUN_OLD_C_FIT)) {
  emc <- make_emc(datCT,designOLBA,type="single")
  emcoCT <- profile_fit_rss_live(fit(emc), label = "fit: OLBA (Old_C)")
}


# C fails, R fine
print(recovery(emccCT,p_vector,selection="alpha"))
print(recovery(emcrCT,p_vector,selection="alpha"))

# save(datCT,emccCT,emcrCT,file="Censoring TEST/LBA.RData")

#### LNR, Parameter recovery, all three only ----

load("Censoring Test/LNR.RData")

# Here is a simple 5 parameter RDM model with contamination that uses Zach's C
designLNR <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LNR,
  formula=list(m~lM,s~1,t0~1,pContaminant~1)
)

# This parameter vector produces reasonable RTs and accuracies with a small
# level of contamination.
p_vector <- sampled_pars(designLNR,doMap = FALSE)
p_vector[] <- c(log(1.5),log(.4), log(1), log(0.3),qnorm(.05))


## This is a hack to get access to my R
LNRR <- function() {
  m <- LNR()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}

# and the design that uses it.
designLNR1 <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LNRR,
  formula=list(m~lM,s~1,t0~1,pContaminant~1)
)

TC=list(LT=.4,UT=2.2,LC=.45,UC=1.8,verbose=TRUE)
datCT <- make_data(p_vector,designLNR,n_trials=10000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(datCT,designLNR)
dadmR <- EMC2:::design_model(datCT,designLNR1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)


# Fitting
emc <- make_emc(datCT,designLNR,type="single")
system.time({emccCT <- fit(emc)})
# Time difference of 9.157945 mins
emc <- make_emc(datCT,designLNR1,type="single")
system.time({emcrCT <- fit(emc)})
# Time difference of 2.199915 mins


# C fails, R fine
print(recovery(emccCT,p_vector,selection="alpha"))
print(recovery(emcrCT,p_vector,selection="alpha"))

# save(datCT,emccCT,emcrCT,file="Censoring TEST/LNR.RData")

dev.off()
### GNG LBA ----

load("Censoring Test/GNG.RData")

# Go/NoGo designs are detected via the use of a "nogo" level in the response.
designLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,
  formula=list(v~lM,B~1,t0~1)
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1, 1.25, log(2), log(0.2))
dat <- make_data(p_vector,designLBA,n_trials=10000)

# All responses that would have been nogo are coded as NA with rt=Inf
tail(dat)
unique(dat$R)

# Response probabilities as expected
tapply(is.na(dat$R),dat$S,mean)

LBAR <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}
designLBA1 <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("go","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBAR,
  formula=list(v~lM,B~1,t0~1)
)


# Compare likelihoods
dadmC <- EMC2:::design_model(dat,designLBA)
dadmR <- EMC2:::design_model(dat,designLBA1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Profiles
print(profile_plot(dat,designLBA,p_vector,n_cores=1,layout=c(2,2)))
print(profile_plot(dat,designLBA1,p_vector,n_cores=1,layout=c(2,2),use_c=FALSE))


# Fitting
emc <- make_emc(dat,designLBA,type="single")
system.time({emcc <- fit(emc)})
# Time difference of 26.37225 secs
emc <- make_emc(dat,designLBA1,type="single")
system.time({emcr <- fit(emc)})
# Time difference of 1.736852 mins

# Both good
print(recovery(emcc,p_vector,selection="alpha"))
print(recovery(emcr,p_vector,selection="alpha"))


### GNG LNR 3 Choice ----

# Race models allow go/nogo with more than two choices, here with an LNR
designLNR <- design(
  factors=list(subjects=1,S=c("left","right","stop")),Rlevels=c("left","right","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LNR,
  formula=list(m~lM,s~1,t0~1)
)

p_vector <- sampled_pars(designLNR,doMap = FALSE)
p_vector[] <- c(log(1.5),log(.3), log(1), log(0.3))

dat3 <- make_data(p_vector,designLNR,n_trials=10000)

# All responses that would have been nogo are coded as NA with rt=Inf
unique(dat3$R)

# Response probabilities as expected, about 70% accuracy
R <- as.character(dat3$R); R[is.na(dat3$R)] <- "nogo"
apply(table(R,S=dat3$S),2,\(x) x/10000)


LNRR <- function() {
  m <- LNR()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(pars, dadm, model, min_ll = min_ll)}
  m$c_name=NULL
  m
}
designLNR1 <- design(
  factors=list(subjects=1,S=c("left","right","stop")),Rlevels=c("left","right","nogo"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LNRR,
  formula=list(m~lM,s~1,t0~1)
)

# Compare likelihoods
dadmC <- EMC2:::design_model(dat3,designLNR)
dadmR <- EMC2:::design_model(dat3,designLNR1)
lfun(p_vector,dadmC)
lfun(p_vector,dadmR,use_c = F)

# Profiles
print(profile_plot(dat3,designLNR,p_vector,n_cores=1,layout=c(2,2)))
print(profile_plot(dat3,designLNR1,p_vector,n_cores=1,layout=c(2,2),use_c=FALSE))


# Fitting
emc <- make_emc(dat3,designLNR,type="single")
system.time({emcc3 <- fit(emc)})
# Time difference of 48.19057 secs
emc <- make_emc(dat3,designLNR1,type="single")
system.time({emcr3 <- fit(emc)})
# Time difference of 5.302002 mins

# Both good
print(recovery(emcc3,p_vector,selection="alpha"))
print(recovery(emcr3,p_vector,selection="alpha"))


### DDM GNG ----

# GNG has been implemented for the DDM in a clumsy fashion and only in R.
# Here is the current implementation which needs to be put into the new GNG
# structure and ported to C.

# For speed do a WDM example. We must specify an upper limit of the response window
TIMEOUT <- 2.5
designDDM <- design(Rlevels = c("left","right"),
              factors=list(subjects=1,S=c("left","right")),
              functions=list(
                # This is the same as UC = 2.5
                TIMEOUT=function(d)rep(TIMEOUT,nrow(d)),
                # specify nogo response level
                Rnogo=function(d)factor(rep("left",nrow(d)),levels=c("left","right")),
                # specify go response level
                Rgo=function(d)factor(rep("right",nrow(d)),levels=c("left","right"))),
                formula=list(v~S,a~1, Z~1, t0~1),
                model=DDMGNG)


p_vector <- sampled_pars(designDDM,doMap = FALSE)
p_vector[] <- c(-.5, 1, log(2), qnorm(.5),log(.3))

# Nogo responses are indicated by rt=NA, the R column is ignored for these trials
datD <- make_data(p_vector,designDDM,n_trials=10000)

# Response probabilities are as expected and there is clear truncation
tapply(is.na(datD$rt),datD$S,mean)
hist(datD$rt,breaks=seq(0,TIMEOUT,.1))

# Compare likelihoods
dadm <- EMC2:::design_model(datD,designDDM)
lfun(p_vector,dadm,use_c = F)

# Profiles
print(profile_plot(datD,designDDM,p_vector,n_cores=1,layout=c(2,3),use_c=FALSE))

# Fitting
emc <- make_emc(datD,designDDM,type="single")
system.time({emcrD <- fit(emc)})
#Time difference of 1.158099 mins

print(recovery(emcrD,p_vector,selection="alpha"))

# save(dat,emcc,emcr,dat3,emcc3,emcr3,datD,emcrD,file="Censoring TEST/GNG.RData")
