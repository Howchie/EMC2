pkgname <- "EMC2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('EMC2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BAwD")
### * BAwD

flush(stderr()); flush(stdout())

### Name: BAwD
### Title: The Ballistic Accumulator with Drive Decay (BAwD)
### Aliases: BAwD

### ** Examples

ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
matchfun <- function(d) d$S == d$lR
# ell is left out of the formula, so it stays at its default of 1 and fixes
# the evidence scale for the lognormal launch.
design_BAwD <- design(data = forstmann, model = BAwD, matchfun = matchfun,
                      formula = list(mu ~ lM, sigma ~ 1, B ~ E, A ~ 1,
                                     t0 ~ 1, k ~ 1),
                      contrasts = list(mu = list(lM = ADmat)))



cleanEx()
nameEx("BOU")
### * BOU

flush(stderr()); flush(stdout())

### Name: BOU
### Title: The Bounded Ornstein-Uhlenbeck Model
### Aliases: BOU

### ** Examples

design_BOU <- design(data = forstmann, model = BOU,
                     formula = list(v ~ 0 + S, a ~ E, beta ~ 1, t0 ~ 1,
                                    s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1),
                     constants = c(s = log(1)))

# Leak and a symmetrically collapsing boundary together
design_BOUc <- design(data = forstmann,
                      model = function() BOU(boundary_collapse = "exponential"),
                      formula = list(v ~ 0 + S, a ~ E, beta ~ 1, t0 ~ 1,
                                     s ~ 1, Z ~ 1, aInf ~ 1, tau ~ 1),
                      constants = c(s = log(1)))



cleanEx()
nameEx("DDM")
### * DDM

flush(stderr()); flush(stdout())

### Name: DDM
### Title: The Diffusion Decision Model
### Aliases: DDM

### ** Examples

design_DDMaE <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).




cleanEx()
nameEx("DDMGNG")
### * DDMGNG

flush(stderr()); flush(stdout())

### Name: DDMGNG
### Title: The GNG (go/nogo) Diffusion Decision Model
### Aliases: DDMGNG

### ** Examples

dGNG <- design(Rlevels = c("left","right"),
               factors=list(subjects=1,S=c("left","right")),
               functions=list(
               TIMEOUT=function(d)rep(2.5,nrow(d)),
               # no go response level
               Rnogo=function(d)factor(rep("left",nrow(d)),levels=c("left","right")),
               # go response level
               Rgo=function(d)factor(rep("right",nrow(d)),levels=c("left","right"))),
               formula=list(v~S,a~1, Z~1, t0~1),
               model=DDMGNG)

p_vector <- sampled_pars(dGNG)



cleanEx()
nameEx("LBA")
### * LBA

flush(stderr()); flush(stdout())

### Name: LBA
### Title: The Linear Ballistic Accumulator model
### Aliases: LBA

### ** Examples

# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
# We now construct our design, with v ~ lM and the contrast for lM the ADmat.
design_LBABE <- design(data = forstmann,model=LBA,matchfun=matchfun,
                       formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                       contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).



cleanEx()
nameEx("LNR")
### * LNR

flush(stderr()); flush(stdout())

### Name: LNR
### Title: The Log-Normal Race Model
### Aliases: LNR

### ** Examples

# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
# We now construct our design, with v ~ lM and the contrast for lM the ADmat.
design_LNRmE <- design(data = forstmann,model=LNR,matchfun=matchfun,
                       formula=list(m~lM + E,s~1,t0~1),
                       contrasts=list(m=list(lM=ADmat)))
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).



cleanEx()
nameEx("MRI")
### * MRI

flush(stderr()); flush(stdout())

### Name: MRI
### Title: GLM model for fMRI data
### Aliases: MRI

### ** Examples

# Create a normal MRI model specification
model_spec <- MRI()

# Access model parameters
model_spec$p_types



cleanEx()
nameEx("MRI_AR1")
### * MRI_AR1

flush(stderr()); flush(stdout())

### Name: MRI_AR1
### Title: Create an AR(1) GLM model for fMRI data
### Aliases: MRI_AR1

### ** Examples

# Create an AR(1) GLM model for fMRI data
model_spec <- MRI_AR1()

# Access model parameters
model_spec$p_types



cleanEx()
nameEx("RDM")
### * RDM

flush(stderr()); flush(stdout())

### Name: RDM
### Title: The Racing Diffusion Model
### Aliases: RDM

### ** Examples

# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
# We now construct our design, with v ~ lM and the contrast for lM the ADmat.
design_RDMBE <- design(data = forstmann,model=RDM,matchfun=matchfun,
                       formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
                       contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).



cleanEx()
nameEx("ROU")
### * ROU

flush(stderr()); flush(stdout())

### Name: ROU
### Title: The Racing Ornstein-Uhlenbeck Model (Leaky Accumulator)
### Aliases: ROU

### ** Examples

# A leaky accumulator design with the same structure as the RDM example.
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
design_ROU <- design(data = forstmann,model=ROU,matchfun=matchfun,
                     formula=list(v~lM,k~1,B~E+lR,A~1,t0~1),
                     contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))

# The same race with an exponentially collapsing threshold.
design_ROUc <- design(data = forstmann,
                      model=function() ROU(boundary_collapse="exponential"),
                      matchfun=matchfun,
                      formula=list(v~lM,k~1,B~E+lR,A~1,t0~1,Binf~1,tau~1),
                      contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))

# The curvature parameterization: stimulus and speed emphasis move the
# reference crossing time, while physical leak and diffusion can be shared.
design_ROUcurv <- design(data = forstmann,
                         model=function() ROU(parameterization="curvature"),
                         matchfun=matchfun,
                         formula=list(tstar~lM+E,k~1,s~1,A~1,t0~1),
                         contrasts=list(tstar=list(lM=ADmat)),
                         constants=c(B=log(1)))

# The equilibrium parameterization, for a deadline design in which theta is free
# to sit below the bound. Give the data a UC column to censor at it.
design_ROUeq <- design(data = forstmann,
                       model=function() ROU(parameterization="equilibrium"),
                       matchfun=matchfun,
                       formula=list(tk~1,theta~lM,chi~1,A~1,t0~1),
                       contrasts=list(theta=list(lM=ADmat)),
                       constants=c(B=log(1)))



cleanEx()
nameEx("SDT")
### * SDT

flush(stderr()); flush(stdout())

### Name: SDT
### Title: Gaussian Signal Detection Theory Model for Binary Responses
### Aliases: SDT

### ** Examples

dprobit <- design(Rlevels = c("left","right"),
           factors=list(subjects=1,S=c("left","right")),
           formula=list(mean ~ 0+S, sd ~ 1,threshold ~ 1),
           matchfun=function(d)d$S==d$lR,
           constants=c(sd=log(1),threshold=0),
           model=SDT)

p_vector <- sampled_pars(dprobit)



cleanEx()
nameEx("align_loadings")
### * align_loadings

flush(stderr()); flush(stdout())

### Name: align_loadings
### Title: Reorder MCMC Samples of Factor Loadings
### Aliases: align_loadings

### ** Examples

# This function works natively with emc objects, but also factor arrays:
# Simulate a small example with 5 variables, 2 factors, and 10 MCMC iterations
set.seed(123)
p <- 5  # Number of variables
q <- 2  # Number of factors
n <- 10 # Number of MCMC iterations

# Create random factor loadings with label switching
lambda <- array(0, dim = c(p, q, n))
for (i in 1:n) {
  # Generate base loadings
  base_loadings <- matrix(rnorm(p*q, 0, 0.5), p, q)
  base_loadings[1:3, 1] <- abs(base_loadings[1:3, 1]) + 0.5  # Strong loadings on factor 1
  base_loadings[4:5, 2] <- abs(base_loadings[4:5, 2]) + 0.5  # Strong loadings on factor 2

  # Randomly switch labels and signs
  if (runif(1) > 0.5) {
    # Switch factor order
    base_loadings <- base_loadings[, c(2, 1)]
  }
  if (runif(1) > 0.5) {
    # Switch sign of factor 1
    base_loadings[, 1] <- -base_loadings[, 1]
  }
  if (runif(1) > 0.5) {
    # Switch sign of factor 2
    base_loadings[, 2] <- -base_loadings[, 2]
  }

  lambda[,,i] <- base_loadings
}

# Align the loadings
result <- align_loadings(lambda = lambda, verbose = TRUE, n_cores = 1)

# Examine the aligned loadings
print(result)




cleanEx()
nameEx("chain_n")
### * chain_n

flush(stderr()); flush(stdout())

### Name: chain_n
### Title: MCMC Chain Iterations
### Aliases: chain_n

### ** Examples

chain_n(samples_LNR)




cleanEx()
nameEx("check")
### * check

flush(stderr()); flush(stdout())

### Name: check.emc
### Title: Convergence Checks for an emc Object
### Aliases: check.emc check

### ** Examples

check(samples_LNR)



cleanEx()
nameEx("compare")
### * compare

flush(stderr()); flush(stdout())

### Name: compare
### Title: Information Criteria and Marginal Likelihoods
### Aliases: compare

### ** Examples




cleanEx()
nameEx("compare_subject")
### * compare_subject

flush(stderr()); flush(stdout())

### Name: compare_subject
### Title: Information Criteria For Each Participant
### Aliases: compare_subject

### ** Examples

# For a broader illustration see `compare`.
# Here we just take two times the same model, but normally one would compare
# different models
compare_subject(list(m0 = samples_LNR, m1 = samples_LNR))



cleanEx()
nameEx("contr.anova")
### * contr.anova

flush(stderr()); flush(stdout())

### Name: contr.anova
### Title: Anova Style Contrast Matrix
### Aliases: contr.anova

### ** Examples

{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.anova),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}



cleanEx()
nameEx("contr.bayes")
### * contr.bayes

flush(stderr()); flush(stdout())

### Name: contr.bayes
### Title: Contrast Enforcing Equal Prior Variance on each Level
### Aliases: contr.bayes

### ** Examples

{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.bayes),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}



cleanEx()
nameEx("contr.decreasing")
### * contr.decreasing

flush(stderr()); flush(stdout())

### Name: contr.decreasing
### Title: Contrast Enforcing Decreasing Estimates
### Aliases: contr.decreasing

### ** Examples

{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.decreasing),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}



cleanEx()
nameEx("contr.increasing")
### * contr.increasing

flush(stderr()); flush(stdout())

### Name: contr.increasing
### Title: Contrast Enforcing Increasing Estimates
### Aliases: contr.increasing

### ** Examples

{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.increasing),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}



cleanEx()
nameEx("convolve_design_matrix")
### * convolve_design_matrix

flush(stderr()); flush(stdout())

### Name: convolve_design_matrix
### Title: Convolve Events with HRF to Construct Design Matrices
### Aliases: convolve_design_matrix

### ** Examples

# Generate a simple example timeseries
ts <- data.frame(
  subjects = rep(1, 100),
  run = rep(1, 100),
  time = seq(0, 99),
  ROI1 = rnorm(100)
)

# Generate example events
events <- data.frame(
  subjects = rep(1, 4),
  run = rep(1, 4),
  onset = c(10, 30, 50, 70),
  duration = rep(0.5, 4),
  event_type = c("hard", "easy", "hard", "easy"),
  modulation = c(1, 1, 1, 1)
)

# Build design matrices
design_matrices <-  convolve_design_matrix(
  timeseries = ts,
  events = events,
  factors = list(difficulty = c("hard", "easy")),
  contrasts = list(difficulty = matrix(c(-1, 1)))
)



cleanEx()
nameEx("credible")
### * credible

flush(stderr()); flush(stdout())

### Name: credible.emc
### Title: Posterior Credible Interval Tests
### Aliases: credible.emc credible

### ** Examples

{
# Run a credible interval test (Bayesian ''t-test'')
credible(samples_LNR, x_name = "m")
# We can also compare between two sets of emc objects

# # Now without a ~ E
# design_null <- design(data = forstmann,model=DDM,
#                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#                            constants=c(s=log(1)))
#
# null_model <- make_emc(forstmann, design_null)
# null_model <- fit(null_model)
# credible(x = null_model, x_name = "a", y = full_model, y_name = "a")
#
# # Or provide custom functions:
# credible(x = full_model, x_fun = function(d) d["a_Eaccuracy"] - d["a_Eneutral"])
}



cleanEx()
nameEx("credint")
### * credint

flush(stderr()); flush(stdout())

### Name: credint.emc.prior
### Title: Posterior Quantiles
### Aliases: credint.emc.prior credint.emc credint

### ** Examples

credint(samples_LNR)



cleanEx()
nameEx("design")
### * design

flush(stderr()); flush(stdout())

### Name: design
### Title: Specify a Design and Model
### Aliases: design

### ** Examples


# load example dataset
dat <- forstmann

# create a function that takes the latent response (lR) factor (d) and returns a logical
# defining the correct response for each stimulus. Here the match is simply
# such that the S factor equals the latent response factor
matchfun <- function(d)d$S==d$lR

# When working with lM and lR, it can be useful to design  an
# "average and difference" contrast matrix. For binary responses, it has a
# simple canonical form
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# Create a design for a linear ballistic accumulator model (LBA) that allows
# thresholds to be a function of E and lR. The final result is a 9 parameter model.
design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
                            formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),
                            constants=c(sv=log(1)))



cleanEx()
nameEx("design_fmri")
### * design_fmri

flush(stderr()); flush(stdout())

### Name: design_fmri
### Title: Create fMRI Design for EMC2 Sampling
### Aliases: design_fmri

### ** Examples

# Generate a simple example timeseries
ts <- data.frame(
 subjects = rep(1, 100),
 run = rep(1, 100),
 time = cumsum(rep(1.38, 100)),
 ROI1 = rnorm(100)
)

# Generate example events
events <- data.frame(
 subjects = rep(1, 4),
 run = rep(1, 4),
 onset = c(10, 30, 50, 70),
 duration = rep(0.5, 4),
 event_type = c("A", "B", "A", "B"),
 modulation = c(1, 1, 1, 1)
)

# Create convolved design matrix
design_matrix <- convolve_design_matrix(
 timeseries = ts,
 events = events,
 factors = list(condition = c("A", "B")),
 hrf_model = "glover"
)

# Create fMRI design for EMC2
fmri_design <- design_fmri(design_matrix, model = MRI_AR1)



cleanEx()
nameEx("ess_summary")
### * ess_summary

flush(stderr()); flush(stdout())

### Name: ess_summary.emc
### Title: Effective Sample Size
### Aliases: ess_summary.emc ess_summary

### ** Examples

ess_summary(samples_LNR, selection = "alpha")



cleanEx()
nameEx("fit")
### * fit

flush(stderr()); flush(stdout())

### Name: fit.emc
### Title: Model Estimation in EMC2
### Aliases: fit.emc fit

### ** Examples




cleanEx()
nameEx("gd_summary")
### * gd_summary

flush(stderr()); flush(stdout())

### Name: gd_summary.emc
### Title: Gelman-Rubin Statistic
### Aliases: gd_summary.emc gd_summary

### ** Examples

gd_summary(samples_LNR, selection = "correlation", stat = "mean", flatten = TRUE)



cleanEx()
nameEx("get_BayesFactor")
### * get_BayesFactor

flush(stderr()); flush(stdout())

### Name: get_BayesFactor
### Title: Bayes Factors
### Aliases: get_BayesFactor

### ** Examples




cleanEx()
nameEx("get_data")
### * get_data

flush(stderr()); flush(stdout())

### Name: get_data.emc
### Title: Get Data
### Aliases: get_data.emc get_data

### ** Examples

get_data(samples_LNR)



cleanEx()
nameEx("get_design")
### * get_design

flush(stderr()); flush(stdout())

### Name: get_design.emc.prior
### Title: Get Design
### Aliases: get_design.emc.prior get_design.emc get_design

### ** Examples

get_design(samples_LNR)



cleanEx()
nameEx("get_pars")
### * get_pars

flush(stderr()); flush(stdout())

### Name: get_pars
### Title: Filter/Manipulate Parameters from emc Object
### Aliases: get_pars

### ** Examples

# E.g. get the group-level mean parameters mapped back to the design
get_pars(samples_LNR, stage = "sample", map = TRUE, selection = "mu")

# Or return the flattened correlation, with 10 iterations per chain
get_pars(samples_LNR, stage = "sample", selection = "correlation", flatten = TRUE, length.out = 10)



cleanEx()
nameEx("get_prior")
### * get_prior

flush(stderr()); flush(stdout())

### Name: get_prior.emc
### Title: Get Prior
### Aliases: get_prior.emc get_prior

### ** Examples

get_prior(samples_LNR)



cleanEx()
nameEx("get_trend_pnames")
### * get_trend_pnames

flush(stderr()); flush(stdout())

### Name: get_trend_pnames
### Title: Get parameter types from trend object
### Aliases: get_trend_pnames

### ** Examples

trend <- make_trend(par_names = "v", cov_names = "trial", kernels = "exp_incr")
get_trend_pnames(trend)




cleanEx()
nameEx("group_design")
### * group_design

flush(stderr()); flush(stdout())

### Name: group_design
### Title: Create Group-Level Design Matrices
### Aliases: group_design

### ** Examples

# Create subject-level design
subj_design <- design(data = forstmann, model = DDM,
                      formula = list(v ~ S, a ~ E, t0 ~ 1),
                      contrasts = list(S = contr.helmert))
# Add some age covariate and roughly demeans
# Demeaning is important to ensure that the interpretation of the group-level intercept
# is the mean of the group (i.e., 'mu' still represents the group-level mean)
forstmann$age <- as.numeric(forstmann$subjects) -mean(as.numeric(forstmann$subjects))
# Create fake group column
forstmann$group <- ifelse(forstmann$subjects %in%
              unique(forstmann$subjects)[seq(1, 19, 2)], "A", "B")

# Create group-level design matrices
group_des <- group_design(
  formula = list(v_S1 ~ age + group, a ~ age),
  data = forstmann,
  subject_design = subj_design,
  contrasts = list(group = contr.bayes)
)
# Then you can make the emc object with
emc <- make_emc(forstmann, subj_design, compress = FALSE, group_design = group_des)



cleanEx()
nameEx("high_pass_filter")
### * high_pass_filter

flush(stderr()); flush(stdout())

### Name: high_pass_filter
### Title: Apply High-Pass Filtering to fMRI Data
### Aliases: high_pass_filter

### ** Examples

# Create a simple example data frame with drift
set.seed(123)
n_frames <- 100
time <- seq(0, 99)

# Create a signal with low-frequency drift
drift <- 0.1 * time
signal <- sin(2 * pi * 0.1 * time) + drift
noise <- rnorm(n_frames, 0, 0.5)
data <- signal + noise

# Create a data frame
df <- data.frame(
  time = time,
  signal = data
)

# Apply high-pass filtering using cosine basis functions
filtered_df <- high_pass_filter(df, high_pass_model = "cosine")



cleanEx()
nameEx("hypothesis")
### * hypothesis

flush(stderr()); flush(stdout())

### Name: hypothesis.emc
### Title: Within-Model Hypothesis Testing
### Aliases: hypothesis.emc hypothesis

### ** Examples

# Here the emc object has an effect parameter (e.g. m),
# that maps onto a certain hypothesis.
# The hypothesis here is that m is different from zero.
# We can test whether there's a group-level effect on m:
hypothesis(samples_LNR, parameter = "m")
# Alternatively we can also test whether two parameters differ from each other
mdiff <- function(p)diff(p[c("m","m_lMd")])
hypothesis(samples_LNR,fun=mdiff)



cleanEx()
nameEx("init_chains")
### * init_chains

flush(stderr()); flush(stdout())

### Name: init_chains
### Title: Initialize Chains
### Aliases: init_chains

### ** Examples




cleanEx()
nameEx("make_data")
### * make_data

flush(stderr()); flush(stdout())

### Name: make_data
### Title: Simulate Data
### Aliases: make_data

### ** Examples

# First create a design
design_DDMaE <- design(factors = list(S = c("left", "right"),
                                           E = c("SPD", "ACC"),
                                           subjects = 1:30),
                            Rlevels = c("left", "right"), model = DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
# Then create a p_vector:
parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_EACC=log(2), t0=log(.2),
              Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))

# Now we can simulate data
data <- make_data(parameters, design_DDMaE, n_trials = 30)

# We can also simulate data based on a specific dataset
design_DDMaE <- design(data = forstmann,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
              t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))

data <- make_data(parameters, design_DDMaE, data = forstmann)



cleanEx()
nameEx("make_emc")
### * make_emc

flush(stderr()); flush(stdout())

### Name: make_emc
### Title: Make an emc Object
### Aliases: make_emc

### ** Examples

dat <- forstmann

# function that takes the lR factor (named diff in the following function) and
# returns a logical defining the correct response for each stimulus. In this
# case the match is simply such that the S factor equals the latent response factor.
matchfun <- function(d)d$S==d$lR

# design an "average and difference" contrast matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# specify design
design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))

# specify priors
pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
           B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
         B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)

# create emc object
LBABE <- make_emc(dat,design_LBABE,type="standard",  prior=prior_LBABE,
                  compress = FALSE)




cleanEx()
nameEx("make_missing")
### * make_missing

flush(stderr()); flush(stdout())

### Name: make_missing
### Title: Add information about missing values to data and modify/filter
###   accordingly.
### Aliases: make_missing

### ** Examples

## Not run: 
##D # First make some data
##D   designRDM <- design(model = RDM,
##D   factors = list(subjects = 1:2, S = c("left", "right")),Rlevels = c("left", "right"),
##D   matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
##D   formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
##D   constants = c(s = log(1)))
##D p_vector <- log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8))
##D dat <- make_data(p_vector, designRDM,n_trials = 10)
##D 
##D # Filter data frame without LT/UC/LT/UT columns (as in most real data files)
##D data <- dat
##D mdata <- make_missing(dat,LT=.7,LC=.75,UC=1.5,UT=1.6,verbose=TRUE)
## End(Not run)



cleanEx()
nameEx("make_random_effects")
### * make_random_effects

flush(stderr()); flush(stdout())

### Name: make_random_effects
### Title: Generate Subject-Level Parameters
### Aliases: make_random_effects

### ** Examples

# First create a design
design_DDMaE <- design(data = forstmann,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
# Then create a group-level means vector:
group_means =c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
               t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
# Now we can create subject-level parameters
subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19)

# We can also define a covariance matrix to simulate from
subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19,
             covariances = diag(.1, length(group_means)))

# The subject level parameters can be used to generate data
make_data(subj_pars, design_DDMaE, n_trials = 10)



cleanEx()
nameEx("make_sem_structure")
### * make_sem_structure

flush(stderr()); flush(stdout())

### Name: make_sem_structure
### Title: Define Structural Equation Model (SEM) Matrices
### Aliases: make_sem_structure

### ** Examples

# Create a design object (simplified from design.R example)
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
matchfun_example <- function(d) d$S==d$lR # Example match function

example_design_obj <- design(
  data = forstmann,
  model= LBA,
  matchfun=matchfun_example,
  formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
  contrasts=list(v=list(lM=ADmat)),
  constants=c(sv=log(1)),
)

# SEM Factor names

# Make a copy of forstmann for example modification
forstmann_mod <- forstmann
set.seed(123) # for reproducibility
subj_trait_values <- stats::setNames(rnorm(length(levels(forstmann_mod$subjects))),
                                    levels(forstmann_mod$subjects))
forstmann_mod$SubjTrait <- subj_trait_values[forstmann_mod$subjects]

my_cov_cols <- c("SubjTrait")

lambda_example_specs <- list(
  Speed = c("v", "v_lMdiff"), # "v" will be fixed to 1
  Caution = c("B", "B_Eneutral", "B_Eaccuracy", "B_lRright", "A") # "B" fixed to 1
)
b_example_specs <- list(Caution = "Speed")
k_example_specs <- list(t0 = "SubjTrait") # "SubjTrait" must be in my_cov_cols
g_example_specs <- list(Speed = "SubjTrait")

sem_settings_definition <- make_sem_structure(
  data = forstmann_mod,
  design = example_design_obj,
  covariate_cols = my_cov_cols,
  lambda_specs = lambda_example_specs,
  b_specs = b_example_specs,
  k_specs = k_example_specs,
  g_specs = g_example_specs
)

print(sem_settings_definition$Lambda_mat)
print(sem_settings_definition$B_mat)
print(sem_settings_definition$K_mat)
print(sem_settings_definition$G_mat)
print(head(sem_settings_definition$covariates))




cleanEx()
nameEx("make_ssd")
### * make_ssd

flush(stderr()); flush(stdout())

### Name: make_ssd
### Title: Create a stop-signal delay generator
### Aliases: make_ssd

### ** Examples

# Fixed SSDs sampled on 25% of trials
ssd_fixed <- make_ssd(values = c(.26, .35, .46), p = rep(.25 / 3, 3), staircase = FALSE)

# Staircase with default parameters operating per subject and stimulus
ssd_stair <- make_ssd(factors = c("subjects", "S"))



cleanEx()
nameEx("make_trend")
### * make_trend

flush(stderr()); flush(stdout())

### Name: make_trend
### Title: Create a trend specification for model parameters
### Aliases: make_trend

### ** Examples

# Put trend on B and v parameters
trend <- make_trend(
  par_names = c("B", "v"),
  cov_names = "strial",
  kernels = c("exp_incr", "poly3"),
  phase = "premap",
  shared = list(shrd = list("B.B0", "v.d1"))
)
get_trend_pnames(trend)


# Using covariate maps

# Covariate maps allow you to specify how trial-by-trial covariates influence
# model parameters for each accumulator. The example below uses a simple data
# frame with two trials. `cov_left` and `cov_right` specify which covariates
# correspond to the left and right accumulators on each trial. `S` indicates the
# correct response, and `cov1`–`cov4` contain the actual covariate values.

data <- data.frame(
  subjects = rep(1, 2),
  S        = c('left', 'right'),
  cov_left = c('1', '4'),
  cov_right= c('3', '2'),
  rt       = c(1.2, 0.8),
  R        = factor(c('left', 'right')),
  cov1     = c(1, NA),
  cov2     = c(NA, 1),
  cov3     = c(NA, 1),
  cov4     = c(1, 1)
)

# A covariate map function must take `dadm` and `cov_names` as inputs and return
# a matrix of size (nrow(dadm), length(cov_names)), coding how each covariate
# contributes to each accumulator.

advantage_map <- function(dadm, cov_names) {

  # Which stimulus does the accumulator correspond to on each trial?
  lS <- paste0('cov', ifelse(dadm$lR == 'left', dadm$cov_left, dadm$cov_right))

  # Which stimulus does the *other* accumulator correspond to?
  lSother <- paste0('cov', ifelse(dadm$lR == 'right', dadm$cov_left, dadm$cov_right))

  # Build indicator matrices
  map_plus1 <- sapply(cov_names, function(col) ifelse(lS     == col,  1, 0))
  map_minus1<- sapply(cov_names, function(col) ifelse(lSother == col, -1, 0))

  map_plus1 + map_minus1
}

# A covariate map function can be supplied to make_trend(), which creates the mapping
# specification for the model for each participant. Here, a single map ('differences') is provided.

trend <- make_trend(
  par_names = 'v',
  kernels   = 'delta',
  bases     = 'lin',
  cov_names = list(c('cov1', 'cov2', 'cov3', 'cov4')),
  maps      = list('differences' = advantage_map),
  at        = 'lR'
)

design_RDM <- design(
  model  = RDM,
  data   = data,
  formula= list(B ~ 1, v ~ 1, t0 ~ 1),
  trend  = trend
)

emc <- make_emc(data, design_RDM, type = 'single')

# The resulting covariate maps for each subject are attached to the `dadm`:
attr(emc[[1]]$data[[1]], 'covariate_maps')
# And to confirm that this mapping is correct, compare with the corresponding `dadm`
emc[[1]]$data[[1]]

# You can also provide multiple covariate maps. Each additional map introduces
# a separate base parameter. For example, the following `sum_map` is suitable
# for RL-ARD–type models:

sum_map <- function(dadm, cov_names) {
  # Which stimulus does the accumulator correspond to on each trial?
  lS <- paste0('cov', ifelse(dadm$lR == 'left', dadm$cov_left, dadm$cov_right))
  # Which stimulus does the *other* accumulator correspond to?
  lSother <- paste0('cov', ifelse(dadm$lR == 'right', dadm$cov_left, dadm$cov_right))

  # Indicator matrices (note: both are added rather than subtracted)
  map_this  <- sapply(cov_names, function(col) ifelse(lS      == col, 1, 0))
  map_other <- sapply(cov_names, function(col) ifelse(lSother == col, 1, 0))

  map_this + map_other
}

trend <- make_trend(
  par_names = 'v',
  kernels   = 'delta',
  bases     = 'lin',
  cov_names = list(c('cov1', 'cov2', 'cov3', 'cov4')),
  maps      = list('differences' = advantage_map,
                   'sums'        = sum_map),
  at        = 'lR'
)

design_RDM <- design(
  model  = RDM,
  data   = data,
  formula= list(B ~ 1, v ~ 1, t0 ~ 1),
  trend  = trend
)

emc <- make_emc(data, design_RDM, type = 'single')

# Now the dadm contains two covariate maps, and the model includes two
# corresponding base parameters (e.g., v.w1 and v.w2):
attr(emc[[1]]$data[[1]], 'covariate_maps')





cleanEx()
nameEx("mapped_pars")
### * mapped_pars

flush(stderr()); flush(stdout())

### Name: mapped_pars
### Title: Parameter Mapping Back to the Design Factors
### Aliases: mapped_pars mapped_pars.emc.design mapped_pars.emc.prior
###   mapped_pars.emc

### ** Examples

# First define a design:
design_DDMaE <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
mapped_pars(design_DDMaE)
# Then create a p_vector:
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
          t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
# This will map the parameters of the p_vector back to the design
mapped_pars(design_DDMaE, p_vector)




cleanEx()
nameEx("model_averaging")
### * model_averaging

flush(stderr()); flush(stdout())

### Name: model_averaging
### Title: Model Averaging
### Aliases: model_averaging

### ** Examples

# First set up some example models (normally these would be alternative models)
samples_LNR2 <- subset(samples_LNR, length.out = 45)
samples_LNR3 <- subset(samples_LNR, length.out = 40)
samples_LNR4 <- subset(samples_LNR, length.out = 35)

# Run compare on them, BayesFactor = F is set for speed.
ICs <- compare(list(S1 = samples_LNR, S2 = samples_LNR2,
                    S3 = samples_LNR3, S4 = samples_LNR4), BayesFactor = FALSE)

# Model averaging can either be done with a vector of ICs:
model_averaging(ICs$BPIC[1:2], ICs$BPIC[2:4])

# Or the output of compare:
model_averaging(ICs[1:2,], ICs[3:4,])




cleanEx()
nameEx("pairs_posterior")
### * pairs_posterior

flush(stderr()); flush(stdout())

### Name: pairs_posterior
### Title: Plot Within-Chain Correlations
### Aliases: pairs_posterior

### ** Examples




cleanEx()
nameEx("parameters")
### * parameters

flush(stderr()); flush(stdout())

### Name: parameters.emc.prior
### Title: Return Data Frame of Parameters
### Aliases: parameters.emc.prior parameters.emc parameters

### ** Examples

# For prior inference:
# First set up a prior
design_DDMaE <- design(data = forstmann,model=DDM,
                       formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                       constants=c(s=log(1)))
# Then set up a prior using make_prior
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
         t0=.4,Z=1,sv=.4,SZ=1)
# Here we left the variance prior at default
prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd)
# Get our prior samples
parameters(prior_DDMaE, N = 100)
# For posterior inference:
# Get 100 samples of the group-level mean (the default)
parameters(samples_LNR, N = 100)
# or from the individual-level parameters and mapped
parameters(samples_LNR, selection = "alpha", map = TRUE)



cleanEx()
nameEx("plot.emc")
### * plot.emc

flush(stderr()); flush(stdout())

### Name: plot.emc
### Title: Plot Function for emc Objects
### Aliases: plot.emc

### ** Examples

plot(samples_LNR)
# Or trace autocorrelation for the second subject:
plot(samples_LNR, subject = 2, selection = "alpha")

# Can also plot the trace of for example the group-level correlation:
plot(samples_LNR, selection = "correlation", col = c("green", "purple", "orange"), lwd = 2)



cleanEx()
nameEx("plot.emc.prior")
### * plot.emc.prior

flush(stderr()); flush(stdout())

### Name: plot.emc.prior
### Title: Plot a prior
### Aliases: plot.emc.prior

### ** Examples




cleanEx()
nameEx("plot_caf")
### * plot_caf

flush(stderr()); flush(stdout())

### Name: plot_caf
### Title: Plot conditional accuracy functions
### Aliases: plot_caf

### ** Examples

# Plot conditional accuracy function for data only,
# NB: the caf_factor must have two levels levels.
# forstmann_speed_accuracy <- forstmann[forstmann$E!="neutral",]
# forstmann_speed_accuracy$E <- droplevels(forstmann_speed_accuracy$E)
# plot_caf(forstmann_speed_accuracy, caf_factor="E",factors="S", smooth_window=10)
#
# Or a list of multiple emc objects ...



cleanEx()
nameEx("plot_cdf")
### * plot_cdf

flush(stderr()); flush(stdout())

### Name: plot_cdf
### Title: Plot Defective Cumulative Distribution Functions
### Aliases: plot_cdf

### ** Examples

# Plot defective CDF for data only
# plot_cdf(forstmann, to_plot = "data")
#
# Plot with posterior predictions
# plot_cdf(samples_LNR, to_plot = c("data","posterior"), n_post=10)
#
# Or a list of multiple emc objects ...



cleanEx()
nameEx("plot_delta")
### * plot_delta

flush(stderr()); flush(stdout())

### Name: plot_delta
### Title: Plot Difference of Cumulative Distribution Functions
### Aliases: plot_delta

### ** Examples

# Plot delta function for data only, not that the delta_factor must have two
# levels.
# fortsmann_speed_accuracy <- forstmann[forstmann$E!="neutral",]
# fortsmann_speed_accuracy$E <- droplevels(fortsmann_speed_accuracy$E)
# plot_delta(fortsmann_speed_accuracy, to_plot = "data")
#
# Plot with posterior predictions
# plot_delta(samples_LNR, to_plot = c("data","posterior"), n_post=10)
#
# Or a list of multiple emc objects ...



cleanEx()
nameEx("plot_density")
### * plot_density

flush(stderr()); flush(stdout())

### Name: plot_density
### Title: Plot Defective Densities
### Aliases: plot_density

### ** Examples

# Plot defective densities for each subject and the factor combination in the design:
plot_density(forstmann)
# or for one subject:
plot_density(forstmann, subject = 1)
# Now collapsing across subjects and using a different defective factor:
plot_density(forstmann, factors = "S", defective_factor = "E")
# Or plot posterior predictives
plot_density(samples_LNR, n_post = 10)



cleanEx()
nameEx("plot_design_fmri")
### * plot_design_fmri

flush(stderr()); flush(stdout())

### Name: plot_design_fmri
### Title: Plot fMRI Design Matrix
### Aliases: plot_design_fmri

### ** Examples

# Example time series
ts <- data.frame(
  subjects = rep(1, 100),
  run      = rep(1, 100),
  time     = seq(0, 99),
  ROI      = rnorm(100)
)
# Create a simple events data frame
events <- data.frame(
  subjects = rep(1, 10),
  run = rep(1, 10),
  onset = seq(0, 90, by = 10),
  condition = rep(c("A", "B"), 5),
  rt = runif(10, 0.5, 1.5),
  accuracy = sample(0:1, 10, replace = TRUE)
)
# Reshape with custom duration for each event_type
reshaped <- reshape_events(events,
                           event_types = c("condition", "accuracy", "rt"),
                           duration = list(condition = 0.5,
                                          accuracy = 0.2,
                                          rt = function(x) x$rt))
design_matrices <-  convolve_design_matrix(
                    timeseries = ts,
                    events = reshaped,
                    covariates = c('accuracy', 'rt'),
                    factors = list(cond = c("condition_A", "condition_B")),
                    contrasts = list(cond = matrix(c(-1, 1))))

# Plot the design matrix
plot_design_fmri(design_matrices)



cleanEx()
nameEx("plot_fmri")
### * plot_fmri

flush(stderr()); flush(stdout())

### Name: plot_fmri
### Title: Plot fMRI peri-stimulus time courses
### Aliases: plot_fmri

### ** Examples

ts <- data.frame(
  subjects = rep(1, 100),
  run      = rep(1, 100),
  time     = seq(0, 99),
  ROI      = rnorm(100)
)
events <- data.frame(
  subjects   = rep(1, 5),
  run        = rep(1, 5),
  onset      = c(10, 30, 50, 70, 90),
  event_type = rep("A", 5),
  modulation = rep(1, 5),
  duration   = rep(0.5, 5)
)
plot_fmri(ts, events = events, event_type = "A")



cleanEx()
nameEx("plot_pars")
### * plot_pars

flush(stderr()); flush(stdout())

### Name: plot_pars
### Title: Plots Density for Parameters
### Aliases: plot_pars

### ** Examples

# Full range of possibilities described in get_pars
plot_pars(samples_LNR)
# Or plot all subjects
plot_pars(samples_LNR, all_subjects = TRUE, col = 'purple')
# Or plot recovery
true_emc <- samples_LNR # This would normally be the data-generating samples
plot_pars(samples_LNR, true_pars = true_emc, true_args = list(col = 'blue'), adjust = 2)



cleanEx()
nameEx("plot_relations")
### * plot_relations

flush(stderr()); flush(stdout())

### Name: plot_relations
### Title: Plot Group-Level Relations
### Aliases: plot_relations

### ** Examples

# For a given set of hierarchical model samples we can make a
# correlation matrix plot.
plot_relations(samples_LNR, only_cred = TRUE, plot_cred = TRUE)
# We can also only plot the correlations where the credible interval does not include zero
plot_relations(samples_LNR, plot_means = TRUE, only_cred = TRUE)




cleanEx()
nameEx("plot_stat")
### * plot_stat

flush(stderr()); flush(stdout())

### Name: plot_stat
### Title: Plot Statistics on Data
### Aliases: plot_stat

### ** Examples

# For example plot the observed and predicted response accuracy
# Can also apply more sophisticated statistics
drt <- function(data) diff(tapply(data$rt,data[,c("E")],mean))
plot_stat(samples_LNR, stat_fun = drt, n_post = 10, stat_name = "RT diff Speed - A/N")




cleanEx()
nameEx("plot_trend")
### * plot_trend

flush(stderr()); flush(stdout())

### Name: plot_trend
### Title: Plots trends over time
### Aliases: plot_trend

### ** Examples

dat <- EMC2:::add_trials(forstmann)
dat$trials2 <- dat$trials/1000

lin_trend <- make_trend(cov_names='trials2',
                        kernels = 'exp_incr',
                        par_names='B',
                        bases='lin',
                        phase = "premap")

design_RDM_lin_B <- design(model=RDM,
                           data=dat,
                           covariates='trials2',   # specify relevant covariate columns
                           matchfun=function(d) d$S==d$lR,
                           transform=list(func=c('B'='identity')),
                           formula=list(B ~ 1, v ~ lM, t0 ~ 1),
                           trend=lin_trend)       # add trend

emc <- make_emc(dat, design=design_RDM_lin_B, compress = FALSE)
p_vector <- c('B'=1, 'v'=1, 'v_lMTRUE'=1, 't0'=0.1, 'B.w'=1, 'B.d_ei'=1)

# Visualize trend
plot_trend(p_vector, emc=emc,
           par_name='B', subject='as1t',
           filter=function(d) d$lR=='right', main='Threshold for right')



cleanEx()
nameEx("predict.emc")
### * predict.emc

flush(stderr()); flush(stdout())

### Name: predict.emc.prior
### Title: Generate Posterior/Prior Predictives
### Aliases: predict.emc.prior predict.emc

### ** Examples




cleanEx()
nameEx("printCompare")
### * printCompare

flush(stderr()); flush(stdout())

### Name: printCompare
### Title: Print a model-comparison table from compare()
### Aliases: printCompare

### ** Examples




cleanEx()
nameEx("prior")
### * prior

flush(stderr()); flush(stdout())

### Name: prior
### Title: Specify Priors for the Chosen Model
### Aliases: prior

### ** Examples

# First define a design for the model
design_DDMaE <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# Then set up a prior using prior
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
                     t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
                     t0=.4,Z=1,sv=.4,SZ=1)
# Here we left the variance prior at default
prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd)
# Also add a group-level variance prior:
pscale <- c(v_Sleft=.6,v_Sright=.6,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
                             t0=.2,Z=.5,sv=.4,SZ=.3)
df <- .4
prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd, A = pscale, df = df)
# If we specify a new design
design_DDMat0E <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~E, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# We can easily update the prior
prior_DDMat0E <- prior(design_DDMat0E, update = prior_DDMaE)



cleanEx()
nameEx("prior_help")
### * prior_help

flush(stderr()); flush(stdout())

### Name: prior_help
### Title: Prior Specification Information
### Aliases: prior_help

### ** Examples

prior_help('diagonal')



cleanEx()
nameEx("profile_plot")
### * profile_plot

flush(stderr()); flush(stdout())

### Name: profile_plot
### Title: Likelihood Profile Plots
### Aliases: profile_plot

### ** Examples




cleanEx()
nameEx("recovery")
### * recovery

flush(stderr()); flush(stdout())

### Name: recovery.emc
### Title: Recovery Plots
### Aliases: recovery.emc recovery

### ** Examples

# Make up some values that resemble posterior samples
# Normally this would be true values that were used to simulate the data
# Make up some values that resemble posterior samples
# Normally this would be true values that were used to simulate the data
pmat <- matrix(rnorm(12, mean = c(-1, -.6, -.4, -1.5), sd = .01), ncol = 4, byrow = TRUE)
# Conventionally this would be created before one makes data with true values
recovery(samples_LNR, pmat, correlation = "pearson", stat = "rmse", selection = "alpha")
# Similarly we can plot recovery of other parameters with a set of true samples
true_samples <- samples_LNR # Normally this would be data-generating samples
recovery(samples_LNR, true_samples, correlation = "pearson", stat = "rmse",
         selection = "correlation", cex = 1.5,
         ci_plot_args = list(lty = 3, length = .2, lwd = 2, col = "brown"))



cleanEx()
nameEx("reshape_events")
### * reshape_events

flush(stderr()); flush(stdout())

### Name: reshape_events
### Title: Reshape events data for fMRI analysis
### Aliases: reshape_events

### ** Examples

# Create a simple events data frame
events <- data.frame(
  subjects = rep(1, 10),
  run = rep(1, 10),
  onset = seq(0, 90, by = 10),
  condition = rep(c("A", "B"), 5),
  rt = runif(10, 0.5, 1.5),
  accuracy = sample(0:1, 10, replace = TRUE)
)

# Reshape with default duration
reshaped1 <- reshape_events(events, event_types = c("condition", "accuracy"))

# Reshape with custom duration for each event_type
reshaped2 <- reshape_events(events,
                           event_types = c("condition", "accuracy", "rt"),
                           duration = list(condition = 0.5,
                                          accuracy = 0.2,
                                          rt = function(x) x$rt))



cleanEx()
nameEx("run_bridge_sampling")
### * run_bridge_sampling

flush(stderr()); flush(stdout())

### Name: run_bridge_sampling
### Title: Estimating Marginal Likelihoods Using WARP-III Bridge Sampling
### Aliases: run_bridge_sampling

### ** Examples




cleanEx()
nameEx("run_emc")
### * run_emc

flush(stderr()); flush(stdout())

### Name: run_emc
### Title: Fine-Tuned Model Estimation
### Aliases: run_emc

### ** Examples




cleanEx()
nameEx("sampled_pars")
### * sampled_pars

flush(stderr()); flush(stdout())

### Name: sampled_pars
### Title: Get Model Parameters from a Design
### Aliases: sampled_pars sampled_pars.emc.design
###   sampled_pars.emc.group_design sampled_pars.emc.prior sampled_pars.emc

### ** Examples

# First define a design
design_DDMaE <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# Then for this design get which cognitive model parameters are sampled:
sampled_pars(design_DDMaE)




cleanEx()
nameEx("softmax")
### * softmax

flush(stderr()); flush(stdout())

### Name: softmax
### Title: Softmax Discrete Choice Model
### Aliases: softmax

### ** Examples

soft <- design(
  Rlevels = c("left","right"),
  factors = list(subjects = 1, S = c("left","right")),
  formula = list(x ~ 0 + S, beta ~ 1),
  matchfun = function(d) d$S == d$lR,
  constants = c(beta = log(1)),
  model = softmax
)

# Sample parameter vector
p_vector <- sampled_pars(soft)




cleanEx()
nameEx("split_timeseries")
### * split_timeseries

flush(stderr()); flush(stdout())

### Name: split_timeseries
### Title: Split fMRI Timeseries Data by ROI Columns
### Aliases: split_timeseries

### ** Examples

# Create a simple example timeseries with multiple ROIs
set.seed(123)
n_frames <- 100

# Create a data frame with multiple ROIs
timeseries <- data.frame(
  subjects = rep(1, n_frames),
  run = rep(1, n_frames),
  time = seq(0, n_frames-1),
  ROI1 = rnorm(n_frames),
  ROI2 = rnorm(n_frames),
  ROI3 = rnorm(n_frames)
)

# Split the timeseries by all ROI columns
split_data <- split_timeseries(timeseries)



cleanEx()
nameEx("subset.emc")
### * subset.emc

flush(stderr()); flush(stdout())

### Name: subset.emc
### Title: Shorten an emc Object
### Aliases: subset.emc

### ** Examples

subset(samples_LNR, length.out = 10)



cleanEx()
nameEx("summary.emc.prior")
### * summary.emc.prior

flush(stderr()); flush(stdout())

### Name: summary.emc.prior
### Title: Summary method for emc.prior objects
### Aliases: summary.emc.prior

### ** Examples

# Take a prior object
prior <- get_prior(samples_LNR)
summary(prior)




cleanEx()
nameEx("trend_help")
### * trend_help

flush(stderr()); flush(stdout())

### Name: trend_help
### Title: Get help information for trend kernels and bases
### Aliases: trend_help

### ** Examples

# Get information about exponential increasing kernel
trend_help(kernel = "exp_incr")

# Get information about linear base
trend_help(base = "lin")

# Return available kernel and base types
trend_help()



cleanEx()
nameEx("update2version")
### * update2version

flush(stderr()); flush(stdout())

### Name: update2version
### Title: Update EMC Objects to the Current Version
### Aliases: update2version

### ** Examples

# Update the model to current version
updated_model <- update2version(samples_LNR)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
