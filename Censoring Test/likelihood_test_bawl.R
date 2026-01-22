# Test of Ballistic Accumulator with Leak

rm(list=ls())
library(EMC2)
source("Censoring Test/test_likelihood_plotfuns_ah.R")


# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- FALSE

# Helper functions
Rtfun <- function(d) {
  out <- d$rt
  out[is.infinite(out)] <- NA
  out
}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)


run_bawl_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,oldC=FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         print_stats=FALSE,
                         sample_file = "samples_LBA.RData",
                         posdrift=TRUE) {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  capture.output(suppressMessages(designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL(posdrift=posdrift), # NB this is just for demonstration; normally use LBA with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM, k~0+lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv_lMFALSE = log(1.2))
  )),file=NULL)

  
  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBA, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(1.6)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(.4)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.15)

  # Drift means for mismatch vs match accumulators
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- .65
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1

  # Sv for match accumulator
  if ("sv_lMFALSE" %in% names(p_vector)) p_vector[["sv_lMFALSE"]] <- log(1.2)
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["sv_lMTRUE"]] <- log(.8)

  # Leakage
  if ("k_lMFALSE" %in% names(p_vector)) p_vector[["k_lMFALSE"]] <- log(.2)
  if ("k_lMTRUE" %in% names(p_vector)) p_vector[["k_lMTRUE"]] <- log(1)
  
  # Contaminant omissions: pContaminant is on probit scale (transformed via pnorm)
  if (estimate_contaminant) {
    if (!("pContaminant" %in% names(p_vector))) {
      stop("estimate_contaminant=TRUE but pContaminant not found in sampled parameters.")
    }
    if (p_contaminant <= 0 || p_contaminant >= 1) stop("p_contaminant must be in (0,1) when estimated.")
    p_vector[["pContaminant"]] <- qnorm(p_contaminant)
  } else {
    # Fixed to 0 contaminant probability (i.e., pnorm(-Inf) = 0)
    if ("pContaminant" %in% names(p_vector)) p_vector[["pContaminant"]] <- qnorm(0)
  }

  dat <- make_data(
    p_vector, designLBA,
    n_trials = n_trials,
    TC=list(UC = UC,UT = UT,LC = LC,LT = LT,LCresponse = LCresponse,
      UCresponse = UCresponse,LCdirection = LCdirection, UCdirection = UCdirection, verbose = TRUE)
  )

  if (print_stats) {

  cat("\n--- BAwL demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x,na.rm=TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\nTruncation rate by stimulus:\n")
  print(1-(tapply(dat$R, dat$S, length)/n_trials))
  
  # Added extra output here. On average I see ~15% censored at UC=3, and ~ 7.5% IOs, and total omissions print around ~23% so it seems like it's working
  cat("\nExpected Failed Starts:\n")
    print(ifelse(posdrift,0,
    pnorm(0,p_vector["v_lMFALSE"],exp(designLBA$constants["sv_lMFALSE"]))*pnorm(0,p_vector["v_lMTRUE"],exp(p_vector["sv_lMTRUE"])))
  )

  cat("\nExpected Weak Finishes (removing failed):\n")
  print((pnorm(p_vector["k_lMFALSE"] * (p_vector["A"] + p_vector["B"]), p_vector["v_lMFALSE"], 
    exp(designLBA$constants["sv_lMFALSE"])) * pnorm(p_vector["k_lMTRUE"] * (p_vector["A"] + p_vector["B"]), p_vector["v_lMTRUE"], exp(p_vector["sv_lMTRUE"])))
    - pnorm(0,p_vector["v_lMFALSE"],exp(designLBA$constants["sv_lMFALSE"]))*pnorm(0,p_vector["v_lMTRUE"],exp(p_vector["sv_lMTRUE"])))

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  #  C _race_cens_trunc
  rtct <- system.time({rtc <- profile_plot_test(dat, designLBA, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural)})

  print(cbind(rtc))
  
  cat(paste0("C likelihood: ",round(rtct[3],2),"\n"))


  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBA, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designLBA, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Test 1: Censored RTs ----

# Defaults
resp <- FALSE
dirn <- TRUE

# Gentle Censor at 10s
res <- run_bawl_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 100,
  UCresponse=resp,
  UCdirection=dirn,
  sample_file = "samples_lba_cens.RData",
  print_stats=TRUE,
  posdrift=TRUE
)

par(mfrow=c(1,1))
hist(res$data$rt,breaks=50,main="RT Histogram with Censoring at 10s",xlab="RT (s)")
print(paste0("Proportion Censored: ",round(mean(res$data$rt>=10),10)))
# Gentle Censor at 10s + posdrift=FALSE
res <- run_bawl_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 100,
  UCresponse=resp,
  UCdirection=dirn,
  sample_file = "samples_lba_cens.RData",
  print_stats=TRUE,
  posdrift=FALSE
)
par(mfrow=c(1,1))
hist(res$data$rt,breaks=50,main="RT Histogram with Censoring at 10s",xlab="RT (s)")
print(paste0("Proportion Censored: ",round(mean(res$data$rt>=10),10)))

# Harsh lower censor
res <- run_bawl_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .9,
  UC = 100, # must be a UC otherwise Intrinsic Omissions aren't handled properly
  LCresponse=resp,
  LCdirection=dirn,
  sample_file = "samples_lba_cens.RData",
  print_stats=TRUE
)

# Combined
res <- run_bawl_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LCresponse=resp,
  UCresponse=resp,
  LCdirection=dirn,
  UCdirection=dirn,
  UC = 1.8,
  LC = .85,n_cores=9,
  sample_file = "samples_lba_cens.RData",
  print_stats=TRUE
)


# Test 2: Contaminant omissions  ------------

res_contam <- run_bawl_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC=100,
  layout=c(2,4),n_cores=9,
  sample_file = "samples_lba_contam.RData",
  print_stats=TRUE
)


# Contaminant with upper and lower censor
res_contam <- run_bawl_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .8,
  UC = 2.25,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_lba_cens_trunc.RData",
  print_stats=TRUE
)

## Sampling Check
# Here is a simple 5 parameter LBAIO model with contamination that uses Zach's C
designLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=function(x){LBA(posdrift=FALSE)}, # ZH: We need to find a way for the design to allow function arguments, it breaks with a Error in model() : could not find function "model" error when there's any arguments passed
  formula=list(v~0+lM,B~1,t0~1,sv~0+lM,pContaminant~1),
  constants=c(A=log(.5)) #
)
# This parameter vector produces reasonable RTs, lowish accuracy and a non-zero percent of intrinsic omissions.
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(.5, .825, log(1.6), log(0.15),log(1), log(0.8), qnorm(.05))

TC=list(LT=0,UT=Inf,LC=0,UC=3,verbose=TRUE)
datCT <- make_data(p_vector,designLBA,n_trials=10000,TC=TC)

# Fitting -- check memory usage and speed equivalence
emc <- make_emc(datCT,designLBA,type="single", compress=TRUE)
emccCT <- fit(emc)
print(recovery(emccCT,p_vector,selection="alpha"))

## Harder fit
p_vector[] <- c(.5, .825, log(1.6), log(0.15),log(1), log(0.8), qnorm(.1))

TC=list(LT=0,UT=Inf,LC=0.65,UC=2.25,verbose=TRUE)
datCT <- make_data(p_vector,designLBA,n_trials=10000,TC=TC)

# Fitting -- check memory usage and speed equivalence
emc <- make_emc(datCT,designLBA,type="single", compress=TRUE)
emccCT <- fit(emc)
print(recovery(emccCT,p_vector,selection="alpha"))
