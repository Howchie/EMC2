## LBA likelihood + recovery demo for EMC2 to illustrate handling of censored and truncated RTs and optional contaminant omissions
# - Simulates 2-choice LBA data with a response deadline (UC) and optional truncation (LT,UT) with censored RTs and omissions coded as rt=Inf and R=NA
# - Runs a simulation with truncated RTs, generating profile plots (R+C likelihood) and a simple parameter-recovery fit, demonstrating the failure of the older likelihood with no censoring/truncation handling
# - Repeats the demo with censoring handling for comparison
# - Runs a simulation with censoring and truncated RTs, generating profile plots (R + C likelihood) and a simple parameter-recovery fit
# - Repeats the demo with contaminant omissions (e.g., 15%)
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


#### LBA ----

LBA_no_cens <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

LBA_cens <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_cens_trunc(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

# Andrew's R function, C versions require a different model MLBA/MRDM/MLNR
# e.g,. designLBA2 in demo
LBA_cens1 <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

run_lba_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,oldC=FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,print_stats=FALSE,
                         sample_file = "samples_LBA.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  capture.output(suppressMessages(designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LBA_cens, LBA_no_cens), # NB this is just for demonstration; normally use LBA with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  )),file=NULL)

  capture.output(suppressMessages(designLBA1 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LBA_cens1, LBA_no_cens), # NB this is just for demonstration; normally use LBA with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  )),file=NULL)

  capture.output(suppressMessages(
    if (!oldC) designLBA2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = MLBA,
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  ) else designLBA2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = OLBA,
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  )),file=NULL)


  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBA, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(2)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(.5)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("v" %in% names(p_vector)) p_vector[["v"]] <- 1
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.25

  # Sv for match accumulator
  if ("sv" %in% names(p_vector)) p_vector[["sv"]] <- log(1)
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["sv_lMTRUE"]] <- log(.8)

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
    UC = UC,
    UT = UT,
    LC = LC,
    LT = LT,
    LCresponse = LCresponse, UCresponse = UCresponse,
    LCdirection = LCdirection, UCdirection = UCdirection,
    force_direction = force_direction, force_response = force_response,
    rtContaminantNA = rtContaminantNA
  )

  if (print_stats) {

  cat("\n--- LBA demo ---\n")

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

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  # R _race_cens_trunc (Zachs)
  rtrt <- system.time({rtr <- profile_plot_test(dat, designLBA, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural)})

  # R _race_missing (mine)
  rtrt1 <- system.time({rtr1 <- profile_plot_test(dat, designLBA1, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R1 Likelihood",natural=natural)})

  #  C _race_cens_trunc
  rtct <- system.time({rtc <- profile_plot_test(dat, designLBA, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural)})

  # C _race_missing
  rtct1 <- system.time({rtc1 <- profile_plot_test(dat, designLBA2, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C1 Likelihood",natural=natural)})

  print(cbind(rtr,max1=rtr1$max,miss1=rtr1$miss,eq=abs(rtr1$miss)-abs(rtr$miss)))
  cat(paste0("R likelihood: ",round(rtrt[3],2),"\n"))
  cat(paste0("R1 likelihood: ",round(rtrt1[3],2),"\n"))
  cat(paste0("Speedup R/R1: ",round(rtrt[3]/rtrt1[3],2),"\n\n"))

  print(cbind(rtc,max1=rtc1$max,miss1=rtc1$miss,eq=abs(rtc1$miss)-abs(rtc$miss)))
  cat(paste0("C likelihood: ",round(rtct[3],2),"\n"))
  cat(paste0("C1 likelihood: ",round(rtct1[3],2),"\n"))
  cat(paste0("Speedup C/C1: ",round(rtct[3]/rtct1[3],2),"\n\n"))

  cat(paste0("\nSpeedup R/C: ",round(rtrt[3]/rtct[3],2),"\n"))
  cat(paste0("Speedup R1/C1: ",round(rtrt1[3]/rtct1[3],2),"\n"))
  cat(paste0("Speedup R1/C: ",round(rtrt1[3]/rtct[3],2),"\n"))

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


# Test 1: Equivalence of R and C++ function no censor/trunc/contam ----
# NB: Here and below cens=FALSE does not call the original Rcpp likelihood
#     (although it is still in particle_ll.cpp), the cens_trunc version is
#     ALWAYS called, whereas in R cens does switch between the orignal and
#     cens_trunc versions.

# max/miss is Zachs, max1/miss1 is mine
res_old <- run_lba_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  sample_file = "samples_lba_old.RData"
)

# C1 is now the old C version
res_old <- run_lba_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  oldC = TRUE,
  sample_file = "samples_lba_old.RData"
)

# Test 2: Censored RTs ----

# Defaults
resp <- FALSE
dirn <- TRUE

# # Opp default
# resp <- T
# dirn <- F
#
# # Opp res
# resp <- T
# dirn <- TRUE
#
# # Opp dirn
# resp <- FALSE
# dirn <- T


res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.5,
  UCresponse=resp,
  UCdirection=dirn,
  sample_file = "samples_lba_cens.RData"
)

res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .9,
  LCresponse=resp,
  LCdirection=dirn,
  sample_file = "samples_lba_cens.RData"
)


res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LCresponse=resp,
  UCresponse=resp,
  LCdirection=dirn,
  UCdirection=dirn,
  UC = 1.8,
  LC = .85,n_cores=9,
  sample_file = "samples_lba_cens.RData"
)


# Test 3: Truncated RTs ------------

res_contam <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.5,
  n_cores=9,
  sample_file = "samples_lba_cens_trunc.RData"
)

res_contam <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .9,
  n_cores=9,
  sample_file = "samples_lba_cens_trunc.RData"
)

res_contam <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.9,
  LT = .875,
  n_cores=9,
  sample_file = "samples_lba_cens_trunc.RData"
)

# Test 4: Truncation and censoring ------------

# MY R BAD!
res_contam <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.8,
  LC = .925,
  UT = 1.9,
  LT = .85,
  n_cores=9,
  sample_file = "samples_lba_cens_trunc.RData"
)

# Test 5: Contaminant omissions  ------------


res_contam <- run_lba_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  layout=c(2,4),
  sample_file = "samples_lba_contam.RData"
)

# My R better
res_contam <- run_lba_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.85,
  LT = .84,
  UT = 1.95,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_lba_cens_trunc.RData"
)

# My R worse
res_contam <- run_lba_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .915,
  LT = .84,
  UT = 1.95,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_lba_cens_trunc.RData"
)

# This case no longer barred by a check, as expected recovery not great but not super terrible
# My R worsete a lot., occasionaly qui
res_contam <- run_lba_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .915,
  UC = 1.85,
  LT = .84,
  UT = 1.95,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_lba_cens_trunc.RData"
)

#### RDM ----

RDM_no_cens <- function() {
  m <- RDM()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

RDM_cens <- function() {
  m <- RDM()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_cens_trunc(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

RDM_cens1 <- function() {
  m <- RDM()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

run_RDM_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,oldC=FALSE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,print_stats=FALSE,
                         sample_file = "samples_RDM.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  capture.output(suppressMessages(designRDM <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, RDM_cens, RDM_no_cens), # NB this is just for demonstration; normally use RDM with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1))
  )),file=NULL)

  capture.output(suppressMessages(designRDM1 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, RDM_cens1, RDM_no_cens), # NB this is just for demonstration; normally use RDM with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1))
  )),file=NULL)

  capture.output(suppressMessages(
    if (!oldC) designRDM2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = MRDM,
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1))
  ) else designRDM2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ORDM,
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1))
  )),file=NULL)


  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designRDM, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(2)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(.5)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("v" %in% names(p_vector)) p_vector[["v"]] <- log(1)
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- log(2)

  # s for match accumulator
  if ("s" %in% names(p_vector)) p_vector[["s"]] <- log(1)
  if ("s_lMTRUE" %in% names(p_vector)) p_vector[["s_lMTRUE"]] <- log(.8)

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
    p_vector, designRDM,
    n_trials = n_trials,
    UC = UC,
    UT = UT,
    LC = LC,
    LT = LT,
    LCresponse = LCresponse, UCresponse = UCresponse,
    LCdirection = LCdirection, UCdirection = UCdirection,
    force_direction = force_direction, force_response = force_response,
    rtContaminantNA = rtContaminantNA
  )

  if (print_stats) {

  cat("\n--- RDM demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x,na.rm=TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\nTruncation rate by stimulus:\n")
  print(1-tapply(dat$R, dat$S, length)/n_trials)

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  rtrt <- system.time({rtr <- profile_plot_test(dat, designRDM, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural)})

  rtrt1 <- system.time({rtr1 <- profile_plot_test(dat, designRDM1, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R1 Likelihood",natural=natural)})

  rtct <- system.time({rtc <- profile_plot_test(dat, designRDM, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural)})

  rtct1 <- system.time({rtc1 <- profile_plot_test(dat, designRDM2, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C1 Likelihood",natural=natural)})

  print(cbind(rtr,max1=rtr1$max,miss1=rtr1$miss,eq=abs(rtr1$miss)-abs(rtr$miss)))
  cat(paste0("R likelihood: ",round(rtrt[3],2),"\n"))
  cat(paste0("R1 likelihood: ",round(rtrt1[3],2),"\n"))
  cat(paste0("Speedup R/R1: ",round(rtrt[3]/rtrt1[3],2),"\n\n"))

  print(cbind(rtc,max1=rtc1$max,miss1=rtc1$miss,eq=abs(rtc1$miss)-abs(rtc$miss)))
  cat(paste0("C likelihood: ",round(rtct[3],2),"\n"))
  cat(paste0("C1 likelihood: ",round(rtct1[3],2),"\n"))
  cat(paste0("Speedup C/C1: ",round(rtct[3]/rtct1[3],2),"\n\n"))

  cat(paste0("\nSpeedup R/C: ",round(rtrt[3]/rtct[3],2),"\n"))
  cat(paste0("Speedup R1/C1: ",round(rtrt1[3]/rtct1[3],2),"\n"))
  cat(paste0("Speedup R1/C: ",round(rtrt1[3]/rtct[3],2),"\n"))

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designRDM, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designRDM, true_pars = p_vector))
}


RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Test 1: Equivalence of R and C++ function no censor/trunc/contam ----

res_old <- run_RDM_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  sample_file = "samples_RDM_old.RData"
)

res_old <- run_RDM_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  oldC = TRUE,
  sample_file = "samples_RDM_old.RData"
)


# Test 2: Censored RTs ----
res_cens <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.6,
  sample_file = "samples_RDM_cens.RData"
)

res_cens <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .85,
  sample_file = "samples_RDM_cens.RData"
)

res_cens <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.8,
  LC = .8,
  sample_file = "samples_RDM_cens.RData"
)


# Test 3: Truncated RTs ------------
res_contam <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.7,n_cores=9,
  sample_file = "samples_RDM_cens_trunc.RData"
)

# My R better
res_contam <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .9,n_cores=9,
  sample_file = "samples_RDM_cens_trunc.RData"
)

# My R better
res_contam <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.75,
  LT = .8,
  n_cores=9,
  sample_file = "samples_RDM_cens_trunc.RData"
)

# Test 4: Truncation and censoring ------------

# My R better
res_contam <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  n_cores=9,
  sample_file = "samples_RDM_cens_trunc.RData"
)

# Test 5: Contaminant omissions  ------------

res_contam <- run_RDM_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  layout=c(2,4),
  sample_file = "samples_RDM_contam.RData"
)

# My R mostly better
res_contam <- run_RDM_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.775,
  UT = 2.2,
  LT = .775,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_rdm_cens_trunc.RData"
)

# My R mostly better
res_contam <- run_RDM_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .825,
  UT = 2.2,
  LT = .775,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_rdm_cens_trunc.RData"
)

# My R mostly better
res_contam <- run_RDM_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  n_cores=9,
  layout=c(2,4),
  sample_file = "samples_rdm_cens_trunc.RData"
)


#### LNR ----

LNR_no_cens <- function() {
  m <- LNR()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

LNR_cens <- function() {
  m <- LNR()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_cens_trunc(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

LNR_cens1 <- function() {
  m <- LNR()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_missing(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

run_LNR_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,oldC=FALSE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,print_stats=FALSE,
                         sample_file = "samples_LNR.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  capture.output(suppressMessages(designLNR <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LNR_cens, LNR_no_cens), # NB this is just for demonstration; normally use LNR with built-in censoring/truncation handling
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
  )),file=NULL)

  capture.output(suppressMessages(designLNR1 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LNR_cens1, LNR_no_cens), # NB this is just for demonstration; normally use LNR with built-in censoring/truncation handling
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
  )),file=NULL)

  capture.output(suppressMessages(if (!oldC) designLNR2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = MLNR,
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
  ) else designLNR2 <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = OLNR,
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
  )),file=NULL)


  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLNR, doMap = FALSE)

  # non-decision time
   if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("m" %in% names(p_vector)) p_vector[["m"]] <- log(1.5)
  if ("m_lMTRUE" %in% names(p_vector)) p_vector[["m_lMTRUE"]] <- log(.45)

  # s for match accumulator
  if ("s" %in% names(p_vector)) p_vector[["s"]] <- log(1)
  if ("s_lMTRUE" %in% names(p_vector)) p_vector[["s_lMTRUE"]] <- log(.8)

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
    p_vector, designLNR,
    n_trials = n_trials,
    UC = UC,
    UT = UT,
    LC = LC,
    LT = LT,
    LCresponse = LCresponse, UCresponse = UCresponse,
    LCdirection = LCdirection, UCdirection = UCdirection,
    force_direction = force_direction, force_response = force_response,
    rtContaminantNA = rtContaminantNA
  )

  if (print_stats) {

  cat("\n--- LNR demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x,na.rm=TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\nTruncation rate by stimulus:\n")
  print(1-tapply(dat$R, dat$S, length)/n_trials)

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  rtrt <- system.time({rtr <- profile_plot_test(dat, designLNR, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural)})

  rtrt1 <- system.time({rtr1 <- profile_plot_test(dat, designLNR1, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R1 Likelihood",natural=natural)})

  rtct <- system.time({rtc <- profile_plot_test(dat, designLNR, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural)})

  rtct1 <- system.time({rtc1 <- profile_plot_test(dat, designLNR2, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C1 Likelihood",natural=natural)})

  print(cbind(rtr,max1=rtr1$max,miss1=rtr1$miss,eq=abs(rtr1$miss)-abs(rtr$miss)))
  cat(paste0("R likelihood: ",round(rtrt[3],2),"\n"))
  cat(paste0("R1 likelihood: ",round(rtrt1[3],2),"\n"))
  cat(paste0("Speedup R/R1: ",round(rtrt[3]/rtrt1[3],2),"\n\n"))

  print(cbind(rtc,max1=rtc1$max,miss1=rtc1$miss,eq=abs(rtc1$miss)-abs(rtc$miss)))
  cat(paste0("C likelihood: ",round(rtct[3],2),"\n"))
  cat(paste0("C1 likelihood: ",round(rtct1[3],2),"\n"))
  cat(paste0("Speedup C/C1: ",round(rtct[3]/rtct1[3],2),"\n\n"))

  cat(paste0("\nSpeedup R/C: ",round(rtrt[3]/rtct[3],2),"\n"))
  cat(paste0("Speedup R1/C1: ",round(rtrt1[3]/rtct1[3],2),"\n"))
  cat(paste0("Speedup R1/C: ",round(rtrt1[3]/rtct[3],2),"\n"))

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLNR, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designLNR, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Test 1: Equivalence of R and C++ function no censor/trunc/contam ----

res_old <- run_LNR_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  sample_file = "samples_LNR_old.RData"
)

res_old <- run_LNR_demo(
  p_contaminant = 0,
  cens=FALSE,oldC=TRUE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  sample_file = "samples_LNR_old.RData"
)

# Test 2: Censored RTs ----
res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.4,
  sample_file = "samples_LNR_cens.RData"
)

res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .45,
  n_cores=9,
  sample_file = "samples_LNR_cens.RData"
)


res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  n_cores=9,
  UC = 1.6,
  LC = .4,
  sample_file = "samples_LNR_cens.RData"
)


# Test 3: Truncated RTs ------------
res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.4,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)


res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .45,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)


res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.5,
  LT = .4,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)

# Test 4: Truncation and censoring ------------

# Uneven but my R a bit worse
res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.5,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)

# Test 5: Contaminant omissions  ------------

res_contam <- run_LNR_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  sample_file = "samples_LNR_contam.RData"
)

# My R worse
res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.5,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

# My R worse
res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

# My R worse
res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.5,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

