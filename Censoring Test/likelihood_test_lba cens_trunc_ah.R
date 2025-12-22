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

run_lba_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,
                         sample_file = "samples_LBA.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LBA_cens, LBA_no_cens), # NB this is just for demonstration; normally use LBA with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  )

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
  print(1-tapply(dat$R, dat$S, length)/n_trials)

  cat("\n")

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  rtr <- system.time({print(profile_plot_test(dat, designLBA, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural))})
  cat(paste0("R likelihood: ",round(rtr[3],2),"\n\n"))

  rtc <- system.time({print(profile_plot_test(dat, designLBA, p_vector, n_cores = n_cores, range=range,
    layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural))})
  cat(paste0("C likelihood: ",round(rtc[3],2),"\n"))
  cat(paste0("\nSpeedup: ",round(rtr[3]/rtc[3],2),"\n"))

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

res_old <- run_lba_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  # UC = Inf,
  sample_file = "samples_lba_old.RData"
)
#           true   max   miss
# B         2.00 1.979  0.021
# v         1.00 1.015 -0.015
# v_lMTRUE  1.25 1.269 -0.019
# A         0.50 0.468  0.032
# t0        0.20 0.192  0.008
# sv_lMTRUE 0.80 0.815 -0.015
# R likelihood: 3.06
#
#           true   max   miss
# B         2.00 1.979  0.021
# v         1.00 1.015 -0.015
# v_lMTRUE  1.25 1.269 -0.019
# A         0.50 0.468  0.032
# t0        0.20 0.192  0.008
# sv_lMTRUE 0.80 0.815 -0.015
# C likelihood: 2.95
#
# Speedup: 1.04

# # SPEED UNDER DEV
#           true   max   miss
# B         2.00 1.979  0.021
# v         1.00 1.015 -0.015
# v_lMTRUE  1.25 1.269 -0.019
# A         0.50 0.468  0.032
# t0        0.20 0.192  0.008
# sv_lMTRUE 0.80 0.815 -0.015
# R likelihood: 2.92
#
#           true   max   miss
# B         2.00 1.979  0.021
# v         1.00 1.015 -0.015
# v_lMTRUE  1.25 1.269 -0.019
# A         0.50 0.468  0.032
# t0        0.20 0.192  0.008
# sv_lMTRUE 0.80 0.815 -0.015
# C likelihood: 2.4
#
# Speedup: 1.22



# Test 2: Censored RTs ----
res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.5,
  sample_file = "samples_lba_cens.RData"
)

res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .9,
  sample_file = "samples_lba_cens.RData"
)


res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
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
  p_contaminant = 0.05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  layout=c(2,4),
  sample_file = "samples_lba_contam.RData"
)

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

# This case no longer barred by a check, as expected recovery not great but not terrible
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

run_RDM_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,
                         sample_file = "samples_RDM.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designRDM <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, RDM_cens, RDM_no_cens), # NB this is just for demonstration; normally use RDM with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1))
  )

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

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)
  rtr <- system.time({print(profile_plot_test(dat, designRDM, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural))})
  cat(paste0("R likelihood: ",round(rtr[3],2),"\n\n"))


  rtc <- system.time({print(profile_plot_test(dat, designRDM, p_vector, n_cores = n_cores, range=range,
    layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural))})
  cat(paste0("C likelihood: ",round(rtc[3],2),"\n"))
  cat(paste0("\nSpeedup: ",round(rtr[3]/rtc[3],2),"\n"))

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
  UC = Inf,
  sample_file = "samples_RDM_old.RData"
)
#          true   max   miss
# B         2.0 1.979  0.021
# v         1.0 1.005 -0.005
# v_lMTRUE  2.0 2.021 -0.021
# A         0.5 0.468  0.032
# t0        0.2 0.192  0.008
# s_lMTRUE  0.8 0.810 -0.010
# R likelihood: 3.02
#
#          true   max   miss
# B         2.0 1.979  0.021
# v         1.0 1.005 -0.005
# v_lMTRUE  2.0 2.021 -0.021
# A         0.5 0.468  0.032
# t0        0.2 0.192  0.008
# s_lMTRUE  0.8 0.810 -0.010
# C likelihood: 2.81
#
# Speedup: 1.07

# SPEED UNDER DEV
#          true   max   miss
# B         2.0 1.979  0.021
# v         1.0 1.005 -0.005
# v_lMTRUE  2.0 2.021 -0.021
# A         0.5 0.468  0.032
# t0        0.2 0.192  0.008
# s_lMTRUE  0.8 0.810 -0.010
# R likelihood: 2.74
#
#          true   max   miss
# B         2.0 1.979  0.021
# v         1.0 1.005 -0.005
# v_lMTRUE  2.0 2.021 -0.021
# A         0.5 0.468  0.032
# t0        0.2 0.192  0.008
# s_lMTRUE  0.8 0.810 -0.010
# C likelihood: 2.39
#
# Speedup: 1.15


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

res_contam <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .9,n_cores=9,
  sample_file = "samples_RDM_cens_trunc.RData"
)

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

# devtools::load_all()
# debug(run_RDM_demo)
# undebug(EMC2:::log_likelihood_race_cens_trunc)
# undebug(profile_plot_test)
# undebug(lfun)
# undebug(EMC2:::calc_ll_R)
# undebug(EMC2:::log_likelihood_race_cens_trunc)



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

run_LNR_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         rtContaminantNA = FALSE,
                         sample_file = "samples_LNR.RData") {
  # if (UC == Inf & p_contaminant>0)
  #   stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLNR <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, LNR_cens, LNR_no_cens), # NB this is just for demonstration; normally use LNR with built-in censoring/truncation handling
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    )
  )

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

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)
  rtr <- system.time({print(profile_plot_test(dat, designLNR, p_vector, n_cores = n_cores, range=range,
    layout = NULL, figure_title = "R Likelihood",natural=natural))})
  cat(paste0("R likelihood: ",round(rtr[3],2),"\n\n"))


  rtc <- system.time({print(profile_plot_test(dat, designLNR, p_vector, n_cores = n_cores, range=range,
    layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural))})
  cat(paste0("C likelihood: ",round(rtc[3],2),"\n"))
  cat(paste0("\nSpeedup: ",round(rtr[3]/rtc[3],2),"\n"))

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
  UC = Inf,
  sample_file = "samples_LNR_old.RData"
)
#            true    max   miss
# m         0.405  0.391  0.014
# m_lMTRUE -0.799 -0.811  0.012
# s         1.000  1.015 -0.015
# s_lMTRUE  0.800  0.814 -0.014
# t0        0.200  0.192  0.008
# R likelihood: 2.09
#
#            true    max   miss
# m         0.405  0.391  0.014
# m_lMTRUE -0.799 -0.811  0.012
# s         1.000  1.015 -0.015
# s_lMTRUE  0.800  0.814 -0.014
# t0        0.200  0.192  0.008
# C likelihood: 2.04
#
# Speedup: 1.02

# SPEED UNDER DEV
#            true    max   miss
# m         0.405  0.391  0.014
# m_lMTRUE -0.799 -0.811  0.012
# s         1.000  1.015 -0.015
# s_lMTRUE  0.800  0.814 -0.014
# t0        0.200  0.192  0.008
# R likelihood: 1.7
#
#            true    max   miss
# m         0.405  0.391  0.014
# m_lMTRUE -0.799 -0.811  0.012
# s         1.000  1.015 -0.015
# s_lMTRUE  0.800  0.814 -0.014
# t0        0.200  0.192  0.008
# C likelihood: 1.65
#
# Speedup: 1.03

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
