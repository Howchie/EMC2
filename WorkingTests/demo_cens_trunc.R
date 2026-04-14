## This script demonstrates censoring and/or truncation handling for all race models (LBA,RDM,LNR) and the DDM (all in C++)
# Also evidences defective cdf/density/stat calls work ok for data that includes omissions

rm(list=ls())
library(EMC2)
source("WorkingTests/test_likelihood_plotfuns_ah.R")


# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- TRUE

# Helper functions
Rtfun <- function(d) {
  out <- d$rt
  out[is.infinite(out)] <- NA
  out
}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)


#### LBA ----
run_lba_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 10000, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,cores_per_chain = 3,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         print_stats=TRUE,
                         label = NULL, posdrift=TRUE) {

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBA(posdrift=posdrift), 
    formula = c(
      list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ 1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(A = log(0.5))
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
  if ("sv_lMFALSE" %in% names(p_vector)) p_vector[["sv_lMFALSE"]] <- log(1)
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["sv_lMTRUE"]] <- log(0.8)

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
      UCresponse = UCresponse,LCdirection = LCdirection, UCdirection = UCdirection)
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
  
  if(!posdrift){
    cat("\nExpected Intrinsic Omissions:\n")
    print(pnorm(0,p_vector["v"],exp(p_vector["sv"]))*pnorm(0,p_vector["v_lMTRUE"],exp(p_vector["sv"])))
  }

  cat("\n")
  }

  # Profile plots. Use cores_per_chain = cores_per_chain on Windows.
  library(parallel)
  par(mfrow=layout)

  rtct <- system.time({rtc <- profile_plot_test(dat, designLBA, p_vector, n_cores = cores_for_chains, range=range,
  layout = NULL, use_c = TRUE,  figure_title = paste("C++:", label), natural=natural)})

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBA, type = "single")
    emc <- fit(emc,fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain=cores_per_chain, cores_for_chains = cores_for_chains), max_tries=30)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designLBA, true_pars = p_vector,emc=emc,pp=post_predict)))
  }

  invisible(list(data = dat, design = designLBA, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

res_no_cens_trunc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  label = "no_cens_trunc"
)

if (RUN_FITS) print(recovery(res_no_cens_trunc$emc, true_pars = res_no_cens_trunc$true_pars))
if (RUN_FITS) plot_cdf(res_no_cens_trunc$dat, post_predict=res_no_cens_trunc$pp, functions=list(Correct=Cfun), defective_factor = "Correct", factors="S")
if (RUN_FITS) plot_stat(res_no_cens_trunc$dat, post_predict=res_no_cens_trunc$pp, factors="S", stat_name = "MeanCorrect",
                        stat_fun = function(d){mean(d$Correct, na.rm = TRUE)}, functions=list(Correct=Cfun))
# Test 2: Censored RTs ----

# upper
res_uc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.5,
  label = "uc"
)
if (RUN_FITS) print(recovery(res_uc$emc, true_pars = res_uc$true_pars))
if (RUN_FITS) plot_cdf(res_uc$dat, post_predict=res_uc$pp, functions=list(Correct=Cfun), defective_factor = "Correct", factors="S")
if (RUN_FITS) plot_stat(res_uc$dat, post_predict=res_uc$pp, factors="S", stat_name = "MeanCorrect",
                        stat_fun = function(d){mean(d$Correct, na.rm = TRUE)}, functions=list(Correct=Cfun))
# lower
res_lc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .9,
  label = "lc"
)
if (RUN_FITS) print(recovery(res_lc$emc, true_pars = res_lc$true_pars))

# both
res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.8,
  LC = .85,
  label = "cens"
)

## Same plot-code can be used anywhere
if (RUN_FITS) print(recovery(res_cens$emc, true_pars = res_cens$true_pars))
if (RUN_FITS) plot_cdf(res_cens$data, post_predict = res_cens$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_density(res_cens$data, post_predict = res_cens$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_stat(res_cens$data, post_predict = res_cens$pp, factors = "S", stat_name = "MeanCorrect_Finite",
                        stat_fun = function(d) mean(d$Correct,na.rm=TRUE), functions = list(Correct = Cfun))
if (RUN_FITS) plot_stat(res_cens$data, post_predict = res_cens$pp, factors = "S", stat_name = "MeanOmissions",
                        stat_fun = function(d) mean(is.na(d$R)), functions = list(Correct = Cfun))
## LBAIO test with UC
res_lbaio <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 3,
  posdrift=FALSE,
  label = "lbaio"
)
if (RUN_FITS) print(recovery(res_lbaio$emc, true_pars = res_lbaio$true_pars))

# LBAIO with moderate lower censoring (to test the "LC means finite response" logic)
res_lbaio_lc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 3,
  LC = 0.8,
  posdrift=FALSE,
  label = "lbaio-uc"
)
if (RUN_FITS) print(recovery(res_lbaio_lc$emc, true_pars = res_lbaio_lc$true_pars))
# Test 3: Truncated RTs ------------

res_ut <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.5,
  label = "ut"
)
if (RUN_FITS) print(recovery(res_ut$emc, true_pars = res_ut$true_pars))
if (RUN_FITS) plot_cdf(res_ut$data, post_predict = res_ut$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_density(res_ut$data, post_predict = res_ut$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_stat(res_ut$data, post_predict = res_ut$pp, factors = "S", stat_name = "MeanCorrect_Finite",
                        stat_fun = function(d) mean(d$Correct,na.rm=TRUE), functions = list(Correct = Cfun))

res_lt <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .9,
  label = "lt"
)
if (RUN_FITS) print(recovery(res_lt$emc, true_pars = res_lt$true_pars))

res_trunc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.5,
  LT = .9,
  label = "trunc"
)
if (RUN_FITS) print(recovery(res_trunc$emc, true_pars = res_trunc$true_pars))

# Test 4: Truncation and censoring ------------

res_cens_trunc <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.8,
  LC = .925,
  UT = 1.9,
  LT = .85,
  label = "cens_trunc"
)
if (RUN_FITS) print(recovery(res_cens_trunc$emc, true_pars = res_cens_trunc$true_pars))
if (RUN_FITS) plot_cdf(res_cens_trunc$data, post_predict = res_cens_trunc$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_density(res_cens_trunc$data, post_predict = res_cens_trunc$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_stat(res_cens_trunc$data, post_predict = res_cens_trunc$pp, factors = "S", stat_name = "MeanCorrect_Finite",
                        stat_fun = function(d) mean(d$Correct,na.rm=TRUE), functions = list(Correct = Cfun))
if (RUN_FITS) plot_stat(res_cens_trunc$data, post_predict = res_cens_trunc$pp, factors = "S", stat_name = "MeanOmissions",
                        stat_fun = function(d) mean(is.na(d$R)), functions = list(Correct = Cfun))

# Test 5: Contaminant omissions  ------------
# test with no UC< this means rt=Inf is by definition a contaminant and gives a clean test of the likelihood
res_contam <- run_lba_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  label = "contam"
)
if (RUN_FITS) print(recovery(res_contam$emc, true_pars = res_contam$true_pars))

# test with finite uc (omissions now come from both censoring AND contamination and the model has to tell them apart)
res_contam_uc <- run_lba_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 2,
  posdrift=TRUE,
  label = "contam_lbaio"
)
if (RUN_FITS) print(recovery(res_contam_uc$emc, true_pars = res_contam_uc$true_pars))

# test with finite ut (I reordered the make_missing function so that truncation occurs before misses are injected so we can still observe them)
res_contam_ut <- run_lba_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UT = 2,
  posdrift=TRUE,
  label = "contam_lbaio"
)
if (RUN_FITS) print(recovery(res_contam_ut$emc, true_pars = res_contam_ut$true_pars))

# test lbaio with omissions (omissions can come from contam, negative drift, or upper censor)
res_contam_lbaio <- run_lba_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 2,
  posdrift=FALSE,
  label = "contam_lbaio"
)
if (RUN_FITS) print(recovery(res_contam_lbaio$emc, true_pars = res_contam_lbaio$true_pars))

# not sure this one makes sense, underestimates pContam but it is pretty heavily chopped up
res_cens_trunc_contam <- run_lba_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .915,
  UC = 1.85,
  LT = .84,
  UT = 1.95,
  label = "cens_trunc_contam"
)
if (RUN_FITS) print(recovery(res_cens_trunc_contam$emc, true_pars = res_cens_trunc_contam$true_pars))

# no truncation for lbaio
res_cens_contam_lbaio <- run_lba_demo(
  p_contaminant = .15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .915,
  UC = 1.85,
  posdrift=FALSE,
  label = "cens_contam_lbaio"
)
if (RUN_FITS) print(recovery(res_cens_contam_lbaio$emc, true_pars = res_cens_contam_lbaio$true_pars))


#### DDM using Henrich & Klauer (2026) ----

run_DDM_demo <- function(
    n_trials = 10000, UC = NULL, UT = NULL, LC = NULL, LT = NULL, range = 1,
    cores_for_chains = 3, layout = c(2, 3), natural = TRUE, cores_per_chain = 3,
    LCresponse = FALSE, UCresponse = FALSE,
    LCdirection = TRUE, UCdirection = TRUE,
    print_stats = TRUE,
    label = NULL
) {

  designDDM <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    model = DDM,
    formula = list(v ~ 0 + S, a ~ 1, t0 ~ 1, Z ~ 1),
    constants = c(s = log(1), sv = log(0), SZ = qnorm(0), st0 = log(0))
  )

  p_vector <- sampled_pars(designDDM, doMap = FALSE)
  if ("v_Sleft" %in% names(p_vector)) p_vector[["v_Sleft"]] <- -1.2
  if ("v_Sright" %in% names(p_vector)) p_vector[["v_Sright"]] <- 1.2
  if ("a" %in% names(p_vector)) p_vector[["a"]] <- log(2)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)
  if ("Z" %in% names(p_vector)) p_vector[["Z"]] <- qnorm(0.6)

  dat <- make_data(
    p_vector, designDDM,
    n_trials = n_trials,
    TC = list(
      UC = UC, UT = UT, LC = LC, LT = LT,
      LCresponse = LCresponse, UCresponse = UCresponse,
      LCdirection = LCdirection, UCdirection = UCdirection
    )
  )

  if (print_stats) {
    cat("\n--- DDM demo ---\n")
    cat("Counts by response (NA = omission):\n")
    print(table(dat$R, useNA = "ifany"))
    cat("\nMean RT by stimulus (finite RTs only):\n")
    print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))
    cat("\nMean Acc by stimulus (defined responses only):\n")
    print(tapply(Cfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))
    cat("\nOmission rate by stimulus:\n")
    print(tapply(is.na(dat$R), dat$S, mean))
    cat("\nTruncation rate by stimulus:\n")
    print(1 - (tapply(dat$R, dat$S, length) / n_trials))
    cat("\n")
  }

  library(parallel)
  par(mfrow = layout)
  rtct <- system.time({
    rtc <- profile_plot_test(
      dat, designDDM, p_vector,
      n_cores = cores_for_chains, range = range,
      layout = NULL, use_c = TRUE,
      figure_title = paste("C++:", label), natural = natural
    )
  })

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designDDM, type = "single")
    emc <- fit(
      emc,
      fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),
      stop_criteria = list(
        sample = list(
          iter = 1000,
          max_gd = 1.10,
          max_flat_loc = 0.5,
          flat_selection = c("alpha", "subj_ll"),
          flat_p1 = 1/3,
          flat_p2 = 1/3,
          max_sample_iter = 5000
        ),
        cores_per_chain = cores_per_chain,
        cores_for_chains = cores_for_chains
      ),
      max_tries = 30
    )
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designDDM, true_pars = p_vector, emc = emc, pp = post_predict)))
  }

  invisible(list(data = dat, design = designDDM, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# 1. uncensored/truncated
res_ddm_uncensored <- run_DDM_demo(
  n_trials = 10000,
  label = "ddm_uncensored"
)
if (RUN_FITS) print(recovery(res_ddm_uncensored$emc, true_pars = res_ddm_uncensored$true_pars))

# 2. upper censor only
res_ddm_uc <- run_DDM_demo(
  n_trials = 10000,
  UC = 1.5,
  label = "ddm_uc"
)

## Same plot-code can be used anywhere
if (RUN_FITS) print(recovery(res_ddm_uc$emc, true_pars = res_ddm_uc$true_pars))
if (RUN_FITS) plot_cdf(res_ddm_uc$data, post_predict = res_ddm_uc$pp, functions = list(Correct = Cfun), defective_factor = "Correct", factors = "S")
if (RUN_FITS) plot_stat(res_ddm_uc$data, post_predict = res_ddm_uc$pp, factors = "S", stat_name = "MeanCorrect_Finite",
                        stat_fun = function(d) mean(d$Correct,na.rm=TRUE), functions = list(Correct = Cfun))
if (RUN_FITS) plot_stat(res_ddm_uc$data, post_predict = res_ddm_uc$pp, factors = "S", stat_name = "MeanOmissions",
                        stat_fun = function(d) mean(is.na(d$R)), functions = list(Correct = Cfun))
# 3. lower censor only
res_ddm_lc <- run_DDM_demo(
  n_trials = 10000,
  LC = 0.5,
  label = "ddm_lc"
)
if (RUN_FITS) print(recovery(res_ddm_lc$emc, true_pars = res_ddm_lc$true_pars))

# 4. both censor bounds
res_ddm_both <- run_DDM_demo(
  n_trials = 10000,
  LC = 0.55,
  UC = 1.3,
  label = "ddm_both"
)
if (RUN_FITS) print(recovery(res_ddm_both$emc, true_pars = res_ddm_both$true_pars))

## TODO fix truncation
# 5. upper trunc only
res_ddm_ut <- run_DDM_demo(
  n_trials = 10000,
  UT = 1.5,
  label = "ddm_ut"
)
if (RUN_FITS) print(recovery(res_ddm_ut$emc, true_pars = res_ddm_ut$true_pars))

# 6. lower trunc only
res_ddm_lt <- run_DDM_demo(
  n_trials = 10000,
  LT = 0.35,
  label = "ddm_lt"
)
if (RUN_FITS) print(recovery(res_ddm_lt$emc, true_pars = res_ddm_lt$true_pars))

# 7. all: trunc + censor
res_ddm_all <- run_DDM_demo(
  n_trials = 10000,
  LT = 0.30,
  LC = 0.55,
  UC = 1.35,
  UT = 1.8,
  label = "ddm_all"
)
if (RUN_FITS) print(recovery(res_ddm_all$emc, true_pars = res_ddm_all$true_pars))


#### RDM ----

run_RDM_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,cores_per_chain = 3,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         print_stats=TRUE,
                         label = NULL) {

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designRDM <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = RDM, # NB this is just for demonstration; normally use RDM with built-in censoring/truncation handling
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
    TC=list(UC = UC,UT = UT,LC = LC,LT = LT,LCresponse = LCresponse,
      UCresponse = UCresponse,LCdirection = LCdirection, UCdirection = UCdirection)
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

  # Profile plots. Use cores_per_chain = cores_per_chain on Windows.
  library(parallel)
  par(mfrow=layout)

  rtct <- system.time({rtc <- profile_plot_test(dat, designRDM, p_vector, n_cores = cores_for_chains, range=range,
  layout = NULL, use_c = TRUE, figure_title = paste("C++:", label),natural=natural)})

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designRDM, type = "single")
    emc <- fit(emc,fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains), max_tries=30)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designRDM, true_pars = p_vector,emc=emc,pp=post_predict)))
  }

  invisible(list(data = dat, design = designRDM, true_pars = p_vector))
}


RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# Test 1 - Uncensored

res_rdm <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  label = "rdm"
)
if (RUN_FITS) print(recovery(res_rdm$emc, true_pars = res_rdm$true_pars))

# Test 2: Censored RTs ----

res_rdm_cens <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.7,
  LC = .6,
  label = "rdm_cens"
)
if (RUN_FITS) print(recovery(res_rdm_cens$emc, true_pars = res_rdm_cens$true_pars))

res_rdm_trunc <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.7,
  LT=.6,
  label = "rdm_trunc"
)
if (RUN_FITS) print(recovery(res_rdm_trunc$emc, true_pars = res_rdm_trunc$true_pars))

res_rdm_cens_trunc <- run_RDM_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  label = "rdm_cens_trunc"
)
if (RUN_FITS) print(recovery(res_rdm_cens_trunc$emc, true_pars = res_rdm_cens_trunc$true_pars))

res_rdm_contam <- run_RDM_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  label = "rdm_contam"
)
if (RUN_FITS) print(recovery(res_rdm_contam$emc, true_pars = res_rdm_contam$true_pars))

res_rdm_cens_trunc_contam <- run_RDM_demo(
  p_contaminant = 0.15,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  label = "rdm_cens_trunc_contam"
)
if (RUN_FITS) print(recovery(res_rdm_cens_trunc_contam$emc, true_pars = res_rdm_cens_trunc_contam$true_pars))

#### LNR ----

run_LNR_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,cores_per_chain = 3,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         print_stats=TRUE,
                         label = NULL) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLNR <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LNR, 
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
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
    TC=list(UC = UC,UT = UT,LC = LC,LT = LT,LCresponse = LCresponse,
      UCresponse = UCresponse,LCdirection = LCdirection, UCdirection = UCdirection)
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

  # Profile plots. Use cores_per_chain = cores_per_chain on Windows.
  library(parallel)
  par(mfrow=layout)

  rtct <- system.time({rtc <- profile_plot_test(dat, designLNR, p_vector, n_cores = cores_for_chains, range=range,
  layout = NULL, use_c = TRUE, figure_title = paste("C++:", label),natural=natural)})

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLNR, type = "single")
    emc <- fit(emc,fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains), max_tries=30)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designLNR, true_pars = p_vector,emc=emc,pp=post_predict)))
  }

  invisible(list(data = dat, design = designLNR, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Test 1 - Uncensored

res_lnr <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  label = "lnr"
)
if (RUN_FITS) print(recovery(res_lnr$emc, true_pars = res_lnr$true_pars))

# Test 2: Censored RTs ----

res_lnr_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.6,
  LC = .85,
  label = "lnr_cens"
)
if (RUN_FITS) print(recovery(res_lnr_cens$emc, true_pars = res_lnr_cens$true_pars))

res_lnr_trunc <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.7,
  LT=.6,
  label = "lnr_trunc"
)
if (RUN_FITS) print(recovery(res_lnr_trunc$emc, true_pars = res_lnr_trunc$true_pars))

res_lnr_cens_trunc <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  label = "lnr_cens_trunc"
)
if (RUN_FITS) print(recovery(res_lnr_cens_trunc$emc, true_pars = res_lnr_cens_trunc$true_pars))

res_lnr_contam <- run_LNR_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  label = "lnr_contam"
)
if (RUN_FITS) print(recovery(res_lnr_contam$emc, true_pars = res_lnr_contam$true_pars))

res_lnr_cens_trunc_contam <- run_LNR_demo(
  p_contaminant = 0.15,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.775,
  LC = .825,
  UT = 2.2,
  LT = .775,
  label = "lnr_cens_trunc_contam"
)
if (RUN_FITS) print(recovery(res_lnr_cens_trunc_contam$emc, true_pars = res_lnr_cens_trunc_contam$true_pars))
save.image("demo_cens_trunc.RData")
