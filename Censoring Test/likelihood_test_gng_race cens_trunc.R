## LBA likelihood + recovery demo for EMC2 to illustrate handling of censored and truncated RTs and optional contaminant omissions
# - Simulates 2-choice LBA Go-No-Go data with a response deadline (UC) and optional truncation (LT,UT) with censored RTs and omissions coded as rt=Inf and R=NA
# - Runs a simulation with truncated RTs, generating profile plots (R+C likelihood) and a simple parameter-recovery fit, demonstrating the failure of the older likelihood with no censoring/truncation handling
# - Repeats the demo with censoring handling for comparison
# - Runs a simulation with censoring and truncated RTs, generating profile plots (R + C likelihood) and a simple parameter-recovery fit
# - Repeats the demo with contaminant omissions (e.g., 15%)

library(EMC2)
source("Censoring Test/test_likelihood_plotfuns.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- FALSE

# Helper functions
Rtfun <- function(d){ d$rt[d$rt==Inf]=NA; d$rt}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)

run_gng_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = Inf, UT=Inf, LC=0, LT=0,
                         cores_for_chains = 3,
                         sample_file = "samples_GNG.RData") {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  designLBA <- design(
    factors=list(subjects=1,S=c("go","go","go","nogo")),
    Rlevels=c("go","nogo"),
    matchfun=matchfun,
    model=LBAGNG,
    formula = c(
      list(B ~ 0 + lM, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ 0 + lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv_lMFALSE = log(1))

  )

  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBA, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B_Sgo" %in% names(p_vector)) p_vector[["B_Sgo"]] <- log(1.0)
  if ("B_Snogo" %in% names(p_vector)) p_vector[["B_Snogo"]] <- log(1.5)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(0.5)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.25)

  # Drift means for mismatch vs match accumulators
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.2
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.1

  # Sv for match accumulator
  if ("sv_lMFALSE" %in% names(p_vector)) p_vector[["sv_lMFALSE"]] <- log(1)
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["sv_lMTRUE"]] <- log(1)

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
    LCdirection = FALSE, # deadline-only censoring -> rt=Inf (not -Inf)
    rtContaminantNA = FALSE # contaminant omissions coded same as censored with Inf
  )
  dat$R[is.na(dat$R)] <- "nogo" # Ensure NAs are coded as "nogo" responses
  cat("\n--- LBA demo ---\n")
  cat("estimate_contaminant =", estimate_contaminant, "\n")
  cat("p_contaminant (true) =", p_contaminant, "\n")
  cat("UC =", UC, " n_trials per cell =", n_trials, "\n\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x,na.rm=TRUE)))

  ## Print some sanity checks to show accuracy and RT by cell (showing there are errors)
  print("Accuracy by Cell: ")
  print(tapply(Cfun(dat),dat$S,function(x){mean(x,na.rm=TRUE)}))
  print("RT by Cell: ")
  print(tapply(Rtfun(dat),dat$S,function(x){mean(x,na.rm=TRUE)}))

  ## Sanity check that there are no "no-go" responses with an RT, and no "go" responses with an undefined RT
  print("RT by Response: ")
  print(tapply(Rtfun(dat),dat[,c("R")],function(x){mean(x,na.rm=TRUE)}))
  print("Proportion of Undefined RTs by Response: ")
  print(tapply(Rtfun(dat),dat[,c("R")],function(x){mean(is.na(x))}))

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  profile_plot_test(dat, designLBA, p_vector, n_cores = 1, range=1,layout = NA, figure_title = "R Likelihood") # R likelihood
  profile_plot_test(dat, designLBA, p_vector, n_cores = 1, range=1,layout = NA, use_c = TRUE, figure_title = "C Likelihood") # C likelihood

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBA, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designLBA, true_pars = p_vector))
}


# --- Test 1: LBAGNG recovery with both R and C likelihoods using censoring/truncation handling to demonstrate equivalence
res_cens <- run_gng_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 500,
  UC = 2.0, # include censoring and stops
  sample_file = "samples_gng_race.RData"
)

# --- Test 2: LBA recovery with truncated RTs and censoring ------------
res_contam <- run_gng_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 500,
  UC = 2.0,
  UT = 3.0,
  LT = 0.2,
  LC = 0.2,
  sample_file = "samples_gng_race_cens_trunc.RData"
)

# --- Test 3: LBA recovery with contaminant omissions (e.g., 15%) ------------
res_contam <- run_gng_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 500,
  UC = 2.0,
  sample_file = "samples_lba_contam.RData"
)
