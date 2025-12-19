## LBA likelihood + recovery demo for EMC2 to illustrate handling of censored and truncated RTs and optional contaminant omissions
# - Simulates 2-choice LBA data with a response deadline (UC) and optional truncation (LT,UT) with censored RTs and omissions coded as rt=Inf and R=NA
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
Rtfun <- function(d) {
  out <- d$rt
  out[is.infinite(out)] <- NA
  out
}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
lba_no_cens <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    log_likelihood_race(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}
  
lba_cens <- function() {
  m <- LBA()
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    log_likelihood_race_cens_trunc(
      pars, dadm, model, min_ll = min_ll
    )
  }
  m
}

run_lba_demo <- function(p_contaminant = 0, cens = TRUE, estimate_contaminant = FALSE,
                         n_trials = 500, UC = 2.0, UT=Inf, LC=0, LT=0,
                         cores_for_chains = 3,
                         sample_file = "samples_LBA.RData") {
  if (UC == Inf & p_contaminant>0) stop("Require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = ifelse(cens, lba_cens, lba_no_cens), # NB this is just for demonstration; normally use LBA with built-in censoring/truncation handling
    formula = c(
      list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ 0 + lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv_lMFALSE = log(1))
  )

  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBA, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(1.0)
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

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

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

# --- Test 1a: Demonstrate equivalence of new C++ function when no censoring present (compares against old R likelihood)
res_old <- run_lba_demo(
  p_contaminant = 0,
  cens=FALSE,
  estimate_contaminant = FALSE,
  n_trials = 500,
  UC = Inf,
  sample_file = "samples_lba_old.RData"
)
# --- Test 1b: Repeat above with both new R and C likelihoods using censoring handling to demonstrate equivalence
res_cens <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 500,
  UC = 2.0,
  sample_file = "samples_lba_cens.RData"
)

# --- Test 2: LBA recovery with truncated RTs and censoring ------------
res_contam <- run_lba_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 500,
  UC = 2.0,
  UT = 3.0,
  LT = 0.2,
  LC = 0.2,
  sample_file = "samples_lba_cens_trunc.RData"
)

# --- Test 3: LBA recovery with contaminant omissions (e.g., 15%) ------------
res_contam <- run_lba_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 500,
  UC = 2.0,
  sample_file = "samples_lba_contam.RData"
)
