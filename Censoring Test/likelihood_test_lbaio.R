## LBAIO likelihood + recovery demo for EMC2
# - Simulates 2-choice LBAIO data with a response deadline (UC) so omissions are encoded as rt=Inf
# - Runs profile plots (R + C likelihood) and a simple parameter-recovery fit
# - Repeats the demo with contaminant omissions (e.g., 15%)
# NB - I couldn't think through a logical situation where truncation could sit with intirnsic omissions logically, so it doesn't support UT<Inf
library(EMC2)
source("Censoring Test/test_likelihood_plotfuns.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- FALSE

# Helper
Rtfun <- function(d) {
  out <- d$rt
  out[is.infinite(out)] <- NA
  out
}

run_lbaio_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                           n_trials = 500, UC = 2.0,
                           cores_for_chains = 3,
                           sample_file = "samples_lbaio.RData") {
  if (UC == Inf) stop("LBAIO demos require a finite UC (deadline) when omissions are present.")

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLBAIO <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBAIO,
    formula = c(
      list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1))
  )

  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBAIO, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(1.0)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(0.4)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.25)

  # Drift means for mismatch vs match accumulators
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.2
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.2

  # Sv for match accumulator
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.2
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.2
  
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
    p_vector, designLBAIO, n_trials = n_trials,
    UC = UC,
    UCresponse = FALSE,      # deadline -> non-responses have unknown choice (R=NA)
    rtContaminantNA = FALSE   # contaminant omissions coded same as censored
  )

  cat("\n--- LBAIO demo ---\n")
  cat("estimate_contaminant =", estimate_contaminant, "\n")
  cat("p_contaminant (true) =", p_contaminant, "\n")
  cat("UC =", UC, " n_trials per cell =", n_trials, "\n\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  profile_plot_test(dat, designLBAIO, p_vector, n_cores = 1, layout = NA, figure_title="R Likelihood") # R likelihood
  profile_plot_test(dat, designLBAIO, p_vector, n_cores = 1, layout = NA, use_c = TRUE, figure_title="C Likelihood") # C likelihood

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBAIO, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, n_post = 20, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designLBAIO, true_pars = p_vector))
}

# --- Test 1: LBAIO recovery without contaminants (pContaminant fixed at 0) ----
res_clean <- run_lbaio_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 2000,
  UC = 2.0,
  sample_file = "samples_lbaio_clean.RData"
)

# --- Test 2: LBAIO recovery with contaminant omissions (e.g., 15%) ------------
res_contam <- run_lbaio_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 500,
  UC = 2.0,
  sample_file = "samples_lbaio_contam15.RData"
)
