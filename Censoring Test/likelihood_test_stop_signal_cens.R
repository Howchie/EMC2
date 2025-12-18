## Stop-signal censoring demo for EMC2
# - Simulates stop-signal data (optionally with a stop-triggered *response* accumulator)
# - Encodes "no response by deadline" as rt = Inf with finite UC
# - Uses the SS censored likelihood (EMC2:::log_likelihood_race_ss_cens_trunc)
# - Runs profile plots and (optionally) a parameter recovery fit
#
# Key idea:
# - The SS model always includes an unobserved STOP process on stop trials (finite SSD).
# - `include_stop_triggered = FALSE` does NOT remove stop signals; it removes the extra
#   stop-triggered *response* accumulator (labelled "st" here). Without that accumulator,
#   when the stop process wins, the observed outcome is "no response" (R=NA, rt=Inf).

library(EMC2)
source("Censoring Test/test_likelihood_plotfuns.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

include_stop_triggered <- FALSE
normalise_trunc <- FALSE
RUN_PROFILE_PLOTS <- TRUE
RUN_FITS <- FALSE

if (include_stop_triggered && normalise_trunc) {
  stop("normalise_trunc=TRUE is not supported with stop-triggered accumulators in the current SS censored likelihood.")
}

# Wrapper model that swaps in the censored SS likelihood.
SSexG_cens <- function(normalise_trunc = FALSE) {
  m <- SSexG()
  m$c_name <- "SSexG"
  m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    EMC2:::log_likelihood_race_ss_cens_trunc(
      pars, dadm, model, min_ll = min_ll, normalise_trunc = normalise_trunc
    )
  }
  m
}

# Stop-triggered indicator per accumulator (lI: 1=ST, 2=GO)
lIfun <- function(d) {
  if (!include_stop_triggered) {
    return(factor(rep(2, nrow(d)), levels = 1:2))
  }
  factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
}

# SSD assignment: finite SSD => stop trial; Inf => go trial
mySSD_function <- function(d) {
  # 50% stop trials split across two SSD values; remaining trials are go trials (SSD=Inf)
  SSD_function(d, SSD = c(.15, .25), pSSD = c(.25, .25))
}

# Define the accumulator set (Rlevels). Include an extra response level for ST if requested.
Rlevels <- if (include_stop_triggered) c("left", "right", "st") else c("left", "right")

matchfun <- function(d) as.character(d$S) == as.character(d$lR)

designSS <- design(
  model = function() SSexG_cens(normalise_trunc = normalise_trunc),
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = Rlevels,
  matchfun = matchfun,
  functions = list(lI = lIfun, SSD = mySSD_function),
  # Allow accumulator-specific finishing-time parameters for go/ST racers
  formula = list(mu ~ 0 + lR, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)

# Set simulation parameters (on transformed scale expected by p_types)
p_vector <- sampled_pars(designSS, doMap = FALSE)

# Reasonable defaults; adjust as needed.
for (nm in names(p_vector)) {
  if (grepl("^mu_", nm)) p_vector[[nm]] <- log(.55)
}
if ("mu_lRst" %in% names(p_vector)) p_vector[["mu_lRst"]] <- log(.35)
if ("muS" %in% names(p_vector)) p_vector[["muS"]] <- log(.30)
if ("sigma" %in% names(p_vector)) p_vector[["sigma"]] <- log(.12)
if ("tau" %in% names(p_vector)) p_vector[["tau"]] <- log(.10)
if ("sigmaS" %in% names(p_vector)) p_vector[["sigmaS"]] <- log(.08)
if ("tauS" %in% names(p_vector)) p_vector[["tauS"]] <- log(.06)
if ("gf" %in% names(p_vector)) p_vector[["gf"]] <- qnorm(.10)
if ("tf" %in% names(p_vector)) p_vector[["tf"]] <- qnorm(.10)

# Deadline (UC) + optional truncation bounds (LT/UT)
UC <- 1.5
LT <- if (normalise_trunc) 0.15 else 0
UT <- if (normalise_trunc) 2 else Inf

dat <- make_data(
  p_vector, designSS, n_trials = 500,
  UC = UC, UCresponse = FALSE,
  LT = LT, UT = UT,
  return_Ffunctions = TRUE
)

# Convert SS "no response" (rt==NA from stop success / go-failure / dual failure) into deadline-censored NR.
dat$rt[is.na(dat$rt)] <- Inf

# Sanity checks / quick summary
summarise_ss_data <- function(d) {
  is_stop <- is.finite(d$SSD)
  is_nr <- is.na(d$R) | is.infinite(d$rt)
  resp <- as.character(d$R)
  resp[is.na(resp)] <- "NR"

  cat("include_stop_triggered =", include_stop_triggered, "\n")
  cat("Rlevels =", paste(Rlevels, collapse = ", "), "\n\n")

  cat("Counts by trial type (finite SSD => stop):\n")
  print(table(trial_type = ifelse(is_stop, "stop", "go")))

  cat("\nCounts by response (NR = no response):\n")
  print(table(response = resp))

  cat("\nCounts by trial type × response:\n")
  print(table(trial_type = ifelse(is_stop, "stop", "go"), response = resp))

  cat("\nNon-response rate by trial type:\n")
  print(tapply(is_nr, ifelse(is_stop, "stop", "go"), mean))
}

summarise_ss_data(dat)

if (isTRUE(RUN_PROFILE_PLOTS)) {
  # Profile plots using the (R) censored SS likelihood embedded in designSS
  profile_plot_test(dat, designSS, p_vector, n_cores = 1, layout = c(3, 3), use_c = FALSE, figure_title = "R Likelihood")
}

if (isTRUE(RUN_FITS)) {
  # Parameter recovery fit
  emc <- make_emc(dat, designSS, type = "single")
  emc <- fit(emc, cores_for_chains = 3, fileName = "samples_ss_cens.RData")

  # Inspect recovery (posterior predictive is optional)
  post_predict <- predict(emc, n_post = 50)
  plot_pars(emc, post_predict = post_predict, n_post = 20, true_pars = p_vector)
}
