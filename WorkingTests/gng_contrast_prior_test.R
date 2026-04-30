rm(list = ls())
library(EMC2)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

UC <- 1.8
LC <- 0.4
UT <- NULL
LT <- NULL
n_trials <- 10000

fit_iter <- 100
max_sample_iter <- 2500
max_tries <- 20
cores_for_chains <- 3
cores_per_chain <- 3

ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "diff"))
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

design_lba_contrast <- design(
  factors = list(subjects = 1, S = c("go", "nogo")),
  Rlevels = c("go", "nogo"),
  matchfun = matchfun,
  functions = list(
    match = function(d) ifelse(d$lM == TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)
  ),
  model = LBA,
  formula = list(v ~ lM, B ~ lR, A ~ 1, t0 ~ 1, sv ~ 1),
  contrasts = list(v = list(lM = ADmat), B = list(lR = ADmat)),
  constants = c(sv = log(1), A = log(0.4))
)

true_pars <- sampled_pars(design_lba_contrast, doMap = FALSE)
if ("v" %in% names(true_pars)) true_pars[["v"]] <- 1.35
if ("v_lMdiff" %in% names(true_pars)) true_pars[["v_lMdiff"]] <- 1.6
if ("B" %in% names(true_pars)) true_pars[["B"]] <- log(0.675)
if ("B_lRdiff" %in% names(true_pars)) true_pars[["B_lRdiff"]] <- -0.23
if ("t0" %in% names(true_pars)) true_pars[["t0"]] <- log(0.2)

cat("True sampled-scale parameters:\n")
print(true_pars)

dat <- make_data(
  true_pars,
  design_lba_contrast,
  n_trials = n_trials,
  TC = list(UC = UC, LC = LC, UT = UT, LT = LT)
)

cat("\nData summary:\n")
print(tapply(is.finite(dat$rt), dat$S, mean))
print(tapply(is.na(dat$R), dat$S, mean))
print(table(S = dat$S, R = addNA(dat$R)))

run_fit <- function(label, prior_list = NULL) {
  emc <- make_emc(dat, design_lba_contrast, type = "single", prior_list = prior_list)
  warns <- character()
  t <- system.time({
    emc <- withCallingHandlers(
      fit(
        emc,
        stop_criteria = list(
          sample = list(
            iter = fit_iter,
            max_sample_iter = max_sample_iter
          ),
          cores_per_chain = cores_per_chain,
          cores_for_chains = cores_for_chains
        ),
        max_tries = max_tries
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })

  cn <- EMC2:::chain_n(emc)
  cat("\n===", label, "===\n")
  cat("elapsed_sec:", unname(t["elapsed"]), "\n")
  print(cn)
  cat("last_stage:", EMC2:::get_last_stage(emc), "\n")
  cat("max_tries_warning:", any(grepl("Max tries reached", warns, fixed = TRUE)), "\n")
  cat("recovery:\n")
  print(recovery(emc, true_pars = true_pars))
  invisible(emc)
}

# 1) default prior
emc_default <- run_fit("default prior")

# 2) informative prior centered at true values
psd <- rep(0.5, length(true_pars))
names(psd) <- names(true_pars)
if ("t0" %in% names(psd)) psd["t0"] <- 0.08
psd[grepl("^v", names(psd)) & grepl("*diff", names(psd))] <- 0.5
psd[grepl("^B", names(psd)) & grepl("*diff", names(psd))] <- 0.5
prior_info <- prior(design_lba_contrast, type = "single", pmean = true_pars, psd = psd)
emc_info <- run_fit("informative prior", prior_list = prior_info)

cat("\nDone.\n")
