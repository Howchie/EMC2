rm(list = ls())
library(EMC2)

# Quick diagnostic comparison for GNG identifiability under censoring.
# Goal: contrast race-style LBA GNG against canonical DDM GNG with the same TC bounds.
#
# Why this is useful:
# - LBA with B~0+lR can trade off threshold split vs drift split when nogo is mostly omissions.
# - DDM has one boundary (a) + start bias (Z), which tends to be better constrained.

UC <- as.numeric(Sys.getenv("GNG_UC", "1.8"))
LC <- as.numeric(Sys.getenv("GNG_LC", "0.4"))
UT_env <- Sys.getenv("GNG_UT", "")
LT_env <- Sys.getenv("GNG_LT", "")
UT <- if (nzchar(UT_env)) as.numeric(UT_env) else NULL
LT <- if (nzchar(LT_env)) as.numeric(LT_env) else NULL
n_trials <- as.integer(Sys.getenv("GNG_N_TRIALS", "10000"))

fit_iter <- as.integer(Sys.getenv("GNG_FIT_ITER", "1000"))
max_sample_iter <- as.integer(Sys.getenv("GNG_MAX_SAMPLE_ITER", "5000"))
max_tries <- as.integer(Sys.getenv("GNG_MAX_TRIES", "30"))
cores_for_chains <- as.integer(Sys.getenv("GNG_CORES_FOR_CHAINS", "3"))
cores_per_chain <- as.integer(Sys.getenv("GNG_CORES_PER_CHAIN", "3"))

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

summarise_data <- function(dat, tag) {
  cat("\n====", tag, "data summary ====\n")
  cat("N:", nrow(dat), "\n")
  cat("Finite RT proportion:", mean(is.finite(dat$rt)), "\n")
  cat("rt==-Inf proportion:", mean(dat$rt == -Inf, na.rm = TRUE), "\n")
  cat("rt==Inf proportion:", mean(dat$rt == Inf, na.rm = TRUE), "\n")
  cat("Observed response proportion:", mean(!is.na(dat$R)), "\n")
  cat("By stimulus finite RT proportion:\n")
  print(tapply(is.finite(dat$rt), dat$S, mean))
  cat("By stimulus omission proportion:\n")
  print(tapply(is.na(dat$R), dat$S, mean))
  cat("Counts by S x R:\n")
  print(table(S = dat$S, R = addNA(dat$R)))
}

fit_with_summary <- function(dat, design_obj, tag) {
  emc <- make_emc(dat, design_obj, type = "single")
  fit_warnings <- character()
  t <- system.time({
    emc <- withCallingHandlers(
      fit(
        emc,
        stop_criteria = list(
          sample = list(
            iter = fit_iter,
            max_gd = 1.10,
            max_sample_iter = max_sample_iter
          ),
          cores_per_chain = cores_per_chain,
          cores_for_chains = cores_for_chains
        ),
        max_tries = max_tries
      ),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })

  cn <- EMC2:::chain_n(emc)
  gd_obj <- tryCatch(EMC2:::gd(emc), error = function(e) NULL)

  cat("\n====", tag, "fit summary ====\n")
  cat("Elapsed (sec):", unname(t["elapsed"]), "\n")
  cat("chain_n:\n")
  print(cn)
  cat("last stage:", EMC2:::get_last_stage(emc), "\n")
  cat("max_tries warning:", any(grepl("Max tries reached", fit_warnings, fixed = TRUE)), "\n")
  if (!is.null(gd_obj)) {
    cat("gd:\n")
    print(gd_obj)
  }

  invisible(list(emc = emc, warnings = unique(fit_warnings)))
}

## 1) LBA GNG with split thresholds (hard case)
matchfun_lba <- function(d) as.numeric(d$S) == as.numeric(d$lR)
design_lba_split <- design(
  factors = list(subjects = 1, S = c("go", "nogo", "nogo", "nogo"), L = c("L", "M")),
  Rlevels = c("go", "nogo"),
  matchfun = matchfun_lba,
  functions = list(
    match = function(d) ifelse(d$lM == TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)
  ),
  model = LBA,
  UC = 3,
  formula = list(B ~ 0 + lR, v ~ 0 + mismatch + S:match:L, A ~ 1, t0 ~ 1, sv ~ 1),
  constants = c(sv = log(1), A = log(0.4))
)

p_lba_split <- sampled_pars(design_lba_split, doMap = FALSE)
if ("B_lRgo" %in% names(p_lba_split)) p_lba_split[["B_lRgo"]] <- log(.75)
if ("B_lRnogo" %in% names(p_lba_split)) p_lba_split[["B_lRnogo"]] <- log(.6)
if ("t0" %in% names(p_lba_split)) p_lba_split[["t0"]] <- log(0.2)
if ("v_mismatch" %in% names(p_lba_split)) p_lba_split[["v_mismatch"]] <- .5
if ("v_Sgo:match:LL" %in% names(p_lba_split)) p_lba_split[["v_Sgo:match:LL"]] <- 2.1
if ("v_Snogo:match:LL" %in% names(p_lba_split)) p_lba_split[["v_Snogo:match:LL"]] <- 2.3
if ("v_Sgo:match:LM" %in% names(p_lba_split)) p_lba_split[["v_Sgo:match:LM"]] <- 1.65
if ("v_Snogo:match:LM" %in% names(p_lba_split)) p_lba_split[["v_Snogo:match:LM"]] <- 1.95

set.seed(123)
dat_lba_split <- make_data(
  p_lba_split,
  design_lba_split,
  n_trials = n_trials,
  TC = list(UC = UC, LC = LC, UT = UT, LT = LT)
)
summarise_data(dat_lba_split, "LBA split-B")
fit_lba_split <- fit_with_summary(dat_lba_split, design_lba_split, "LBA split-B")

## 2) LBA GNG with shared threshold (usually easier)
design_lba_shared <- design(
  factors = list(subjects = 1, S = c("go", "nogo", "nogo", "nogo"), L = c("L", "M")),
  Rlevels = c("go", "nogo"),
  matchfun = matchfun_lba,
  functions = list(
    match = function(d) ifelse(d$lM == TRUE, 1, 0),
    mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)
  ),
  model = LBA,
  UC = 3,
  formula = list(B ~ 1, v ~ 0 + mismatch + S:match:L, A ~ 1, t0 ~ 1, sv ~ 1),
  constants = c(sv = log(1), A = log(0.4))
)

p_lba_shared <- sampled_pars(design_lba_shared, doMap = FALSE)
if ("B" %in% names(p_lba_shared)) p_lba_shared[["B"]] <- log(.68)
if ("t0" %in% names(p_lba_shared)) p_lba_shared[["t0"]] <- log(0.2)
if ("v_mismatch" %in% names(p_lba_shared)) p_lba_shared[["v_mismatch"]] <- .5
if ("v_Sgo:match:LL" %in% names(p_lba_shared)) p_lba_shared[["v_Sgo:match:LL"]] <- 2.1
if ("v_Snogo:match:LL" %in% names(p_lba_shared)) p_lba_shared[["v_Snogo:match:LL"]] <- 2.3
if ("v_Sgo:match:LM" %in% names(p_lba_shared)) p_lba_shared[["v_Sgo:match:LM"]] <- 1.65
if ("v_Snogo:match:LM" %in% names(p_lba_shared)) p_lba_shared[["v_Snogo:match:LM"]] <- 1.95

set.seed(123)
dat_lba_shared <- make_data(
  p_lba_shared,
  design_lba_shared,
  n_trials = n_trials,
  TC = list(UC = UC, LC = LC, UT = UT, LT = LT)
)
summarise_data(dat_lba_shared, "LBA shared-B")
fit_lba_shared <- fit_with_summary(dat_lba_shared, design_lba_shared, "LBA shared-B")

## 3) Canonical DDM GNG (C++ path: DDM + Rlevels includes 'nogo')
design_ddm <- design(
  factors = list(subjects = 1, S = c("go", "nogo")),
  Rlevels = c("go", "nogo"),
  model = DDM,
  formula = list(v ~ 0 + S, a ~ 1, t0 ~ 1, Z ~ 1),
  constants = c(s = log(1), sv = log(0), SZ = qnorm(0), st0 = log(0))
)

p_ddm <- sampled_pars(design_ddm, doMap = FALSE)
if ("v_Sgo" %in% names(p_ddm)) p_ddm[["v_Sgo"]] <- -1.5
if ("v_Snogo" %in% names(p_ddm)) p_ddm[["v_Snogo"]] <- 1.5
if ("a" %in% names(p_ddm)) p_ddm[["a"]] <- log(2.0)
if ("t0" %in% names(p_ddm)) p_ddm[["t0"]] <- log(0.2)
if ("Z" %in% names(p_ddm)) p_ddm[["Z"]] <- qnorm(0.5)

set.seed(123)
dat_ddm <- make_data(
  p_ddm,
  design_ddm,
  n_trials = n_trials,
  TC = list(UC = UC, LC = LC, UT = UT, LT = LT)
)
summarise_data(dat_ddm, "DDM")
fit_ddm <- fit_with_summary(dat_ddm, design_ddm, "DDM")

cat("\nDone. Compare the three fit summaries above.\n")
