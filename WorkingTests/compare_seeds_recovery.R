## Compare recovery stats across seeds for UC=Inf (no censoring/truncation) LBA fit
## Checks whether seed=123 is anomalously bad vs other seeds

rm(list = ls())
library(EMC2)

# Reconstruct the design and true_pars exactly as in run_lba_demo()
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designLBA <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = LBA(posdrift = TRUE),
  formula = list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM),
  constants = c(A = log(0.4))
)

p_vector <- sampled_pars(designLBA, doMap = FALSE)
p_vector[["B"]]          <- log(1.2)
p_vector[["A"]]          <- log(0.4)
p_vector[["t0"]]         <- log(0.15)
p_vector[["v_lMFALSE"]]  <- 0.4
p_vector[["v_lMTRUE"]]   <- 1.2
p_vector[["sv_lMFALSE"]] <- log(1.2)
p_vector[["sv_lMTRUE"]]  <- log(1)

# A is a constant in the design — exclude from true_pars passed to recovery()
free_pars <- p_vector[names(p_vector) != "A"]

cat("True parameters (free only):\n")
print(free_pars)
cat("\n")

# Helper to load an emc object from a saved RData file and report recovery
report_recovery <- function(rdata_path, label, true_pars) {
  if (!file.exists(rdata_path)) {
    cat(sprintf("[%s] File not found: %s\n\n", label, rdata_path))
    return(invisible(NULL))
  }
  env <- new.env()
  load(rdata_path, envir = env)
  # fit() saves the emc object; find it
  obj_names <- ls(env)
  emc_obj <- NULL
  for (nm in obj_names) {
    obj <- get(nm, envir = env)
    if (inherits(obj, "emc")) { emc_obj <- obj; break }
    # Sometimes stored as a list with $emc
    if (is.list(obj) && !is.null(obj$emc) && inherits(obj$emc, "emc")) {
      emc_obj <- obj$emc; break
    }
  }
  if (is.null(emc_obj)) {
    # Try direct assignment: fit() often saves 'emc' as the object name
    if (exists("emc", envir = env)) emc_obj <- get("emc", envir = env)
  }
  if (is.null(emc_obj)) {
    cat(sprintf("[%s] Could not find emc object in %s (found: %s)\n\n",
                label, rdata_path, paste(obj_names, collapse=", ")))
    return(invisible(NULL))
  }
  cat(sprintf("=== %s ===\n", label))
  tryCatch({
    rec <- recovery(emc_obj, true_pars = true_pars)
    print(rec)
  }, error = function(e) cat("recovery() error:", conditionMessage(e), "\n"))
  cat("\n")
  invisible(emc_obj)
}

wd <- "/data/work/EMC2_dev_oo"

# Seed 123, UC=Inf — the original failing case
report_recovery(file.path(wd, "samples_ss_no_cens_trunc.RData"),
                "seed=123, UC=Inf (original label)", free_pars)

report_recovery(file.path(wd, "samples_ss_a_inf_seed123.RData"),
                "seed=123, UC=Inf (a_inf label)", free_pars)

# Seed 123, UC=100 — the comparison that "worked"
report_recovery(file.path(wd, "samples_ss_no_cens_trunc_uc100.RData"),
                "seed=123, UC=100 (original label)", free_pars)

report_recovery(file.path(wd, "samples_ss_b_100_seed123.RData"),
                "seed=123, UC=100 (b_100 label)", free_pars)

# Seed 456, UC=Inf — alternative seed
report_recovery(file.path(wd, "samples_ss_seed456.RData"),
                "seed=456, UC=Inf", free_pars)

# Run additional seeds if needed (will generate data + fit if not cached)
run_seed_recovery <- function(seed, label, wd) {
  fname <- file.path(wd, sprintf("samples_ss_seed%d_ucinf.RData", seed))
  if (file.exists(fname)) {
    report_recovery(fname, sprintf("seed=%d, UC=Inf (cached)", seed), p_vector)
    return(invisible(NULL))
  }

  cat(sprintf("=== seed=%d, UC=Inf (running fit) ===\n", seed))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  dat <- make_data(
    p_vector, designLBA,
    n_trials = 10000,
    TC = list(UC = Inf, UT = Inf, LC = 0, LT = 0,
              LCresponse = FALSE, UCresponse = FALSE,
              LCdirection = TRUE, UCdirection = TRUE)
  )
  emc <- make_emc(dat, designLBA, type = "single", rt_resolution = NULL)
  emc <- fit(emc, fileName = fname,
             stop_criteria = list(
               sample = list(iter = 1000, max_gd = 1.10,
                             max_flat_loc = 0.5,
                             flat_selection = c("alpha", "subj_ll"),
                             flat_p1 = 1/3, flat_p2 = 1/3,
                             max_sample_iter = 5000)
             ), cores_per_chain = 3, cores_for_chains = 3, max_tries = 30)
  tryCatch({
    rec <- recovery(emc, true_pars = free_pars)
    print(rec)
  }, error = function(e) cat("recovery() error:", conditionMessage(e), "\n"))
  cat("\n")
}

# Try seeds 456 and 789 with UC=Inf
run_seed_recovery(456, "seed456_ucinf", wd)
run_seed_recovery(789, "seed789_ucinf", wd)
