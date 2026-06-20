# WorkingTests/test_SBC.R
# --------------------------------------------------------------------------
# Demonstrates recovering / resuming an interrupted SBC run for a SINGLE
# (non-hierarchical) fit, using the LNR design from the example batch script.
#
# Scenarios covered:
#   (A) A run is interrupted (node reboot / deliberate abort). You still have
#       this script (design + prior) and the intermediate .rds files. Use
#       recover_sbc() to assemble and inspect the replicates completed so far.
#   (B) Worst case: only the per-replicate rep_<i>.rds files survive (the
#       prior_samples.rds and the .RData are gone). For a single fit the
#       results are still fully recoverable for inspection (ranks/coverage/bias
#       are baked into the rep files), but the run can no longer be resumed.
#   (C) Resume to the target: re-issue the original run_sbc() call with the same
#       fileName; it runs only the missing replicates.
#
# An interruption cannot be forced from within a script, so scenario (A) is
# emulated by running the per-replicate worker for a SUBSET of replicates,
# leaving a genuine partial _temp directory behind.
#
# Settings are scaled down for a quick demonstration. Set FAST <- FALSE (and/or
# edit replicates/trials/cores) to approach the full batch
# (replicates = 500, trials = 100, cores_per_chain = 3).
# --------------------------------------------------------------------------

rm(list = ls())
library(EMC2)

FAST       <- TRUE
replicates <- if (FAST) 12 else 500
n_partial  <- if (FAST) 5  else 200    # replicates "completed" before the crash
trials     <- 100
cores      <- 3

# --- design + prior (from the batch script) --------------------------------

matchfun <- function(d) d$S == d$lR
UC <- UT <- function(d) quantile(d$rt, .9)
LC <- LT <- function(d) quantile(d$rt, .1)

design_LNR <- design(factors = list(subjects = 1, S = c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun, TC = list(rt_resolution = 1/60, LT = LT, UT = UT),
                     formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                     model = LNR)

prior_LNR <- prior(design_LNR, type = "single",
                   pmean = c(-.7, -.5, log(1), log(.2)),
                   psd   = c(.2, .1, .1, .05))

fileName <- "SBC_TLNR.RData"
temp_dir <- paste0(tools::file_path_sans_ext(fileName), "_temp")

# start clean
unlink(temp_dir, recursive = TRUE)
unlink(fileName)

# --------------------------------------------------------------------------
# Emulate an interrupted run: draw the prior, lay down the temp directory the
# way run_sbc() would, then run the worker for the first n_partial replicates.
# --------------------------------------------------------------------------
dir.create(temp_dir, showWarnings = FALSE)
prior_alpha <- parameters(prior_LNR, N = replicates, selection = "alpha")
saveRDS(prior_alpha, file.path(temp_dir, "prior_samples.rds"))
save(prior_alpha, file = fileName)

# dots mirror what SBC_single() builds internally
dots <- EMC2:::add_defaults(
  list(cores_per_chain = 1),
  max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
  stop_criteria = list(min_es = 100, max_gd = 1.1, selection = c("alpha", "mu", "Sigma")),
  cores_per_chain = 1
)
dots$verbose <- FALSE

cat("\n== Emulating a crash after", n_partial, "of", replicates, "replicates ==\n")
for (rep in seq_len(n_partial))
  EMC2:::run_SBC_subject(rep, design_LNR, prior_alpha, trials, prior_LNR, dots, temp_dir)

# --------------------------------------------------------------------------
# Scenario (A): inspect the partial results without running anything more.
# --------------------------------------------------------------------------
cat("\n== Scenario A: recover_sbc() on the partial run ==\n")
# Returns what run_sbc() would have returned at n_partial replicates. The first
# argument is the `temp_dir` holding the rep files; the second is the `design`.
# `fileName` is optional (as in run_sbc): supply it to also save the result in
# run_sbc()'s format. Here we save it under a separate name so the original run
# files are left untouched.
recovered <- "TLNR_partial_recovered.RData"
SBC_partial <- recover_sbc(temp_dir, design_LNR, fileName = recovered,
                           prior_in = prior_LNR)
cat("recovered_reps:", attr(SBC_partial, "recovered_reps"), "\n")
stopifnot(nrow(SBC_partial$rank$alpha) == n_partial)
# the saved file holds SBC (+ prior_alpha), exactly like a normal single run
local({ e <- new.env(); load(recovered, envir = e)
        cat("Saved objects:", paste(ls(e), collapse = ", "), "\n") })

pdf("TLNR_partial.pdf")
plot_sbc_ecdf(SBC_partial)
dev.off()
cat("Partial ECDF written to TLNR_partial.pdf\n")

# --------------------------------------------------------------------------
# Scenario (B): worst case -- only rep_<i>.rds survive. Move the prior + .RData
# aside and confirm a single run is still recoverable for inspection.
# --------------------------------------------------------------------------
cat("\n== Scenario B: worst case (only rep_<i>.rds files) ==\n")
ps_bak <- file.path(temp_dir, "prior_samples.rds")
file.rename(ps_bak, paste0(ps_bak, ".bak"))
fn_bak <- paste0(fileName, ".bak")
file.rename(fileName, fn_bak)

# the prior is gone; recover_sbc still finds rep_<i>.rds in temp_dir and recovers
# for inspection (prints "Inspect-only").
SBC_worst <- recover_sbc(temp_dir, design_LNR, fileName = "TLNR_worst_recovered.RData",
                         prior_in = prior_LNR)
stopifnot(nrow(SBC_worst$rank$alpha) == n_partial)
stopifnot(identical(unname(SBC_worst$rank$alpha), unname(SBC_partial$rank$alpha)))
cat("Worst-case recovery matches the full-prior recovery (ranks identical).\n")

# restore so we can resume
file.rename(paste0(ps_bak, ".bak"), ps_bak)
file.rename(fn_bak, fileName)

# --------------------------------------------------------------------------
# Scenario (C): resume to the target by re-issuing the original run_sbc() call.
# It reads the completed replicates from temp_dir and runs only the rest.
# --------------------------------------------------------------------------
cat("\n== Scenario C: resume to", replicates, "replicates ==\n")
SBC_TLNR <- run_sbc(design_LNR, prior_LNR, replicates = replicates, trials = trials,
                    fileName = fileName, cores_per_chain = cores)
stopifnot(nrow(SBC_TLNR$rank$alpha) == replicates)
cat("Resumed run completed with", nrow(SBC_TLNR$rank$alpha), "replicates.\n")

# tail of the batch script
save(SBC_TLNR, file = "TLNR.RData")
pdf("TLNR.pdf")
plot_sbc_ecdf(SBC_TLNR)
dev.off()
cat("\nDone. Full ECDF written to TLNR.pdf\n")
