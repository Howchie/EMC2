# Burn from a saved post-preburn state; the build under test is EMC_LIB.
#   CASE=freeA15 ARM=<label> SEED=1 [STOP=fixed] [EMC2_OPTS=...] Rscript benchmarks/mode-start/burn_arms.R
# STOP=fixed runs 500 burn iterations; otherwise the default burn stopping rule
# (mean R-hat < 1.1 on alpha and mu, max_tries 10) decides.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
# EMC2_OPTS: extra options() for this run, e.g. "emc2.sampler_profile = TRUE".
if (nzchar(Sys.getenv("EMC2_OPTS")))
  eval(parse(text = paste0("options(", Sys.getenv("EMC2_OPTS"), ")")))
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
case <- Sys.getenv("CASE"); arm <- Sys.getenv("ARM", "default")
seed <- as.integer(Sys.getenv("SEED", "1"))
fixed <- identical(Sys.getenv("STOP"), "fixed")
emc <- readRDS(file.path(out, paste0("preburn_", case, ".rds")))
set.seed(1000L * seed + 7L)
# The chains draw from random streams stored in the preburn object; seed 1 keeps
# them, other seeds replace them with streams from set.seed().
if (seed != 1L) for (j in seq_along(emc)) emc[[j]]$rng <- NULL
cores <- as.integer(Sys.getenv("CORES_PER_CHAIN", "2"))
stop_criteria <- if (fixed) list(iter = 500) else
  list(mean_gd = 1.1, omit_mpsrf = TRUE, selection = c("alpha", "mu"))
t <- system.time(emc <- suppressWarnings(run_emc(emc, stage = "burn",
  stop_criteria = stop_criteria, cores_for_chains = 3, cores_per_chain = cores,
  max_tries = if (fixed) 50 else 10, trim = FALSE, verbose = FALSE)))
n_burn <- chain_n(emc)[1L, "burn"]
rh <- function(sel, stat) tryCatch(suppressWarnings(as.numeric(gd_summary(emc, selection = sel,
  stage = "burn", filter = floor(n_burn / 2), stat = stat, stat_only = TRUE, omit_mpsrf = TRUE))),
  error = function(e) NA_real_)
res <- data.frame(case = case, arm = arm, seed = seed, stop = if (fixed) "fixed" else "default",
  burn_iter = n_burn, elapsed = t[["elapsed"]], cpu = t[["user.self"]] + t[["user.child"]],
  mean_rhat_alpha = rh("alpha", "mean"), mean_rhat_mu = rh("mu", "mean"),
  max_rhat_alpha = rh("alpha", "max"), max_rhat_mu = rh("mu", "max"),
  eff_ready = isTRUE(emc[[1]]$burn_eff_ready))
# Local random-walk step size per chain x subject.
eps <- unlist(lapply(emc, function(ch) vapply(attr(ch$samples, "pm_settings"),
  function(p) p[[1L]]$epsilon[1L], numeric(1))))
res$eps_min <- min(eps); res$eps_median <- stats::median(eps); res$eps_max <- max(eps)
# Stuck subjects: per subject-parameter R-hat over the second half of burn, and
# chains that moved on under 30% of the last 100 iterations.
idx <- emc[[1L]]$samples$idx
half <- seq.int(idx - floor(n_burn / 2) + 1L, idx)
a <- simplify2array(lapply(emc, function(ch) ch$samples$alpha[, , half, drop = FALSE]))
n <- length(half)
rhat <- apply(a, c(1, 2), function(m) {
  m <- matrix(m, n)
  W <- mean(apply(m, 2, stats::var)); B <- n * stats::var(colMeans(m))
  sqrt(((n - 1) / n * W + B / n) / W)
})
last <- seq.int(idx - 99L, idx)
moved <- sapply(emc, function(ch) apply(ch$samples$alpha[, , last, drop = FALSE], 2,
  function(m) mean(rowSums(abs(diff(t(m)))) > 0)))
res$n_pairs_rhat_1.1 <- sum(rhat > 1.1, na.rm = TRUE)
res$n_pairs <- length(rhat)
res$n_subj_rhat_1.1 <- sum(apply(rhat, 2, max, na.rm = TRUE) > 1.1)
res$n_subj_stuck <- sum(apply(moved, 1, min) < 0.3)
res$n_subj <- ncol(rhat)
print(res)
saveRDS(res, file.path(out, sprintf("burnarm_%s_%s_%s_%d.rds", case, arm, res$stop, seed)))
