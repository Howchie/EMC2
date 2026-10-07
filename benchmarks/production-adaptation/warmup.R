# Burn + adapt from a saved preburn state with the default stopping rules, so
# every production arm of one (case, seed) starts from the same state.
#   CASE=wide15 SEED=1 EMC_LIB=<lib> Rscript benchmarks/production-adaptation/warmup.R
# Preburn states: benchmarks/mode-start/results/preburn_<case>.rds (setup.R).
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
cat("EMC2 from", find.package("EMC2"), "\n")
src <- Sys.getenv("PREBURN_DIR", "benchmarks/mode-start/results")
out <- Sys.getenv("PADAPT_OUT", "benchmarks/production-adaptation/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
case <- Sys.getenv("CASE")
seed <- as.integer(Sys.getenv("SEED", "1"))
cores <- as.integer(Sys.getenv("CORES_PER_CHAIN", "2"))
emc <- readRDS(file.path(src, paste0("preburn_", case, ".rds")))
# INDEP=1 runs the chains fully independently from burn onwards.
if (Sys.getenv("INDEP", "0") == "1") options(emc2.independent_chains = TRUE)
# GROUP_MOVE=ensemble runs the ensemble group update without the shear.
options(emc2.group_move = Sys.getenv("GROUP_MOVE", "shear"))
set.seed(1000L * seed + 7L)
# Seed 1 keeps the chain streams stored with the preburn; other seeds replace them.
if (seed != 1L) for (j in seq_along(emc)) emc[[j]]$rng <- NULL

cpu <- function() { p <- proc.time(); sum(p[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) }
stage_run <- function(emc, stage, stop_criteria) {
  c0 <- cpu(); t0 <- proc.time()[["elapsed"]]
  emc <- suppressWarnings(run_emc(emc, stage = stage, stop_criteria = stop_criteria,
    cores_for_chains = 3, cores_per_chain = cores, step_size = 100,
    max_tries = 20, verbose = FALSE))
  list(emc = emc, cpu = cpu() - c0, elapsed = proc.time()[["elapsed"]] - t0)
}
# The defaults fit() uses for these stages.
# BURN_RULE=default is fit()'s mean R-hat < 1.1 on each of alpha and mu;
# BURN_RULE=strict keeps that for alpha but needs max R-hat < 1.2 on mu.
burn_rule <- Sys.getenv("BURN_RULE", "default")
burn_sc <- switch(burn_rule,
  default = list(mean_gd = 1.1),
  strict = list(mean_gd = c(alpha = 1.1), max_gd = c(mu = 1.2)),
  stop("unknown BURN_RULE"))
b <- stage_run(emc, "burn", c(burn_sc, list(omit_mpsrf = TRUE, selection = c("alpha", "mu"))))
a <- stage_run(b$emc, "adapt", list(min_unique = 150))
emc <- a$emc
# Burn's final R-hat check (mean per selection under the default rule).
blog <- Filter(function(e) identical(e$stage, "burn"), emc[[1]]$convergence_log)
bst <- blog[[length(blog)]]$gd$stats
info <- data.frame(case = case, seed = seed,
  burn_converged = !is.null(bst) && all(bst$value < bst$threshold),
  burn_rule = burn_rule,
  burn_rhat_alpha = bst$value[bst$selection == "alpha"],
  burn_rhat_mu = bst$value[bst$selection == "mu"],
  burn_iter = chain_n(emc)[1, "burn"], adapt_iter = chain_n(emc)[1, "adapt"],
  burn_cpu = b$cpu, adapt_cpu = a$cpu, burn_elapsed = b$elapsed, adapt_elapsed = a$elapsed)
print(info)
attr(emc, "warmup_info") <- info
saveRDS(emc, file.path(out, sprintf("warmup_%s_%d.rds", case, seed)))
