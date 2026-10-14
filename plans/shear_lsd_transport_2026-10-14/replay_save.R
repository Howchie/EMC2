# Replay the continuation's forced production rebuild on the lbaE checkpoint
# and report each chain's final plan.
args <- commandArgs(TRUE); ckpt <- args[1]; seed <- as.integer(args[2])
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
patch <- Sys.getenv("PATCH"); if (nzchar(patch)) { e <- new.env(); sys.source(patch, e); for (n in ls(e, all.names = TRUE)) { f <- get(n, e); if (is.function(f)) environment(f) <- ns; if (exists(n, ns, inherits = FALSE)) { unlockBinding(n, ns); assign(n, f, ns); lockBinding(n, ns) } else assign(n, f, globalenv()) } }
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(ckpt)), quiet = TRUE)
class(emc) <- "emc"
base <- local({ set.seed(seed, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) {
  base <- parallel::nextRNGStream(base)
  emc[[j]]$rng <- ns$.emc_chain_rng(base, emc[[j]]$n_subjects)
  emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions
}
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample",
  stop_criteria = list(sample = list(iter = 1)), step_size = 1L, max_tries = 1L,
  cores_for_chains = length(emc), cores_per_chain = 1, trim = FALSE,
  particle_factor = 50, verbose = FALSE, verboseProgress = FALSE)))
saveRDS(emc, Sys.getenv("OUT")); quit(save="no")
for (c in seq_along(emc)) {
  h <- emc[[c]]$alloc$history; n <- length(h)
  for (b in max(1, n - 2):n) { e <- h[[b]]; qs <- e$particle_shared; if (is.null(qs)) qs <- e$particle_rate
    cat(sprintf("chain %d b%d it %d pass2=%s choice %s [%s] regime %s: own q %.2f shared q %.2f; r>0 %d, mean r %.2f; starved (q<1, r=0) %d\n",
      c, b, e$iteration, e$particle_fixed, if (is.null(e$choice)) "-" else e$choice, if (is.null(e$options)) "" else paste(sprintf("%s=%.4g", names(e$options), e$options), collapse=" "), e$regime, mean(e$particle_rate), mean(qs), sum(e$rate > 0), mean(e$rate), sum(qs < 1 & e$rate == 0))) }
}
