# Applied ridge test on the N-back witness (rdm7s_mass_noA_t0R_UT0_vs2_mt_d42):
# no group update ("baseline") vs the ensemble group update ("ensemble").
#
#   Rscript applied.R <lib> <arm> <seed> <emc0.rds> <outdir>
#
# <emc0.rds> is the unfitted emc from the project's search script
# (RDM_DRYRUN=1 RDM_DRYRUN_SAVE=...), built with the same <lib>. Package
# defaults for preburn/burn/adapt, then a fixed 1,500-iteration sample stage so
# ESS is comparable. 3 chains x 9 cores.
args <- commandArgs(TRUE)
lib <- args[1]; arm <- args[2]; seed <- as.integer(args[3])
emc0 <- args[4]; outdir <- args[5]
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
setwd("/data/work/PM/NirvanaHons_Nback")      # custom trend kernels live here
# Arms: "baseline" = no group update; "ensemble" = ensemble group update.
stopifnot(arm %in% c("baseline", "ensemble"))
options(emc2.ensemble = arm == "ensemble")
emc <- readRDS(emc0)
set.seed(seed)
n_sample <- as.integer(Sys.getenv("NBACK_SAMPLE", "1500"))
timing <- c()
for (stage in c("preburn", "burn", "adapt", "sample")) {
  crit <- EMC2:::get_stop_criteria(stage, if (stage == "sample") list(iter = n_sample) else NULL,
                                   emc[[1]]$type)
  t0 <- proc.time()[["elapsed"]]
  emc <- run_emc(emc, stage = stage, stop_criteria = crit, cores_for_chains = 3L,
                 cores_per_chain = 9L, verbose = FALSE)
  timing[stage] <- proc.time()[["elapsed"]] - t0
  cat(format(Sys.time(), "%T"), arm, "seed", seed, stage, "done", round(timing[stage]), "s\n")
}
emc <- EMC2:::restore_duplicates(emc)
idx <- which(emc[[1]]$samples$stage == "sample")
arr <- function(f) {
  x <- lapply(emc, function(ch) f(ch$samples))
  out <- array(NA_real_, c(length(idx), length(emc), nrow(x[[1]])),
               dimnames = list(NULL, NULL, rownames(x[[1]])))
  for (ch in seq_along(x)) out[, ch, ] <- t(x[[ch]])
  out
}
mu <- arr(function(s) s$theta_mu[, idx])
logsd <- arr(function(s) { v <- apply(s$theta_var[, , idx], 3, diag); 0.5 * log(v) })
ens <- tryCatch(EMC2:::ensemble_diagnostics(emc), error = function(e) NULL)
alpha_ess <- sapply(seq_len(dim(emc[[1]]$samples$alpha)[2]), function(sj)
  sapply(seq_len(dim(emc[[1]]$samples$alpha)[1]), function(p) posterior::ess_bulk(
    sapply(emc, function(ch) ch$samples$alpha[p, sj, idx]))))
saveRDS(list(arm = arm, seed = seed, lib = find.package("EMC2"), timing = timing,
             iters = chain_n(emc)[1, ], mu = mu, logsd = logsd, ens = ens, alpha_ess = alpha_ess,
             subj_ll = sapply(emc, function(ch) colSums(ch$samples$subj_ll[, idx]))),
        file.path(outdir, sprintf("%s-seed%d.rds", arm, seed)))
save(emc, file = file.path(outdir, sprintf("emc-%s-seed%d.RData", arm, seed)))
