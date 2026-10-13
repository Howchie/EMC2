# Calibration run for the allocator audit: continue a sample-stage checkpoint
# for ITERS production transitions with rebuilds at the offsets REBUILDS
# (after a forced rebuild at the start), recording every rebuild's inputs,
# measurements and predictions (trace.R) and every draw.
#   EMC_LIB=<lib> R_LIBS=<lib> Rscript benchmarks/allocator/calibrate.R <ckpt.rds> <out_dir> <seed>
# Environment: ITERS (default 500), REBUILDS (default "100,300"), CORES_PER_CHAIN
# (default 9), ARM_R (an R file defining arm_options() and/or arm_hook(emc), as
# in benchmarks/subject-hmc/continue_arm.R), SAVE_STATES (default 1: save the
# state add_proposals() received at each rebuild, for replay and bootstrap).
args <- commandArgs(TRUE)
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
ckpt <- args[1]; out_dir <- args[2]; seed <- as.integer(args[3])
iters <- as.integer(Sys.getenv("ITERS", "500"))
offsets <- as.integer(strsplit(Sys.getenv("REBUILDS", "100,300"), ",")[[1]])
cores_per_chain <- as.integer(Sys.getenv("CORES_PER_CHAIN", "9"))
ns <- asNamespace("EMC2")
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "trace.R"))
arm_options <- function() list()
arm_hook <- function(emc) emc
if (nzchar(Sys.getenv("ARM_R"))) source(Sys.getenv("ARM_R"), local = TRUE)
opts <- arm_options()
if (length(opts)) do.call(options, opts)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(ckpt)), quiet = TRUE)
class(emc) <- "emc"
n_subjects <- emc[[1]]$n_subjects
p <- emc[[1]]$n_pars
base <- local({ set.seed(seed, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
n_tr <- emc[[1]]$sample_transitions
if (is.null(n_tr)) n_tr <- sum(emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)] == "sample")
for (j in seq_along(emc)) {
  base <- parallel::nextRNGStream(base)
  emc[[j]]$rng <- ns$.emc_chain_rng(base, n_subjects)
  emc[[j]]$sample_transitions <- n_tr
  emc[[j]]$sample_proposals_next <- n_tr
  emc[[j]]$ws1_trace <- NULL
}
# Production options travel with the checkpoint, so the schedule and any arm
# option go into every chain's record.
prod_opts <- opts[intersect(names(opts), names(ns$.EMC_PRODUCTION_OPTION_DEFAULTS))]
prod_opts$emc2.sample_rebuilds <- c(n_tr + offsets, n_tr + iters, 100L * (n_tr + iters))
do.call(options, prod_opts)
for (j in seq_along(emc)) {
  po <- emc[[j]]$production_options
  if (is.null(po)) po <- ns$.EMC_PRODUCTION_OPTION_DEFAULTS
  po[names(prod_opts)] <- prod_opts
  emc[[j]]$production_options <- po
}
cat("Rebuilds at transitions", n_tr, n_tr + offsets, "; run ends at", n_tr + iters, "\n")
emc <- arm_hook(emc)
ws1_install_trace(out_dir, save_states = Sys.getenv("SAVE_STATES", "1") == "1")
n0 <- emc[[1]]$samples$idx
t0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample",
  stop_criteria = list(sample = list(iter = iters)), step_size = 100L, max_tries = 1L,
  cores_for_chains = length(emc), cores_per_chain = cores_per_chain, trim = FALSE,
  verbose = FALSE, verboseProgress = FALSE)))
el <- proc.time() - t0
emc <- ns$restore_duplicates(emc)
n1 <- emc[[1]]$samples$idx
stopifnot(n1 - n0 == iters)
# The final state is the next rebuild's input.
saveRDS(list(emc = emc, stage = "sample", n_cores = cores_per_chain, n_blocks = NULL,
             block_iter = 100L, seed = NULL),
        file.path(out_dir, sprintf("pre_%02d.rds", get("ws1_k", globalenv()) + 1L)),
        compress = FALSE)
draws <- lapply(emc, function(ch) {
  k <- n0:n1
  list(mu = ch$samples$theta_mu[, k, drop = FALSE],
       lsd = 0.5 * log(apply(ch$samples$theta_var[, , k, drop = FALSE], 3, diag)),
       alpha = ch$samples$alpha[, , k, drop = FALSE])
})
saveRDS(list(ckpt = ckpt, seed = seed, iters = iters, offsets = offsets, n_tr = n_tr,
             n0 = n0, n1 = n1, elapsed = el[["elapsed"]], draws = draws,
             history = lapply(emc, function(ch) ch$alloc$history),
             par_names = emc[[1]]$par_names),
        file.path(out_dir, "run.rds"))
cat("done", el[["elapsed"]], "s\n")
