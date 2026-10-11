args <- commandArgs(trailingOnly = TRUE)
pkg <- args[[1L]]
label <- args[[2L]]
out_file <- args[[3L]]
continuation <- 200L
cores_per_chain <- 4L
pkgload::load_all(pkg, quiet = TRUE, compile = FALSE)

x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
options(emc2.subject_hmc = FALSE)
for (j in seq_along(x)) {
  x[[j]]$production_options$emc2.subject_hmc <- FALSE
  # Force the same fresh sample-stage proposal rebuild at the starting state
  # in both arms, so the adapted arm calibrates and the baseline arm uses its
  # pre-change collapsed transport/allocation.
  x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
}

started <- proc.time()[["elapsed"]]
x <- EMC2:::add_proposals(x, "sample", n_cores = length(x) * cores_per_chain, n_blocks = 1L,
                          block_iter = continuation)
rebuild_seconds <- proc.time()[["elapsed"]] - started
rounds <- vapply(x, function(ch) ch$collapsed$plan$rounds, numeric(1))
scales <- vapply(x, function(ch) if (is.null(ch$collapsed$plan$scale)) 1 else
  ch$collapsed$plan$scale, numeric(1))
calibrations <- lapply(x, function(ch) ch$collapsed$plan$calibration)
message(label, " initial collapsed rounds: ", paste(rounds, collapse = ","),
        "; scales: ", paste(signif(scales, 4), collapse = ","),
        "; cores/chain: ", cores_per_chain,
        "; calibration sec: ", paste(vapply(calibrations, function(z)
          if (is.null(z)) NA_real_ else z$elapsed, numeric(1)), collapse = ","),
        "; rebuild sec: ", round(rebuild_seconds, 1))

started <- proc.time()[["elapsed"]]
fit <- EMC2::fit(x, stage = "sample",
  stop_criteria = list(sample = list(iter = continuation)),
  step_size = 100L, cores_for_chains = length(x), cores_per_chain = cores_per_chain,
  particle_factor = 10, max_tries = 1L, trim = FALSE,
  verbose = FALSE, verboseProgress = FALSE)
sample_seconds <- proc.time()[["elapsed"]] - started
saveRDS(fit, paste0(out_file, ".fit.rds"))

source(file.path(dirname(out_file), "summarise_ess.R"))
summary <- .summarise_ess(fit, label, continuation, rebuild_seconds,
                          sample_seconds, cores_per_chain)
saveRDS(summary, out_file)
cat("\n", label, " elapsed", round(summary$total_seconds, 1),
    " sec; rounds", paste(summary$collapsed_rounds, collapse = ","), "\n")
print(summary$key_ess, row.names = FALSE)
cat("group mu ESS bulk range:", range(summary$group_mu_bulk),
    " median:", median(summary$group_mu_bulk), "\n")
cat("52091 alpha ESS bulk range:", range(summary$subject_52091_bulk),
    " median:", median(summary$subject_52091_bulk), "\n")
