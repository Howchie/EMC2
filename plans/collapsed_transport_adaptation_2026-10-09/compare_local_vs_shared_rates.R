pkgload::load_all("/data/work/EMC2_dev_oo", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
options(emc2.subject_hmc = FALSE)
for (j in seq_along(x)) {
  x[[j]]$production_options$emc2.subject_hmc <- FALSE
  x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
}

started <- proc.time()[["elapsed"]]
x <- EMC2:::add_proposals(x, "sample", n_cores = 12L, n_blocks = 1L,
                          block_iter = 200L)
rebuild_seconds <- proc.time()[["elapsed"]] - started
local_rates <- vapply(x, function(ch) ch$collapsed$plan$rounds, numeric(1))
scales <- vapply(x, function(ch) ch$collapsed$plan$scale, numeric(1))
stopifnot(length(unique(local_rates)) > 1L)
shared <- unserialize(serialize(x, NULL, version = 3L))
shared_rate <- max(local_rates)
for (j in seq_along(shared)) shared[[j]]$collapsed$plan$rounds <- shared_rate

sample_arm <- function(start, label) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(start, stage = "sample",
    stop_criteria = list(sample = list(iter = 200L)),
    step_size = 100L, cores_for_chains = 3L, cores_per_chain = 4L,
    particle_factor = 10, max_tries = 1L, trim = FALSE,
    verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  list(summary = .summarise_ess(fit, label, 200L,
    rebuild_seconds = rebuild_seconds, sample_seconds = elapsed,
    cores_per_chain = 4L), elapsed = elapsed)
}
source(file.path(artifact_dir, "summarise_ess.R"))
# The arms share the same calibrated state, scale, and RNG seeds. Only the
# collapsed round policy differs; each arm gets its own 200-draw continuation.
local <- sample_arm(x, "chain-local rates")
global <- sample_arm(shared, "shared maximum rate")
result <- list(local = local$summary, shared = global$summary,
  local_rates = local_rates, shared_rate = shared_rate, scales = scales,
  rebuild_seconds = rebuild_seconds,
  local_sample_seconds = local$elapsed, shared_sample_seconds = global$elapsed,
  cores_for_chains = 3L, cores_per_chain = 4L)
saveRDS(result, file.path(artifact_dir, "local_vs_shared_rates.rds"))
cat("local rates", paste(local_rates, collapse = ","),
    "shared rate", shared_rate, "scales", paste(signif(scales, 5), collapse = ","),
    "rebuild sec", round(rebuild_seconds, 1),
    "local sample sec", round(local$elapsed, 1),
    "shared sample sec", round(global$elapsed, 1), "\n")
for (arm in c("local", "shared")) {
  z <- result[[arm]]
  cat("\n", arm, "collapsed acceptance",
      paste(round(z$collapsed_block_acceptance, 3), collapse = ","),
      "group ESS median", median(z$group_mu_bulk),
      "52091 ESS median", median(z$subject_52091_bulk), "\n")
  print(z$key_ess, row.names = FALSE)
}
