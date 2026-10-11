pkgload::load_all("/data/work/EMC2_dev_oo", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit <- EMC2:::restore_custom_kernel_pointers(
  readRDS(file.path(artifact_dir, "adapted_local_ess.rds.fit.rds")), quiet = TRUE)
rounds <- vapply(fit, function(ch) ch$collapsed$plan$rounds, numeric(1))
scales <- vapply(fit, function(ch) ch$collapsed$plan$scale, numeric(1))
started <- proc.time()[["elapsed"]]
extended <- EMC2::fit(fit, stage = "sample",
  stop_criteria = list(sample = list(iter = 200L)),
  step_size = 100L, cores_for_chains = 3L, cores_per_chain = 4L,
  particle_factor = 10, max_tries = 1L, trim = FALSE,
  verbose = FALSE, verboseProgress = FALSE)
sample_seconds <- proc.time()[["elapsed"]] - started
stopifnot(identical(rounds, vapply(extended, function(ch)
  ch$collapsed$plan$rounds, numeric(1))),
  identical(scales, vapply(extended, function(ch)
    ch$collapsed$plan$scale, numeric(1))))
saveRDS(extended, file.path(artifact_dir, "adapted_local400_ess.rds.fit.rds"))
source(file.path(artifact_dir, "summarise_ess.R"))
first <- readRDS(file.path(artifact_dir, "adapted_local_ess.rds"))
summary <- .summarise_ess(extended, "adapted, chain-local, 400 draws", 400L,
  rebuild_seconds = first$rebuild_seconds,
  sample_seconds = first$sample_seconds + sample_seconds,
  cores_per_chain = 4L)
summary$continuation_seconds <- sample_seconds
saveRDS(summary, file.path(artifact_dir, "adapted_local400_ess.rds"))
cat("continued sec", round(sample_seconds, 1),
    "total sec", round(summary$total_seconds, 1),
    "collapsed acceptance", paste(round(summary$collapsed_block_acceptance, 3), collapse = ","),
    "\n")
print(summary$key_ess, row.names = FALSE)
cat("group ESS median", median(summary$group_mu_bulk),
    "52091 ESS median", median(summary$subject_52091_bulk), "\n")
