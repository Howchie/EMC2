pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(artifact_dir, "summarise_ess.R"))
source(file.path(artifact_dir, "previous_budget_solver.R"))
setup_only_solver <- previous_budget_solver()
x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")), quiet = TRUE)
options(emc2.subject_hmc = FALSE)
for (j in seq_along(x)) {
  x[[j]]$production_options$emc2.subject_hmc <- FALSE
  x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
}
.captured <- list()
trace(".emc_alloc_solve_units", where = asNamespace("EMC2"), print = FALSE,
  exit = quote({
    .GlobalEnv$.captured[[length(.GlobalEnv$.captured) + 1L]] <-
      list(units = u, solution = returnValue())
  }))
started <- proc.time()[["elapsed"]]
x <- EMC2:::add_proposals(x, "sample", n_cores = 12L, n_blocks = 1L, block_iter = 200L)
rebuild_seconds <- proc.time()[["elapsed"]] - started
untrace(".emc_alloc_solve_units", where = asNamespace("EMC2"))
stopifnot(length(.captured) == length(x))
starts <- list(corrected = x, geometric = unserialize(serialize(x, NULL)),
               no_collapsed = unserialize(serialize(x, NULL)))
solutions <- lapply(.captured, function(z) {
  u <- z$units
  old <- u
  old$joint$a <- old$joint$g
  old$joint$fixed_cost <- 0
  none <- u
  none$joint <- NULL
  list(corrected = setup_only_solver(u),
       geometric = setup_only_solver(old),
       no_collapsed = EMC2:::.emc_alloc_solve_units(none))
})
for (arm in names(starts)) {
  for (j in seq_along(x)) {
    sol <- solutions[[j]][[arm]]
    starts[[arm]][[j]]$collapsed$plan$rounds <- sol$joint
    if (EMC2:::.emc_shear_active(starts[[arm]][[j]]))
      starts[[arm]][[j]]$shear$plan$rounds <- sol$x[seq_len(2L * x[[j]]$n_pars)]
  }
  cat(arm, "rates", paste(vapply(starts[[arm]], function(ch)
    ch$collapsed$plan$rounds, 0), collapse = ","), "\n")
}
for (j in seq_along(x)) {
  u <- .captured[[j]]$units
  cat("chain", j, "scale", x[[j]]$collapsed$plan$scale,
      "marginal cost", u$joint$cost, "setup cost", u$joint$fixed_cost,
      "saturating coordinates", sum(u$joint$a < u$joint$g - 1e-8),
      "of", length(u$joint$g), "\n")
}
saveRDS(list(starts = starts, captured = .captured, solutions = solutions),
        file.path(fit_dir, "setup_saturation_start.rds"))
# The two repetitions start at identical states but use distinct chain streams.
# Reversed arm order in the second repetition reduces systematic timing bias.
results <- list()
for (replicate in 1:2) {
  order <- if (replicate == 1L) names(starts) else rev(names(starts))
  for (arm in order) {
    start <- starts[[arm]]
    if (replicate > 1L) for (j in seq_along(start)) {
      advance <- function(seed) {
        for (k in seq_len(length(start))) seed <- parallel::nextRNGStream(seed)
        seed
      }
      start[[j]]$rng$chain <- advance(start[[j]]$rng$chain)
      start[[j]]$rng$gibbs <- advance(start[[j]]$rng$gibbs)
      start[[j]]$rng$subjects <- lapply(start[[j]]$rng$subjects, advance)
    }
    started <- proc.time()[["elapsed"]]
    fit <- EMC2::fit(start, stage = "sample",
      stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
      cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
      max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started
    key <- paste(arm, replicate, sep = "_")
    saveRDS(fit, file.path(fit_dir, paste0("setup_saturation_", key, ".fit.rds")))
    z <- .summarise_ess(fit, key, 200L, rebuild_seconds, elapsed, 4L)
    results[[key]] <- z
    saveRDS(results, file.path(artifact_dir, "setup_saturation_ess.rds"))
    cat(key, "seconds", round(elapsed, 1), "group ESS", median(z$group_mu_bulk),
        "subject ESS", median(z$subject_52091_bulk), "\n")
    flush.console()
  }
}
