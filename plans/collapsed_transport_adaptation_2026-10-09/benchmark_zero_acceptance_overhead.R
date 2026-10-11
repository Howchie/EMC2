# Measure the cost of leaving the original unit-scale collapsed proposal
# active after its exact acceptance has collapsed. Compare with the same
# calibrated LBA start and all other moves unchanged, with no collapsed rounds.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("plans/collapsed_transport_adaptation_2026-10-09/summarise_ess.R")

fit_dir <- "/tmp/collapsed_transport_followup"
out_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
n_draws <- as.integer(Sys.getenv("ZERO_ACCEPT_DRAWS", "100"))
prepared <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))
base <- prepared$start
options(emc2.subject_hmc = FALSE, emc2.collapsed_move = TRUE)
for (j in seq_along(base)) {
  base[[j]]$production_options$emc2.subject_hmc <- FALSE
  base[[j]]$sample_proposals_next <- Inf
  base[[j]]$collapsed$plan$rounds <- 0
}

unit_scale <- unserialize(serialize(base, NULL))
ch <- unit_scale[[2L]]
plan <- ch$collapsed$plan
ctx <- EMC2:::.emc_collapsed_context(ch)
phi <- EMC2:::.emc_collapsed_to_phi(
  ch$samples$theta_mu[, ch$samples$idx],
  as.matrix(ch$samples$theta_var[, , ch$samples$idx]), ctx$chart)
plan$scale <- 1
plan$t0_ratio <- 1
plan$t0_mean_ratio <- 1
plan$t0_hyper_ratio <- 1
plan$mean_idx <- integer()
plan$mean_ratio <- 1
plan$length <- 3
C <- EMC2:::.emc_collapsed_preconditioner(plan$C, 1, 1, plan$t0_idx,
  1, plan$t0_mean_idx)
plan$eps <- EMC2:::.emc_with_preserved_rng(
  EMC2:::.emc_collapsed_tune_eps(phi, C, plan$eps, ctx, plan$length))
plan$calibrate <- FALSE
plan$calibration_screen <- FALSE
plan$rounds <- 4
plan$block <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
ch$collapsed$plan <- plan
unit_scale[[2L]] <- ch

cgroup_path <- file.path("/sys/fs/cgroup", sub("^0::", "",
  readLines("/proc/self/cgroup")[1L]))
cpu_seconds <- function() {
  x <- tryCatch(readLines(file.path(cgroup_path, "cpu.stat")),
                error = function(e) character())
  u <- grep("^usage_usec", x, value = TRUE)
  if (!length(u)) return(NA_real_)
  as.numeric(sub("usage_usec ", "", u)) / 1e6
}

arms <- list(unit_scale_four_rounds = unit_scale,
             no_collapsed = unserialize(serialize(base, NULL)))
results <- list()
for (arm in names(arms)) {
  cpu0 <- cpu_seconds()
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arms[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  wall <- proc.time()[["elapsed"]] - started
  cpu <- cpu_seconds() - cpu0
  summary <- .summarise_ess(fit, arm, n_draws, 0, wall, 4L)
  summary$cgroup_cpu_seconds <- cpu
  summary$collapsed_attempts <- vapply(fit, function(x)
    x$collapsed$plan$block$n, integer(1))
  summary$collapsed_acceptance_probability <- vapply(fit, function(x) {
    b <- x$collapsed$plan$block
    if (b$n > 0L) b$acc / b$n else NA_real_
  }, numeric(1))
  summary$collapsed_accepted_rounds <- vapply(fit, function(x)
    x$collapsed$plan$block$moved, integer(1))
  results[[arm]] <- summary
  cat(arm, "wall", round(wall, 1), "CPU", round(cpu, 1),
    "collapsed attempts", paste(summary$collapsed_attempts, collapse = "/"),
    "mean exact acceptance probability",
    paste(round(summary$collapsed_acceptance_probability, 6), collapse = "/"),
    "accepted rounds", paste(summary$collapsed_accepted_rounds, collapse = "/"),
    "\n")
  saveRDS(results, file.path(out_dir, "zero_acceptance_overhead.rds"))
  rm(fit)
  gc()
}
