pkgload::load_all(".", quiet = TRUE, compile = FALSE)

start_file <- Sys.getenv("T30_STAGE_START",
  "benchmarks/subject-hmc/results_funnel/T30_seed1_emc.rds")
out_file <- Sys.getenv("T30_STAGE_OUT",
  "/tmp/t30_stage_transition_ess.rds")
n_draws <- as.integer(Sys.getenv("T30_STAGE_DRAWS", "200"))
n_cores <- as.integer(Sys.getenv("T30_STAGE_CORES", "12"))

start <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(start_file)), quiet = TRUE)
stopifnot(length(start) == 3L)
options(emc2.collapsed_move = TRUE, emc2.subject_hmc = FALSE,
        emc2.sampler_profile = FALSE)
for (j in seq_along(start)) {
  start[[j]]$collapsed <- NULL
  start[[j]]$sample_proposals_next <- start[[j]]$sample_transitions
}

adapt_started <- proc.time()[["elapsed"]]
adapted <- EMC2:::add_proposals(start, stage = "adapt", n_cores = n_cores,
  n_blocks = 1L, block_iter = 200L)
adapt_seconds <- proc.time()[["elapsed"]] - adapt_started
stopifnot(all(vapply(adapted, function(ch)
  is.list(ch$collapsed$plan$calibration), logical(1))))

sample_setup <- function(screen) {
  x <- unserialize(serialize(adapted, NULL))
  for (j in seq_along(x)) {
    x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
    if (!screen) x[[j]]$collapsed$plan$calibration <- NULL
  }
  tic <- proc.time()[["elapsed"]]
  x <- EMC2:::add_proposals(x, stage = "sample", n_cores = n_cores,
    n_blocks = 1L, block_iter = n_draws)
  elapsed <- proc.time()[["elapsed"]] - tic
  for (j in seq_along(x)) {
    x[[j]]$production_options <- EMC2:::.EMC_PRODUCTION_OPTION_DEFAULTS
    x[[j]]$production_options$emc2.subject_hmc <- FALSE
    x[[j]]$production_options$emc2.collapsed_move <- TRUE
  }
  list(emc = x, elapsed = elapsed,
       rounds = vapply(x, function(ch) ch$collapsed$plan$rounds, 0.0),
       pilots = vapply(x, function(ch) ch$collapsed$plan$calibration$pilots, 0.0),
       mode = vapply(x, function(ch) ch$collapsed$plan$calibration$search_mode,
                     character(1)))
}

setup <- list(full = sample_setup(FALSE), screen = sample_setup(TRUE))

summarise <- function(fit, label, sample_seconds, sample_setup_seconds) {
  p <- fit[[1L]]$n_pars
  nsub <- fit[[1L]]$n_subjects
  pn <- fit[[1L]]$par_names
  subject_ids <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  gd <- array(NA_real_, c(n_draws, length(fit), 2L * p),
    dimnames = list(NULL, NULL, c(paste0("mu:", pn), paste0("logsd:", pn))))
  ad <- array(NA_real_, c(n_draws, length(fit), p * nsub),
    dimnames = list(NULL, NULL,
      unlist(lapply(subject_ids, function(id) paste0(id, ":", pn)))))
  for (j in seq_along(fit)) {
    ch <- fit[[j]]
    idx <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"),
                n_draws)
    stopifnot(length(idx) == n_draws)
    gd[, j, seq_len(p)] <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    gd[, j, p + seq_len(p)] <- t(vapply(idx, function(i)
      0.5 * log(diag(ch$samples$theta_var[, , i])), numeric(p)))
    ad[, j, ] <- t(matrix(ch$samples$alpha[, , idx], nrow = p * nsub))
  }
  ge <- posterior::summarise_draws(posterior::as_draws_array(gd),
    posterior::ess_bulk, posterior::rhat)
  ae <- posterior::summarise_draws(posterior::as_draws_array(ad),
    posterior::ess_bulk, posterior::rhat)
  total_seconds <- adapt_seconds + sample_setup_seconds + sample_seconds
  group_logsd_rows <- match(paste0("logsd:", pn), ge$variable)
  group_ess_col <- grep("ess_bulk", names(ge), value = TRUE)
  subject_ess_col <- grep("ess_bulk", names(ae), value = TRUE)
  list(label = label, adapt_setup_seconds = adapt_seconds,
    sample_setup_seconds = sample_setup_seconds,
    sample_seconds = sample_seconds, total_seconds = total_seconds,
    rounds = setup[[label]]$rounds, pilots = setup[[label]]$pilots,
    mode = setup[[label]]$mode,
    group = ge, subject = ae,
    group_logsd_ess_per_second = median(ge[[group_ess_col]][group_logsd_rows]) /
      total_seconds,
    subject_median_ess_per_second = median(ae[[subject_ess_col]]) / total_seconds)
}

results <- list(adapt_setup_seconds = adapt_seconds)
screen_arm <- unserialize(serialize(setup$full$emc, NULL))
for (j in seq_along(screen_arm))
  screen_arm[[j]]$collapsed$plan <- setup$screen$emc[[j]]$collapsed$plan
arms <- list(full = setup$full$emc, screen = screen_arm)
for (label in c("full", "screen")) {
  tic <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arms[[label]], stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - tic
  results[[label]] <- summarise(fit, label, elapsed,
                                setup[[label]]$elapsed)
  saveRDS(results, out_file)
  cat(label, "mode", paste(unique(results[[label]]$mode), collapse = ","),
      "pilots", paste(results[[label]]$pilots, collapse = ","),
      "rounds", paste(round(results[[label]]$rounds, 2), collapse = ","),
      "total seconds", round(results[[label]]$total_seconds, 2),
      "group log-SD ESS/sec",
      round(results[[label]]$group_logsd_ess_per_second, 3),
      "subject median ESS/sec",
      round(results[[label]]$subject_median_ess_per_second, 3), "\n")
  rm(fit)
  gc()
}
saveRDS(results, out_file)
