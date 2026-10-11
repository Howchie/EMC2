pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
start_file <- Sys.getenv("T30_START",
  "benchmarks/subject-hmc/results_funnel/T30_seed1_emc.rds")
out_file <- Sys.getenv("T30_OUT", "/tmp/t30_collapsed_ess.rds")
use_hmc <- identical(Sys.getenv("T30_HMC", "0"), "1")
n_draws <- as.integer(Sys.getenv("T30_DRAWS", "200"))

start <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(start_file)), quiet = TRUE)
stopifnot(length(start) == 3L)

options(emc2.collapsed_move = TRUE, emc2.subject_hmc = FALSE,
        emc2.sampler_profile = FALSE)
for (j in seq_along(start))
  start[[j]]$sample_proposals_next <- start[[j]]$sample_transitions

setup_started <- proc.time()[["elapsed"]]
start <- EMC2:::add_proposals(start, stage = "sample", n_cores = 12L,
                              n_blocks = 1L, block_iter = 200L)
setup_seconds <- proc.time()[["elapsed"]] - setup_started
allocated_rounds <- vapply(start, function(ch) ch$collapsed$plan$rounds, 0.0)

for (j in seq_along(start)) {
  start[[j]]$production_options$emc2.subject_hmc <- use_hmc
  start[[j]]$production_options$emc2.collapsed_move <- TRUE
}

make_arm <- function(enabled) {
  out <- unserialize(serialize(start, NULL))
  for (j in seq_along(out)) {
    out[[j]]$collapsed$plan$rounds <- if (enabled) allocated_rounds[j] else 0
    out[[j]]$collapsed$plan$last <- EMC2:::.emc_collapsed_record(NULL,
      out[[j]]$n_pars)
    out[[j]]$collapsed$plan$pool <- EMC2:::.emc_collapsed_record(NULL,
      out[[j]]$n_pars)
    out[[j]]$collapsed$plan$block <- EMC2:::.emc_collapsed_record(NULL,
      out[[j]]$n_pars)
    out[[j]]$production_options$emc2.collapsed_move <- enabled
  }
  out
}

summarise <- function(fit, label, elapsed) {
  nd <- n_draws
  nchain <- length(fit)
  p <- fit[[1L]]$n_pars
  nsub <- fit[[1L]]$n_subjects
  pn <- fit[[1L]]$par_names
  subject_ids <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  gd <- array(NA_real_, c(nd, nchain, 2L * p),
    dimnames = list(NULL, NULL, c(paste0("mu:", pn), paste0("logsd:", pn))))
  ad <- array(NA_real_, c(nd, nchain, p * nsub),
    dimnames = list(NULL, NULL,
      unlist(lapply(subject_ids, function(id) paste0(id, ":", pn)))))
  for (j in seq_len(nchain)) {
    ch <- fit[[j]]
    idx <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), nd)
    stopifnot(length(idx) == nd)
    gd[, j, seq_len(p)] <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    lsd <- vapply(idx, function(i) 0.5 * log(diag(ch$samples$theta_var[, , i])),
                  numeric(p))
    gd[, j, p + seq_len(p)] <- t(lsd)
    ad[, j, ] <- t(matrix(ch$samples$alpha[, , idx], nrow = p * nsub))
  }
  ge <- posterior::summarise_draws(posterior::as_draws_array(gd),
                                    posterior::ess_bulk,
                                    posterior::ess_tail,
                                    posterior::rhat)
  ae <- posterior::summarise_draws(posterior::as_draws_array(ad),
                                    posterior::ess_bulk,
                                    posterior::ess_tail,
                                    posterior::rhat)
  names(ge)[2:4] <- c("ess_bulk", "ess_tail", "rhat")
  names(ae)[2:4] <- c("ess_bulk", "ess_tail", "rhat")
  list(label = label, sample_seconds = elapsed,
       setup_seconds = if (label == "collapsed") setup_seconds else 0,
       total_seconds = elapsed + if (label == "collapsed") setup_seconds else 0,
       group = ge, subject = ae,
       collapsed_acceptance = vapply(fit, function(ch) {
         b <- ch$collapsed$plan$block
         if (is.null(b$n) || b$n == 0L) NA_real_ else b$acc / b$n
       }, 0.0),
       collapsed_attempts = vapply(fit, function(ch)
         ch$collapsed$plan$block$n, 0.0))
}

results <- list(allocated_rounds = allocated_rounds,
                setup_seconds = setup_seconds)
for (label in c("collapsed", "control")) {
  enabled <- identical(label, "collapsed")
  options(emc2.collapsed_move = enabled, emc2.subject_hmc = use_hmc)
  tic <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(make_arm(enabled), stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - tic
  results[[label]] <- summarise(fit, label, elapsed)
  saveRDS(results, out_file)
  cat(label, "sample seconds", round(elapsed, 2),
      "total seconds", round(results[[label]]$total_seconds, 2),
      "collapsed acceptance",
      paste(round(results[[label]]$collapsed_acceptance, 3), collapse = ","),
      "\n")
  rm(fit)
  gc()
}
saveRDS(results, out_file)
