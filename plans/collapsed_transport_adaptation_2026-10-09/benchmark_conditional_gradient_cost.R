pkgload::load_all(".", quiet = TRUE)

archive_root <- Sys.getenv(
  "EMC2_FULLC_ROOT",
  ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final")
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
cases <- data.frame(
  case = c("T30", "corr15", "Forstmann", "LBA"),
  file = c("T30_fullc_1_emc.rds", "corr15_fullc_1_emc.rds",
           "forstmann_fullc_2_emc.rds", NA_character_),
  chain = c(1L, 1L, 2L, 1L),
  stringsAsFactors = FALSE)
n_eval <- as.integer(Sys.getenv("EMC2_GRAD_EVALS", "250"))
n_reps <- as.integer(Sys.getenv("EMC2_GRAD_REPS", "6"))

load_case <- function(i) {
  if (cases$case[i] == "LBA") {
    path <- Sys.getenv("EMC2_LBA_RDS", "/tmp/lba_sampling_audit/snapshot.rds")
    fit <- EMC2:::restore_custom_kernel_pointers(
      EMC2:::restore_duplicates(readRDS(path)), quiet = TRUE)
  } else {
    path <- file.path(archive_root, cases$file[i])
    fit <- EMC2:::restore_custom_kernel_pointers(
      EMC2:::restore_duplicates(readRDS(path)), quiet = TRUE)
  }
  fit[[cases$chain[i]]]
}

measure_case <- function(i) {
  chain <- load_case(i)
  ctx_fixed <- EMC2:::.emc_collapsed_context(chain)
  ctx_free <- EMC2:::.emc_collapsed_context(
    chain, fixed_idx = integer(),
    fixed_value = matrix(numeric(), 0L, chain$n_subjects))
  it <- chain$samples$idx
  cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, it],
    as.matrix(chain$samples$theta_var[, , it]), ctx_fixed$fits,
    chain$prior, ctx_fixed$spec)
  stopifnot(!is.null(cur), is.finite(cur$lp))
  phi <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx_fixed$chart)
  stopifnot(is.finite(EMC2:::collapsed_surrogate_cpp(phi, ctx_fixed$cpp)$lp),
    is.finite(EMC2:::collapsed_surrogate_cpp(phi, ctx_free$cpp)$lp))

  # Warm both paths, then alternate their order to reduce drift from CPU load.
  for (j in seq_len(10L)) {
    EMC2:::collapsed_surrogate_cpp(phi, ctx_fixed$cpp)
    EMC2:::collapsed_surrogate_cpp(phi, ctx_free$cpp)
  }
  timing <- matrix(NA_real_, n_reps, 2L,
    dimnames = list(NULL, c("conditioned", "unconditioned")))
  for (r in seq_len(n_reps)) {
    order <- if (r %% 2L) c("conditioned", "unconditioned") else
      c("unconditioned", "conditioned")
    for (arm in order) {
      model <- if (arm == "conditioned") ctx_fixed$cpp else ctx_free$cpp
      started <- proc.time()[["elapsed"]]
      for (j in seq_len(n_eval))
        EMC2:::collapsed_surrogate_cpp(phi, model)
      timing[r, arm] <- (proc.time()[["elapsed"]] - started) / n_eval
    }
  }
  data.frame(case = cases$case[i], subjects = chain$n_subjects,
    group_parameters = chain$n_pars, fixed_subject_coefficients =
      length(ctx_fixed$fixed_idx), evaluations_per_rep = n_eval,
    reps = n_reps, conditioned_median_ms = 1000 * median(timing[, "conditioned"]),
    unconditioned_median_ms = 1000 * median(timing[, "unconditioned"]),
    ratio = median(timing[, "conditioned"]) / median(timing[, "unconditioned"]),
    conditioned_iqr_ms = 1000 * IQR(timing[, "conditioned"]),
    unconditioned_iqr_ms = 1000 * IQR(timing[, "unconditioned"]))
}

result <- do.call(rbind, lapply(seq_len(nrow(cases)), measure_case))
write.csv(result, file.path(artifact_dir,
  "conditional_gradient_cost_current.csv"), row.names = FALSE)
print(result, row.names = FALSE)
