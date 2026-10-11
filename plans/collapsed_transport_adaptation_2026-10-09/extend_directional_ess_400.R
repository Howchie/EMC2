pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
options(emc2.subject_hmc = FALSE)

continue_one <- function(label, input, output) {
  x <- EMC2:::restore_custom_kernel_pointers(readRDS(input), quiet = TRUE)
  before <- vapply(x, function(ch)
    sum(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), integer(1))
  stopifnot(length(unique(before)) == 1L, before >= 200L)
  started <- proc.time()[["elapsed"]]
  x <- EMC2::fit(x, stage = "sample",
    stop_criteria = list(iter = 200L), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE,
    verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  after <- vapply(x, function(ch)
    sum(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), integer(1))
  stopifnot(all(after - before == 200L))
  saveRDS(x, output)
  summary <- .summarise_ess(x, label, 400L,
    rebuild_seconds = if (label == "fixed_scale") 45.601 else 150.075,
    sample_seconds = if (label == "fixed_scale") 399.101 + elapsed
                    else 410.555 + elapsed,
    cores_per_chain = 4L)
  par_names <- dimnames(x[[1L]]$samples$theta_mu)[[1L]]
  rm(x)
  gc()
  list(summary = summary, extension_seconds = elapsed,
       par_names = par_names)
}

fixed <- continue_one("fixed_scale",
  file.path(fit_dir, "safe_fit_ess_safe_fit_split_t0.fit.rds"),
  file.path(fit_dir, "safe_fit_ess_fixed_400.fit.rds"))
saveRDS(fixed, file.path(artifact_dir, "fixed_ess_400_summary.rds"))
cat("fixed extension seconds", round(fixed$extension_seconds, 1), "\n")
flush.console()

directional <- continue_one("directional_scale",
  file.path(fit_dir, "safe_fit_ess_directional_200.fit.rds"),
  file.path(fit_dir, "safe_fit_ess_directional_400.fit.rds"))
saveRDS(directional,
  file.path(artifact_dir, "directional_ess_400_summary.rds"))
cat("directional extension seconds", round(directional$extension_seconds, 1), "\n")
flush.console()

fixed_summary <- fixed$summary
directional_summary <- directional$summary
vg <- match("v.g", fixed$par_names)
table <- do.call(rbind, lapply(list(fixed = fixed_summary,
                                    directional = directional_summary),
  function(z) data.frame(
    arm = z$label,
    draws_per_chain = z$draws_per_chain,
    first_200_sample_seconds = if (z$label == "fixed_scale") 399.101
                               else 410.555,
    extension_seconds = if (z$label == "fixed_scale")
      fixed$extension_seconds else directional$extension_seconds,
    rebuild_seconds = z$rebuild_seconds,
    sample_seconds = z$sample_seconds,
    total_seconds = z$total_seconds,
    median_group_ess = median(z$group_mu_bulk),
    median_group_ess_per_total_second =
      median(z$group_mu_bulk) / z$total_seconds,
    median_group_ess_per_sample_second =
      median(z$group_mu_bulk) / z$sample_seconds,
    median_52091_ess = median(z$subject_52091_bulk),
    median_52091_ess_per_total_second =
      median(z$subject_52091_bulk) / z$total_seconds,
    v_g_group_ess = z$group_mu_bulk[vg],
    v_g_ess_per_total_second = z$group_mu_bulk[vg] / z$total_seconds,
    median_key_ess = median(z$key_ess$bulk),
    max_key_rhat = max(z$key_ess$rhat)
  )))
write.csv(table,
  file.path(artifact_dir, "directional_ess_400_summary.csv"),
  row.names = FALSE)
print(table, row.names = FALSE)
