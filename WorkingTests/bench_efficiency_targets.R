# Short, direct native likelihood baseline for the architecture-efficiency plan.
#
# Examples:
#   Rscript WorkingTests/bench_efficiency_targets.R
#   MODEL=PCOUNTER TRIALS=600 CELLS=4 PARTICLES=25 REPS=9 \
#     Rscript WorkingTests/bench_efficiency_targets.R
#   MODEL=PCOUNTER TRIALS=600 CELLS=255,256,257,512 PARTICLES=5 \
#     BRANCH=varying PROFILE=1 Rscript WorkingTests/bench_efficiency_targets.R
#
# This intentionally does not call fit(), run_emc(), or a worker pool.  The
# optional profiling pass is a separate explanatory run and must not be used as
# the speed number.

lib <- Sys.getenv("EMC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "WorkingTests/bench_efficiency_targets.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
helper <- file.path(root, "tests", "testthat", "helper-audit-harness.R")
if (!file.exists(helper)) stop("cannot locate audit helper at ", helper)
source(helper, local = .GlobalEnv)

read_ints <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, default)
  vals <- suppressWarnings(as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
  vals <- vals[is.finite(vals) & vals >= minimum]
  if (!length(vals)) as.integer(default) else unique(vals)
}
read_cells <- function() {
  raw <- trimws(Sys.getenv("CELLS", "1"))
  if (!nzchar(raw) || tolower(raw) %in% c("null", "none")) return(list(NULL))
  vals <- suppressWarnings(as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
  vals <- vals[is.finite(vals) & vals >= 1L]
  if (!length(vals)) list(1L) else lapply(unique(vals), as.integer)
}
read_bool <- function(name, default = FALSE) {
  raw <- tolower(trimws(Sys.getenv(name, if (default) "1" else "0")))
  raw %in% c("1", "true", "t", "yes", "y", "on")
}

model_tokens <- trimws(strsplit(Sys.getenv("MODEL", "PCOUNTER"), ",", fixed = TRUE)[[1L]])
model_names <- c(BAWL = "BAwL", PCOUNTER = "PCOUNTER", RDM = "RDM", DDM = "DDM",
                 LBA = "LBA")
models <- unname(model_names[toupper(model_tokens)])
models <- models[!is.na(models) & nzchar(models)]
if (!length(models)) models <- "PCOUNTER"
trials <- read_ints("TRIALS", "200")
particles <- read_ints("PARTICLES", "5")
cells <- read_cells()
branches <- trimws(strsplit(Sys.getenv("BRANCH", "varying"), ",", fixed = TRUE)[[1L]])
branches <- branches[nzchar(branches)]
if (!length(branches)) branches <- "varying"
reps <- max(7L, read_ints("REPS", "7")[[1L]])
profile <- read_bool("PROFILE", FALSE)

fmt <- function(x) {
  if (length(x) != 1L || is.na(x)) "NA" else format(x, digits = 17, scientific = TRUE, trim = TRUE)
}

run_case <- function(model, nt, nc, np, branch) {
  cell_label <- if (is.null(nc)) "null" else as.character(nc)
  label <- sprintf("model=%s trials=%d cells=%s particles=%d branch=%s",
                   model, nt, cell_label, np, branch)
  # A one-cell workload is the intercept-only design, not a rank-deficient
  # covariate with one observed level.
  fixture_cells <- if (!is.null(nc) && identical(as.integer(nc), 1L)) NULL else nc
  cat("\nCASE ", label, "\n", sep = "")
  fx <- suppressMessages(audit_fixture(
    model = model, n_trials = nt, n_particles = np,
    cells = fixture_cells,
    branch = if (identical(model, "PCOUNTER")) branch else NULL,
    seed = 91026L, rt_resolution = NULL))
  cat(sprintf("  fixture dadm_rows=%d proposal_rows=%d proposal_cols=%d\n",
              nrow(fx$dadm), nrow(fx$prop), ncol(fx$prop)))

  # Speed numbers are always from the unprofiled pass.
  fast <- audit_time_direct(fx, reps = reps, profile = FALSE)
  cat(sprintf("  speed median_s=%s mad_s=%s cpu_median_s=%s cpu_total_s=%s",
              fmt(fast$median_seconds), fmt(fast$mad_seconds),
              fmt(fast$cpu_seconds), fmt(fast$cpu_total_seconds)))
  cat(sprintf(" checksum=%s native_calls=%d all_finite=%s invalid_floor=%s checksum_reproducible=%s\n",
              fmt(fast$likelihood_checksum), fast$native_calls,
              fast$all_finite, fast$invalid_floor, fast$checksum_reproducible))
  cat(sprintf("  speed kernel_rows=%s kernel_cells=%s kernel_seconds=%s\n",
              fmt(fast$kernel_rows), fmt(fast$kernel_cells), fmt(fast$kernel_seconds)))

  # One direct native call on one particle is a smoke check, not another
  # sampler step and not part of the speed median.
  one <- fx
  one$prop <- fx$prop[1L, , drop = FALSE]
  one_tm <- system.time(one_ll <- audit_ll_direct(one))
  one_n_likelihood <- length(attr(one$dadm, "expand"))
  if (!one_n_likelihood) one_n_likelihood <- nrow(one$dadm)
  one_floor <- any(is.finite(one_ll) &
                   abs(as.numeric(one_ll) - one_n_likelihood * log(1e-10)) <= 1e-12)
  cat(sprintf("  particle_step cpu_s=%s checksum=%s native_calls=1 invalid_floor=%s\n",
              fmt(one_tm[["user.self"]] + one_tm[["sys.self"]]),
              fmt(sum(as.numeric(one_ll))), one_floor))

  if (profile) {
    explain <- audit_time_direct(fx, reps = reps, profile = TRUE)
    cat(sprintf("  explanation(profile=on) median_s=%s mad_s=%s cpu_median_s=%s\n",
                fmt(explain$median_seconds), fmt(explain$mad_seconds),
                fmt(explain$cpu_seconds)))
    cat(sprintf("  explanation kernel_rows=%s kernel_cells=%s kernel_seconds=%s kernel_calls=%s\n",
                fmt(explain$kernel_rows), fmt(explain$kernel_cells),
                fmt(explain$kernel_seconds), fmt(explain$kernel_calls)))
    if (is.data.frame(explain$pcounter) && nrow(explain$pcounter)) {
      pc <- explain$pcounter[1L, , drop = FALSE]
      cat(sprintf(paste0("  pcounter d_calls=%d p_calls=%d logS_calls=%d",
                         " rows=%d prep_s=%s rt_s=%s",
                         " sv_zero=%d gamma_zero=%d omega_zero=%d",
                         " k_min=%d k_max=%d k_mean=%s stirling=%d rising=%d\n"),
                  pc$dpcounter_calls, pc$ppcounter_calls,
                  pc$pcounter_logS_calls,
                  pc$dpcounter_rows + pc$ppcounter_rows + pc$pcounter_logS_rows,
                  fmt(pc$preparation_seconds), fmt(pc$rt_sum_seconds),
                  pc$sv_zero, pc$gamma_zero, pc$omega_zero,
                  pc$k_min, pc$k_max,
                  fmt(if (pc$k_observations > 0) pc$k_sum / pc$k_observations else NA),
                  pc$stirling_terms, pc$rising_terms))
    }
  }
  invisible(fast)
}

for (model in models) for (nt in trials) for (nc in cells) for (np in particles)
  for (branch in branches) run_case(model, nt, nc, np, branch)
