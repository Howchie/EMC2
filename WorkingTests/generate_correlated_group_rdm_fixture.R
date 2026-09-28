# Generate an RDM fixture with correlated subject-level coefficients.
#
# This script generates data only; it does not initialize or run a fit.
# The large draw checks that make_random_effects() reproduces the requested
# coefficient correlation matrix before the 30-subject fixture is generated.
#
# Run with:
#   EMC_LIB=/path/to/EMC2 Rscript WorkingTests/generate_correlated_group_rdm_fixture.R
# Optional settings: CORR_CHECK_N, CORR_TRIALS_PER_CELL, CORR_FIXTURE_OUT,
# CORR_N_SUBJECTS (default 30) and CORR_STRENGTH (default 1): a multiplier on
# every partial correlation. 0 gives a diagonal covariance; any value in
# [0, 1/0.63) keeps the matrix positive definite by construction.
# CORR_SD_SCALE (default 1) multiplies every group standard deviation.

lib <- Sys.getenv("EMC_LIB", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))

n_subjects <- as.integer(Sys.getenv("CORR_N_SUBJECTS", "30"))
corr_strength <- as.numeric(Sys.getenv("CORR_STRENGTH", "1"))
# Uniform multiplier on every group SD (default 1).
sd_scale <- as.numeric(Sys.getenv("CORR_SD_SCALE", "1"))
# A = 0 by default (pure Wald race): a free start-point range adds an A/B/t0
# trade-off inside every subject that stalls every sampler. CORR_FREE_A=1 keeps it.
free_A <- identical(Sys.getenv("CORR_FREE_A", "0"), "1")
rdm_formula <- if (free_A) list(v ~ lM * E, B ~ E * lR, A ~ 1, t0 ~ 1) else
  list(v ~ lM * E, B ~ E * lR, t0 ~ 1)
rdm_constants <- if (free_A) c(s = log(1)) else c(s = log(1), A = log(0))
n_check <- as.integer(Sys.getenv("CORR_CHECK_N", "5000"))
n_trials <- as.integer(Sys.getenv("CORR_TRIALS_PER_CELL", "50"))
out_file <- Sys.getenv(
  "CORR_FIXTURE_OUT",
  file.path("WorkingTests", "correlated_group_rdm_fixture.rds")
)

matchfun <- function(d) d$S == d$lR
design_rdm <- design(
  factors = list(
    subjects = seq_len(n_subjects),
    E = c("speed", "neutral", "accuracy"),
    S = c("left", "right")
  ),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RDM,
  formula = rdm_formula,
  constants = rdm_constants
)

invisible(capture.output(group_means <- sampled_pars(design_rdm)))
mean_values <- c(
  v = log(2.0),
  v_lMTRUE = log(1.50),
  v_Eneutral = log(0.95),
  v_Eaccuracy = log(1.10),
  "v_lMTRUE:Eneutral" = log(0.90),
  "v_lMTRUE:Eaccuracy" = log(1.20),
  B = log(1.10),
  B_Eneutral = log(1.05),
  B_Eaccuracy = log(1.15),
  B_lRright = log(1.08),
  "B_Eneutral:lRright" = log(0.95),
  "B_Eaccuracy:lRright" = log(1.05),
  A = log(0.15),
  t0 = log(0.25)
)
mean_values <- mean_values[names(mean_values) %in% names(group_means)]
if (!setequal(names(group_means), names(mean_values))) {
  stop("Unexpected sampled parameter names: ", paste(names(group_means), collapse = ", "))
}
group_means[] <- mean_values[names(group_means)]

# Construct a positive-definite correlation matrix from mixed-sign partial
# correlations, then scale it to heterogeneous coefficient standard deviations.
p <- length(group_means)
partial <- matrix(0, p, p)
magnitudes <- c(0.18, 0.27, 0.36, 0.45, 0.54, 0.63)
for (i in seq_len(p)) {
  if (i > 1L) for (j in seq_len(i - 1L)) {
    k <- (7L * i + 11L * j) %% length(magnitudes) + 1L
    sign <- if ((2L * i + 3L * j) %% 2L == 0L) 1 else -1
    partial[i, j] <- sign * magnitudes[k] * corr_strength
  }
}
chol_corr <- diag(p)
for (i in 2:p) {
  remaining <- 1
  for (j in seq_len(i - 1L)) {
    chol_corr[i, j] <- partial[i, j] * remaining
    remaining <- remaining * sqrt(1 - partial[i, j]^2)
  }
  chol_corr[i, i] <- remaining
}
intended_cor <- chol_corr %*% t(chol_corr)
dimnames(intended_cor) <- list(names(group_means), names(group_means))

sd_values <- c(
  v = 0.12, v_lMTRUE = 0.15, v_Eneutral = 0.09, v_Eaccuracy = 0.13,
  "v_lMTRUE:Eneutral" = 0.17, "v_lMTRUE:Eaccuracy" = 0.14,
  B = 0.10, B_Eneutral = 0.08, B_Eaccuracy = 0.12, B_lRright = 0.11,
  "B_Eneutral:lRright" = 0.14, "B_Eaccuracy:lRright" = 0.09,
  A = 0.10, t0 = 0.07
)
sd_values <- sd_values[names(group_means)] * sd_scale
stopifnot(!anyNA(sd_values))
intended_cov <- diag(sd_values) %*% intended_cor %*% diag(sd_values)
dimnames(intended_cov) <- dimnames(intended_cor)
if (min(eigen(intended_cor, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  stop("Intended correlation matrix is not positive definite")
}

set.seed(20260927)
check_effects <- make_random_effects(
  design_rdm, group_means, n_subj = n_check, covariances = intended_cov
)
empirical_cor <- cor(check_effects)
cor_error <- empirical_cor - intended_cor
max_abs_error <- max(abs(cor_error))
rmse_error <- sqrt(mean(cor_error^2))
if (max_abs_error > 0.08) {
  stop("Large-draw correlation check failed: maximum absolute error = ",
       signif(max_abs_error, 3))
}

set.seed(20260928)
subject_effects <- make_random_effects(
  design_rdm, group_means, n_subj = n_subjects, covariances = intended_cov
)
simulated_data <- make_data(
  subject_effects, design_rdm, n_trials = n_trials, verbose = FALSE
)
if (!all(is.finite(simulated_data$rt))) {
  stop("Generated fixture contains non-finite response times")
}
if (length(unique(simulated_data$subjects)) != n_subjects ||
    length(unique(simulated_data$E)) != 3L || length(unique(simulated_data$S)) != 2L) {
  stop("Generated fixture does not contain all requested subjects and conditions")
}

fixture <- list(
  data = simulated_data,
  design = design_rdm,
  group_means = group_means,
  intended_covariance = intended_cov,
  intended_correlation = intended_cor,
  subject_effects = subject_effects,
  correlation_check = list(
    n_subjects = n_check,
    empirical_correlation = empirical_cor,
    error = cor_error,
    max_absolute_error = max_abs_error,
    rmse = rmse_error
  ),
  settings = list(
    model = "RDM",
    n_subjects = n_subjects,
    trials_per_design_cell = n_trials,
    corr_strength = corr_strength,
    sd_scale = sd_scale,
    formula = rdm_formula,
    constants = rdm_constants,
    seed_check = 20260927L,
    seed_fixture = 20260928L
  )
)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, out_file)

cat(sprintf("Saved data-only fixture: %s\n", normalizePath(out_file, mustWork = FALSE)))
cat(sprintf("Subjects: %d; trials: %d; conditions: 3; response thresholds: 2\n",
            n_subjects, nrow(simulated_data)))
cat(sprintf("Subject effects: %d coefficients; %d-draw correlation check; max |error| %.4f; RMSE %.4f\n",
            p, n_check, max_abs_error, rmse_error))
cat(sprintf("Intended correlations range [%.3f, %.3f]; covariance minimum eigenvalue %.5g\n",
            min(intended_cor[upper.tri(intended_cor)]),
            max(intended_cor[upper.tri(intended_cor)]),
            min(eigen(intended_cov, symmetric = TRUE, only.values = TRUE)$values)))
