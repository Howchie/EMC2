#!/usr/bin/env Rscript

# Benchmark the isolated adaptive-GK candidate against fixed Gauss-Hermite
# rules on a vector of shared-latent conditional-normal likelihoods.  This is
# a proxy for the BAwLcorr integration shape; it deliberately does not alter
# or call the production BAwLcorr likelihood.

suppressPackageStartupMessages({
  library(Rcpp)
  library(statmod)
})

source_file <- file.path("tools", "adaptive-gk-probe.cpp")
cache <- file.path(tempdir(), paste0("emc2-adaptive-gk-bench-", Sys.getpid()))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
Rcpp::sourceCpp(source_file, cacheDir = cache, rebuild = TRUE,
                showOutput = FALSE, verbose = FALSE)

n <- as.integer(Sys.getenv("EMC2_ADAPTIVE_GK_BENCH_N", "5000"))
reps <- as.integer(Sys.getenv("EMC2_ADAPTIVE_GK_BENCH_REPS", "5"))
rel_tol <- as.numeric(Sys.getenv("EMC2_ADAPTIVE_GK_BENCH_REL_TOL", "1e-4"))
if (!is.finite(n) || n < 1L || !is.finite(reps) || reps < 1L)
  stop("EMC2_ADAPTIVE_GK_BENCH_N and _REPS must be positive integers")
if (!is.finite(rel_tol) || rel_tol <= 0)
  stop("EMC2_ADAPTIVE_GK_BENCH_REL_TOL must be positive")

set.seed(20260716)
rho <- runif(n, .2, .95)
means <- rnorm(n, sd = .35)
loadings <- .8 * sqrt(rho)
residual_sd <- .8 * sqrt(1 - rho)
latent_truth <- rnorm(n, sd = .9)
observations <- means + loadings * latent_truth +
  residual_sd * rnorm(n)
exact <- dnorm(observations, means,
               sqrt(loadings^2 + residual_sd^2), log = TRUE)

# The candidate integrates over x = Phi(z).  These breaks cover the tails
# without imposing a production policy for where a future BAwLcorr call would
# seed its partition.
breaks <- pnorm(seq(-7, 7, by = 1))

candidate <- function() {
  adaptive_gk_latent_normal(
    means, loadings, residual_sd, observations,
    abs_tol = 1e-11, rel_tol = rel_tol, limit = 128,
    initial_breaks = breaks)
}

gh <- function(n_nodes) {
  rule <- statmod::gauss.quad(n_nodes, kind = "hermite")
  gh_latent_normal(means, loadings, residual_sd, observations,
                   rule$nodes, rule$weights)
}

# Warm compilation and the candidate's allocation paths.
candidate_result <- candidate()
if (!candidate_result$converged)
  warning("adaptive candidate hit its interval limit")

candidate_times <- vapply(seq_len(reps), function(i)
  system.time(candidate())[["elapsed"]], numeric(1))

nodes <- c(40L, 80L, 200L)
gh_times <- lapply(nodes, function(k) vapply(seq_len(reps), function(i)
  system.time(gh(k))[["elapsed"]], numeric(1)))

error_summary <- function(x) {
  e <- abs(x - exact)
  c(median = median(e), p95 = unname(quantile(e, .95)),
    p99 = unname(quantile(e, .99)), max = max(e),
    above_1e4 = sum(e > 1e-4))
}

error_rows <- c(list(error_summary(candidate_result$log_integral)),
                lapply(nodes, function(k) error_summary(gh(k))))
error_rows <- as.data.frame(do.call(rbind, error_rows))

rows <- data.frame(
  method = c("adaptive_gk", paste0("gauss_hermite_", nodes)),
  median_seconds = c(median(candidate_times),
                     vapply(gh_times, median, numeric(1))),
  eval_points = c(candidate_result$eval_points, nodes),
  error_rows,
  check.names = FALSE
)
print(rows, row.names = FALSE, digits = 6)
cat(sprintf("adaptive_rel_tol = %.3g\n", rel_tol))
cat(sprintf("adaptive_converged = %s\n", candidate_result$converged))
cat(sprintf("adaptive_intervals = %d\n", candidate_result$intervals))
for (i in seq_along(nodes)) {
  cat(sprintf("speedup_vs_gh_%d = %.3fx\n", nodes[i],
              rows$median_seconds[i + 1L] / rows$median_seconds[1L]))
}
