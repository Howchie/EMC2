#!/usr/bin/env Rscript

compare_wald_logspace <- function() {
  lib_loc <- Sys.getenv("EMC2_OPTIMIZED_LIB", "/data/tmp/r_libs/optimized_version")
  stopifnot(requireNamespace("EMC2", lib.loc = lib_loc, quietly = TRUE))
  pigt_hybrid <- getExportedValue("EMC2", "pigt")
  pigt_log <- getExportedValue("EMC2", "pigt_log")
  pigt_old <- getExportedValue("EMC2", "pigt_old")
  digt_hybrid <- getExportedValue("EMC2", "digt")
  digt_log <- getExportedValue("EMC2", "digt_log")
  digt_old <- getExportedValue("EMC2", "digt_old")

  pigt_exact <- function(t, k, l, a, threshold = 1e-10) {
    if (t <= 0) return(0)
    if (a < threshold) {
      mu <- k / l
      lambda <- k * k
      p1 <- pnorm(sqrt(lambda / t) * (1 + t / mu), lower.tail = FALSE)
      p2 <- pnorm(sqrt(lambda / t) * (1 - t / mu), lower.tail = FALSE)
      return(exp(exp(log(2 * lambda) - log(mu)) + log(p1)) + p2)
    }

    sqt <- sqrt(t)
    lgt <- log(t)
    if (l < threshold) {
      t5a <- 2 * pnorm((k + a) / sqt, lower.tail = TRUE) - 1
      t5b <- 2 * pnorm((-k - a) / sqt, lower.tail = TRUE) - 1
      t6a <- -0.5 * ((k + a) * (k + a) / t - log(2) - log(pi) + lgt) - log(a)
      t6b <- -0.5 * ((k - a) * (k - a) / t - log(2) - log(pi) + lgt) - log(a)
      cdf <- 1 + exp(t6a) - exp(t6b) + ((-k + a) * t5a - (k - a) * t5b) / (2 * a)
    } else {
      t1a <- exp(-0.5 * (k - a - t * l)^2 / t)
      t1b <- exp(-0.5 * (a + k - t * l)^2 / t)
      t1 <- exp(0.5 * (lgt - log(2) - log(pi))) * (t1a - t1b)
      t2a <- exp(2 * l * (k - a) + pnorm(- (k - a + t * l) / sqt, lower.tail = TRUE, log.p = TRUE))
      t2b <- exp(2 * l * (k + a) + pnorm(- (k + a + t * l) / sqt, lower.tail = TRUE, log.p = TRUE))
      t2 <- a + (t2b - t2a) / (2 * l)
      t4a <- 2 * pnorm((k + a) / sqt - sqt * l, lower.tail = TRUE) - 1
      t4b <- 2 * pnorm((k - a) / sqt - sqt * l, lower.tail = TRUE) - 1
      t4 <- 0.5 * (t * l - a - k + 0.5 / l) * t4a + 0.5 * (k - a - t * l - 0.5 / l) * t4b
      cdf <- 0.5 * (t4 + t2 + t1) / a
    }
    if (cdf < 0 || is.nan(cdf)) 0 else cdf
  }

  digt_exact <- function(t, k, l, a, threshold = 1e-10) {
    if (t <= 0) return(0)
    if (a < threshold) {
      lambda <- k * k
      if (l == 0) {
        e <- -0.5 * lambda / t
      } else {
        mu <- k / l
        e <- - (lambda / (2 * t)) * ((t * t) / (mu * mu) - 2 * t / mu + 1)
      }
      return(exp(e + 0.5 * log(lambda) - 0.5 * log(2 * t * t * t * pi)))
    }

    if (l < threshold) {
      term <- exp(-(k - a) * (k - a) / (2 * t)) - exp(-(k + a) * (k + a) / (2 * t))
      pdf <- exp(-0.5 * (log(2) + log(pi) + log(t)) + log(term) - log(2) - log(a))
    } else {
      sqt <- sqrt(t)
      t1a <- -(a - k + t * l) * (a - k + t * l) / (2 * t)
      t1b <- -(a + k - t * l) * (a + k - t * l) / (2 * t)
      t1 <- 2^(-1/2) * (exp(t1a) - exp(t1b)) / (sqrt(pi) * sqt)
      t2a <- 2 * pnorm((-k + a) / sqt + sqt * l, lower.tail = TRUE) - 1
      t2b <- 2 * pnorm((k + a) / sqt - sqt * l, lower.tail = TRUE) - 1
      t2 <- exp(log(0.5) + log(l)) * (t2a + t2b)
      pdf <- exp(log(t1 + t2) - log(2) - log(a))
    }
    if (pdf < 0 || is.nan(pdf)) 0 else pdf
  }

  eval_grid <- function(f, grid) {
    vapply(grid, function(x) do.call(f, x), numeric(1))
  }

  check_grid <- function(f_hybrid, f_log, f_old, f_exact, grid, label) {
    hybrid <- eval_grid(f_hybrid, grid)
    logv <- eval_grid(f_log, grid)
    old <- eval_grid(f_old, grid)
    exact <- eval_grid(f_exact, grid)
    mat <- cbind(hybrid = hybrid, log = logv, old = old, exact = exact)
    list(
      label = label,
      max_abs_hybrid_exact = max(abs(mat[, "hybrid"] - mat[, "exact"])),
      max_abs_log_exact = max(abs(mat[, "log"] - mat[, "exact"])),
      max_abs_old_exact = max(abs(mat[, "old"] - mat[, "exact"])),
      worst = mat[which.max(abs(mat[, "hybrid"] - mat[, "exact"])), , drop = FALSE]
    )
  }

  grid <- list(
    list(t = 0.3241393266620424, k = 15.10219804282587, l = 0.2940211321757806, a = 0.4204380070021307),
    list(t = 0.12, k = 1.4, l = 0.2, a = 0.35),
    list(t = 0.4, k = 2.0, l = 0.9, a = 0.25),
    list(t = 1.2, k = 3.0, l = -0.8, a = 0.4),
    list(t = 2.0, k = 1.1, l = 1.7, a = 0.15)
  )

  p_res <- check_grid(pigt_hybrid, pigt_log, pigt_old, pigt_exact, grid, "pigt")
  d_res <- check_grid(digt_hybrid, digt_log, digt_old, digt_exact, grid, "digt")

  list(pigt = p_res, digt = d_res)
}

if (sys.nframe() == 0) {
  res <- compare_wald_logspace()
  print(res)
}
