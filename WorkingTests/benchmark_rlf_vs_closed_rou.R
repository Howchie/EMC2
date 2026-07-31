# Benchmark the batched RLF RT kernel against the closed-form RDM/Wald
# calculation and the ROU Fokker--Planck kernel.
#
# Run from the package root after installing EMC2, for example:
#
#   Rscript WorkingTests/benchmark_rlf_vs_closed_rou.R
#
# The scaling benchmark has two regimes:
#   shared_key   -- all rows share one numerical solve (cache amortisation)
#   distinct_key -- every row has a distinct parameter tuple (particle-like)
#
# The discrepancy benchmark uses the common Brownian limit RLF(alpha = 2) and
# ROU(k = 0), with the RDM/Wald result as the closed-form oracle.  Set
# EMC2_RLF_BENCH_OUT to a directory to write CSV files as well as printing the
# results, e.g.:
#
#   EMC2_RLF_BENCH_OUT=/tmp/rlf-bench \
#     Rscript WorkingTests/benchmark_rlf_vs_closed_rou.R

if (!requireNamespace("EMC2", quietly = TRUE)) {
  stop("Install the EMC2 package before running this benchmark.")
}
library(EMC2)

particle_counts <- c(1L, 16L, 256L, 1024L)

# These match the shipped model defaults used for the cost comparison.
rlf_grid <- list(
  nx = 160L, dt = 0.016, tgrade = 1,
  richardson = TRUE, richardson_ratio = 1.25
)
rou_grid <- list(nx = 512L, dt = 0.004, grade = 8, tgrade = 32)

make_inputs <- function(n, distinct = FALSE) {
  i <- seq_len(n)
  v <- if (distinct) 1.35 + 0.0008 * i else rep(1.5, n)
  B <- if (distinct) 0.95 + 0.0004 * (i %% 11) else rep(1, n)
  list(
    rt = rep(0.8, n),
    v = v,
    B = B,
    A = rep(0.2, n),
    t0 = rep(0.1, n),
    s = rep(1, n),
    alpha = rep(2, n)
  )
}

call_rlf <- function(x, grid = rlf_grid) {
  EMC2:::rlf_pdf_cdf_vec(
    x$rt, x$v, x$B, x$A, x$t0, x$s, x$alpha,
    nx = grid$nx, dt_target = grid$dt, tgrade = grid$tgrade,
    adaptive = FALSE, explicit_inverse = TRUE,
    sparse_output = TRUE, simd_batch = FALSE,
    richardson = grid$richardson,
    richardson_ratio = grid$richardson_ratio
  )
}

call_rou <- function(x, grid = rou_grid) {
  EMC2:::rou_pdf_cdf_vec(
    x$rt, x$v, rep(0, length(x$rt)), x$B, x$A, x$t0, x$s,
    nx = grid$nx, dt_target = grid$dt, grade = grid$grade,
    tgrade = grid$tgrade
  )
}

# dWald/pWald are the vectorised C++ closed-form reference.  Calling both is
# comparable with the RLF/ROU wrappers, which return PDF and CDF together from
# one numerical solve.
call_closed_form <- function(x) {
  list(
    pdf = EMC2:::dWald(x$rt, x$v, x$B, x$A, x$t0, FALSE),
    cdf = EMC2:::pWald(x$rt, x$v, x$B, x$A, x$t0, FALSE)
  )
}

elapsed_ms <- function(fun, repetitions) {
  values <- vapply(seq_len(3L), function(...) {
    invisible(gc(FALSE))
    started <- Sys.time()
    for (j in seq_len(repetitions)) invisible(fun())
    as.numeric(difftime(Sys.time(), started, units = "secs")) * 1000 /
      repetitions
  }, numeric(1))
  unname(stats::median(values))
}

scaling_benchmark <- function() {
  rows <- list()
  for (distinct in c(FALSE, TRUE)) {
    regime <- if (distinct) "distinct_key" else "shared_key"
    for (n in particle_counts) {
      x <- make_inputs(n, distinct = distinct)
      pde_repetitions <- if (n <= 16L) 10L else if (n <= 256L) 3L else 1L
      closed_repetitions <- if (n <= 16L) 1000L else if (n <= 256L) 500L else 200L

      rlf_ms <- elapsed_ms(function() call_rlf(x), pde_repetitions)
      rou_ms <- elapsed_ms(function() call_rou(x), pde_repetitions)
      closed_ms <- elapsed_ms(function() call_closed_form(x), closed_repetitions)

      rlf_check <- call_rlf(x)
      rou_check <- call_rou(x)
      rows[[length(rows) + 1L]] <- data.frame(
        regime = regime,
        n = n,
        RLF_ms = rlf_ms,
        ROU_ms = rou_ms,
        closed_form_ms = closed_ms,
        RLF_solves = rlf_check$n_solves,
        ROU_solves = rou_check$n_solves,
        RLF_us_per_row = 1000 * rlf_ms / n,
        ROU_us_per_row = 1000 * rou_ms / n,
        closed_form_us_per_row = 1000 * closed_ms / n,
        RLF_over_closed = rlf_ms / closed_ms,
        ROU_over_closed = rou_ms / closed_ms,
        row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}

error_metrics <- function(x, oracle) {
  c(
    max_abs_pdf = max(abs(x$pdf - oracle$pdf)),
    rms_pdf = sqrt(mean((x$pdf - oracle$pdf)^2)),
    max_abs_cdf = max(abs(x$cdf - oracle$cdf)),
    rms_cdf = sqrt(mean((x$cdf - oracle$cdf)^2))
  )
}

discrepancy_benchmark <- function() {
  n <- 200L
  x <- make_inputs(n)
  x$rt <- seq(0.11, 3, length.out = n)
  oracle <- call_closed_form(x)

  rlf <- call_rlf(x)
  rou <- call_rou(x)
  data.frame(
    comparison = c("RLF_alpha2_vs_RDM", "ROU_k0_vs_RDM"),
    rbind(error_metrics(rlf, oracle), error_metrics(rou, oracle)),
    row.names = NULL
  )
}

resolution_benchmark <- function() {
  n <- 200L
  x <- make_inputs(n)
  x$rt <- seq(0.11, 3, length.out = n)
  oracle <- call_closed_form(x)
  rows <- list()

  for (grid in list(
    list(model = "RLF", nx = 100L, dt = 0.01),
    list(model = "RLF", nx = 128L, dt = 0.004),
    list(model = "RLF", nx = 200L, dt = 0.005),
    list(model = "RLF", nx = 400L, dt = 0.0025),
    list(model = "RLF", nx = 800L, dt = 0.00125),
    list(model = "ROU", nx = 256L, dt = 0.008),
    list(model = "ROU", nx = 512L, dt = 0.004),
    list(model = "ROU", nx = 1024L, dt = 0.002)
  )) {
    if (grid$model == "RLF") {
      fun <- function() call_rlf(x, list(
        nx = grid$nx, dt = grid$dt, tgrade = 1,
        richardson = FALSE, richardson_ratio = 1.25
      ))
    } else {
      fun <- function() call_rou(x, list(
        nx = grid$nx, dt = grid$dt, grade = 8, tgrade = 32
      ))
    }
    result <- fun()
    errors <- error_metrics(result, oracle)
    rows[[length(rows) + 1L]] <- do.call(
      data.frame,
      c(
        list(
          model = grid$model,
          nx = grid$nx,
          dt = grid$dt,
          batch_ms = elapsed_ms(fun, 2L)
        ),
        as.list(errors),
        list(solves = if (!is.null(result$n_solves)) result$n_solves else NA_integer_),
        list(row.names = NULL)
      )
    )
  }
  do.call(rbind, rows)
}

scaling <- scaling_benchmark()
discrepancy <- discrepancy_benchmark()
resolution <- resolution_benchmark()

cat("SCALING (one PDF+CDF batch)\n")
print(scaling, row.names = FALSE, digits = 5)
cat("\nBROWNIAN-LIMIT DISCREPANCY\n")
print(discrepancy, row.names = FALSE, digits = 7)
cat("\nRESOLUTION SWEEP (200 RT queries)\n")
print(resolution, row.names = FALSE, digits = 7)

output_dir <- Sys.getenv("EMC2_RLF_BENCH_OUT", "")
if (nzchar(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(scaling, file.path(output_dir, "scaling.csv"), row.names = FALSE)
  utils::write.csv(discrepancy, file.path(output_dir, "discrepancy.csv"), row.names = FALSE)
  utils::write.csv(resolution, file.path(output_dir, "resolution.csv"), row.names = FALSE)
  cat("\nWrote CSV files to ", normalizePath(output_dir), "\n", sep = "")
}
