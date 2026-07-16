probe <- local({
  env <- new.env(parent = globalenv())
  path <- normalizePath(testthat::test_path("..", "..", "tools",
                                            "adaptive-gk-probe.cpp"),
                        mustWork = TRUE)
  cache <- file.path(tempdir(), paste0("emc2-adaptive-gk-", Sys.getpid()))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(path, env = env, cacheDir = cache, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE)
  env
})

gaussian_log_integral <- function(center, width, amplitude = 0) {
  amplitude + log(width) + 0.5 * log(2 * pi) +
    log(pnorm((1 - center) / width) - pnorm(-center / width))
}

test_that("adaptive GK integrates a shared vector of log-space targets", {
  centers <- c(.2, .5, .8)
  widths <- c(.15, .05, .2)
  amplitudes <- c(0, -100, 3)
  got <- probe$adaptive_gk_gaussian(centers, widths, amplitudes,
                                    abs_tol = 1e-12, rel_tol = 1e-10,
                                    limit = 128, initial_breaks = numeric())
  expected <- gaussian_log_integral(centers, widths, amplitudes)

  expect_true(got$converged)
  expect_equal(got$log_integral, expected, tolerance = 1e-9)
  expect_equal(got$callback_points, got$eval_points, tolerance = 0)
  expect_gt(got$intervals, 1)
})

test_that("initial breaks expose narrow peaks to the adaptive rule", {
  center <- .947
  width <- .0025
  amplitude <- -50
  got <- probe$adaptive_gk_gaussian(
    center, width, amplitude,
    abs_tol = 1e-12, rel_tol = 1e-9, limit = 128,
    initial_breaks = c(.90, .94, .947, .954, .99))
  expected <- gaussian_log_integral(center, width, amplitude)

  expect_true(got$converged)
  expect_equal(got$log_integral, expected, tolerance = 1e-7)
})

test_that("the latent-factor candidate matches its analytic normal mixture", {
  means <- c(0, .2, -.1)
  loadings <- c(.5, 1, -.7)
  residual_sd <- c(.5, .2, .8)
  observations <- c(.1, .6, -.4)
  # The breaks are on x = Phi(z), not on z itself.  This gives the candidate
  # a fair chance to see tails while retaining one shared partition.
  breaks <- pnorm(seq(-6, 6, by = 1))
  got <- probe$adaptive_gk_latent_normal(
    means, loadings, residual_sd, observations,
    abs_tol = 1e-12, rel_tol = 1e-8, limit = 128,
    initial_breaks = breaks)
  expected <- dnorm(observations, means,
                    sqrt(loadings^2 + residual_sd^2), log = TRUE)

  expect_true(got$converged)
  expect_equal(got$log_integral, expected, tolerance = 1e-7)
  expect_equal(got$callback_points, got$eval_points, tolerance = 0)
})

test_that("invalid adaptive tolerances fail explicitly", {
  expect_error(
    probe$adaptive_gk_gaussian(.5, .1, 0, abs_tol = 0,
                               rel_tol = 1e-9, limit = 128,
                               initial_breaks = numeric()),
    "tolerances must be positive")
})
