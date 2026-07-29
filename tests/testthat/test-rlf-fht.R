# Race Levy Flight (RLF) model (src/model_RLF.h, src/rlf_diffusion.cpp).
#
# Deterministic and simulation checks for the Option A SFFPE solver and CMS path simulator.

tg <- seq(0.1, 2.0, by = 0.05)

wald_pdf <- function(t, mu, b) b / sqrt(2 * pi * t^3) * exp(-(b - mu * t)^2 / (2 * t))
wald_cdf <- function(t, mu, b) {
  pnorm((mu * t - b) / sqrt(t)) + exp(2 * mu * b) * pnorm(-(b + mu * t) / sqrt(t))
}

test_that("RLF PDE solver reduces to Wald form when alpha == 2.0", {
  # When alpha = 2.0, the space-fractional operator reduces to standard Brownian motion.
  for (cs in list(c(1, 1, 1), c(2, 1, 1), c(0.5, 1, 1.5))) {
    v <- cs[1]; sigma <- cs[2]; b0 <- cs[3]
    r <- EMC2:::rlf_fht_pdf_cdf_vec(tg, v, sigma, 2.0, b0, 0, 250L, 500L)

    expect_lt(max(abs(r$pdf - wald_pdf(tg, v / sigma, b0 / sigma))), 0.25)
    expect_lt(max(abs(r$cdf - wald_cdf(tg, v / sigma, b0 / sigma))), 0.035)
  }
})

test_that("RLF PDE solver produces valid PDF and CDF for heavy-tailed noise (alpha = 1.7)", {
  r <- EMC2:::rlf_fht_pdf_cdf_vec(tg, 1.0, 1.0, 1.7, 1.0, 0.0, 200L, 400L)

  # CDF should be monotonically increasing and strictly inside [0, 1]
  expect_true(all(diff(r$cdf) >= -1e-6))
  expect_true(all(r$cdf >= 0 & r$cdf <= 1.0))
  expect_true(all(r$pdf >= 0))
})

test_that("RLF path simulator generates valid hit times and matches PDE solver CDF", {
  set.seed(42)
  sim_hits <- EMC2:::simulate_rlf_hit_times_cpp(20000L, 1.0, 1.0, 1.7, 1.0, 0.0, 2.0, 0.001, 42L)

  # Calculate empirical CDF at t=2.0
  emp_cdf_t2 <- mean(!is.na(sim_hits) & sim_hits <= 2.0)

  # PDE solver CDF at t=2.0
  r <- EMC2:::rlf_fht_pdf_cdf_vec(c(2.0), 1.0, 1.0, 1.7, 1.0, 0.0, 250L, 500L)
  pde_cdf_t2 <- r$cdf[1]

  # Should agree within Monte Carlo & discretization sampling margin (~5%)
  expect_lt(abs(emp_cdf_t2 - pde_cdf_t2), 0.06)
})
