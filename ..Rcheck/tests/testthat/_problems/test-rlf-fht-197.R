# Extracted from test-rlf-fht.R:197

# prequel ----------------------------------------------------------------------
tg <- seq(0.1, 2.0, by = 0.05)
wald_pdf <- function(t, mu, b) {
  b / sqrt(2 * pi * t^3) * exp(-(b - mu * t)^2 / (2 * t))
}
wald_cdf <- function(t, mu, b) {
  pnorm((mu * t - b) / sqrt(t)) +
    exp(2 * mu * b) * pnorm(-(b + mu * t) / sqrt(t))
}

# test -------------------------------------------------------------------------
hits <- EMC2:::simulate_rlf_hit_times_cpp(
    20000L, 1, 1, 1.7, 1, 0, 2, 0.001, 42L
  )
empirical <- mean(!is.na(hits) & hits <= 2)
pde <- EMC2:::rlf_fht_pdf_cdf_vec(
    2, 1, 1, 1.7, 1, 0, 250L, 500L
  )
expect_equal(empirical, pde$cdf[1], tolerance = 0.015)
expect_lt(pde$lower_boundary_pressure, 0.005)
