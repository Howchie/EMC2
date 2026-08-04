# Extracted from test-rlf-fht.R:95

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
r <- EMC2:::rlf_fht_pdf_cdf_vec(
    seq(0.01, 2, length.out = 201),
    1, 1, 1.7, 1, 0, 80L, 100L
  )
expect_gt(r$domain_refinements, 0L)
expect_gt(r$spatial_refinements, 0L)
expect_lt(r$x_lo, -12)
expect_gt(r$nx_used, 80L)
expect_gt(r$nt_used, 100L)
expect_lte(r$domain_pdf_error, 0.01)
expect_lte(r$domain_cdf_error, 0.002)
expect_lte(r$spatial_pdf_error, 0.05)
expect_lte(r$spatial_cdf_error, 0.0125)
expect_true(all(diff(r$cdf) >= -1e-10))
expect_true(all(r$cdf >= 0 & r$cdf <= 1 + 1e-12))
expect_true(all(r$pdf >= 0))
expect_gte(r$min_density, -1e-10)
expect_lt(r$mismatch, 1e-10)
