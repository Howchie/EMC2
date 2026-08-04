# Extracted from test-rlf-fht.R:83

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
    tg, 1, 1, 1.7, 1, 0, 250L, 500L
  )
expect_false(r$refinement_checked)
expect_true(r$refinement_skipped)
expect_equal(r$domain_refinements, 0L)
expect_equal(r$spatial_refinements, 0L)
expect_true(is.nan(r$domain_cdf_error))
expect_true(is.nan(r$spatial_cdf_error))
