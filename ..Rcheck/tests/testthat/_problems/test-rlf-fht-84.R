# Extracted from test-rlf-fht.R:84

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
