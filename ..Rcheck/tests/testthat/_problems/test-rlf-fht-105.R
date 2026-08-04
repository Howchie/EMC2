# Extracted from test-rlf-fht.R:105

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
expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(
      seq(0.001, 0.3, length.out = 101),
      10, 0.05, 1.99, 1, 0, 120L, 120L
    ),
    "automatic spatial refinement failed to converge"
  )
