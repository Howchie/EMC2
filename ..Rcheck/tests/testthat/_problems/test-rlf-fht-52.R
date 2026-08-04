# Extracted from test-rlf-fht.R:52

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
tg2 <- seq(0.1, 2.0, by = 0.05)
mismatch <- vapply(c(500L, 1000L, 2000L), function(nt) {
    EMC2:::rlf_fht_pdf_cdf_vec(tg2, 1, 1, 2, 1, 0, 250L, nt)$mismatch
  }, numeric(1))
expect_lt(mismatch[1], 1e-4)
expect_gt(mismatch[1] / mismatch[2], 3.5)
