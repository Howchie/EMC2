# Extracted from test-rlf-fht.R:30

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
for (cs in list(c(1, 1, 1), c(2, 1, 1), c(0.5, 1, 1.5))) {
    v <- cs[1]
    sigma <- cs[2]
    b0 <- cs[3]
    r <- EMC2:::rlf_fht_pdf_cdf_vec(
      tg, v, sigma, 2.0, b0, 0, 250L, 500L
    )

    expect_lt(
      max(abs(r$pdf - wald_pdf(tg, v / sigma, b0 / sigma))),
      0.18
    )
    expect_lt(
      max(abs(r$cdf - wald_cdf(tg, v / sigma, b0 / sigma))),
      0.022
    )
    expect_lt(r$mismatch, 1e-10)
  }
