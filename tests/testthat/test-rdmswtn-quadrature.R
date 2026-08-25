skip_model_validation()

test_that("infinite-time SWTN/RDMSWTN masses match analytic hit-mass identities", {
  mu <- -0.4
  sv <- 0.5
  s <- 1.0
  b <- 1.3
  A <- 0.35
  v <- 1.1

  H_point <- function(mu, sv, s, b) {
    alpha <- 2 * b / s^2
    pnorm(mu / sv) +
      exp(alpha * mu + 0.5 * alpha^2 * sv^2) *
      pnorm(-(mu + alpha * sv^2) / sv)
  }

  expect_equal(
    EMC2:::pswtn(Inf, mu, b, sv = sv, s = s, posdrift = FALSE),
    H_point(mu, sv, s, b),
    tolerance = 1e-10
  )

  H_spv <- integrate(function(r) H_point(mu, sv, s, r),
                     lower = b - A, upper = b)$value / A

  expect_equal(
    EMC2:::prdmswtn(Inf, mu, b, A, sv = sv, s = s, posdrift = FALSE),
    H_spv,
    tolerance = 1e-8
  )

  expect_lt(
    EMC2:::prdmswtn(Inf, v, b, A, sv = sv, s = s,
                    lambda_k = 0.5, posdrift = TRUE),
    1
  )

  expect_equal(
    EMC2:::prdmswtn(Inf, v, b, A, sv = sv, s = s,
                    lambda_k = 0, posdrift = TRUE),
    1,
    tolerance = 1e-12
  )
})
