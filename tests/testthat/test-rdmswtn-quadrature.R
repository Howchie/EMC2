test_that("prdmswtn matches direct threshold quadrature after integration swap", {
  skip_if_not_installed("statmod")

  t <- 0.85
  b <- 1.6
  v <- 1.1
  A <- 0.35
  sv <- 0.55
  s <- 1.0
  lambda <- 0.25

  gl <- statmod::gauss.quad(80, "legendre")
  center <- b - 0.5 * A
  half_width <- 0.5 * A

  # Integrate defective per-threshold CDFs, then normalise outside (posdrift=FALSE inside).
  direct_defective <- sum(gl$weights * vapply(gl$nodes, function(node) {
    threshold <- center + half_width * node
    EMC2:::pswtn(t, v, threshold, sv = sv, s = s, lambda_k = lambda, posdrift = FALSE)
  }, numeric(1))) * half_width / A

  expect_equal(
    EMC2:::prdmswtn(t, v, b, A, sv = sv, s = s, lambda_k = lambda,
                    n_gauss_nodes = 20, posdrift = FALSE),
    direct_defective,
    tolerance = 1e-6
  )

  # posdrift=TRUE: normalised result should equal defective / hit_mass
  p_pos  <- EMC2:::prdmswtn(t, v, b, A, sv = sv, s = s, lambda_k = lambda,
                             n_gauss_nodes = 20, posdrift = TRUE)
  p_inf  <- EMC2:::prdmswtn(Inf, v, b, A, sv = sv, s = s, lambda_k = 0,
                             n_gauss_nodes = 20, posdrift = FALSE)
  expect_equal(p_pos, direct_defective / p_inf, tolerance = 1e-5)
})

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
