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
