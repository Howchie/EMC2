test_that("prdmswtn matches direct threshold quadrature after integration swap", {
  skip_if_not_installed("statmod")

  t <- 0.85
  b <- 1.6
  v <- 1.1
  A <- 0.35
  sv <- 0.55
  s <- 1.0
  c <- 0.15
  lambda <- 0.25

  gl <- statmod::gauss.quad(80, "legendre")
  center <- b - 0.5 * A
  half_width <- 0.5 * A

  direct <- sum(gl$weights * vapply(gl$nodes, function(node) {
    threshold <- center + half_width * node
    EMC2:::pswtn(t, v, threshold, sv = sv, s = s, c = c, lambda = lambda)
  }, numeric(1))) * half_width / A

  expect_equal(
    EMC2:::prdmswtn(t, v, b, A, sv = sv, s = s, c = c, lambda = lambda, n_gauss_nodes = 20),
    direct,
    tolerance = 1e-8
  )
})
