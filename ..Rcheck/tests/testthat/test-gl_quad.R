# Pins the shared Gauss-Legendre cache (src/gl_quad.h::gl_get_rule) that all
# C++ quadratures now read (RDMSWTN/GBM hit-mass, logical-rules 31-node pass,
# global-kill omission mass, SS stop-success GL route). Reference values come
# from statmod::gauss.quad, which the removed implementations
# (get_gl_nodes_weights / hardcoded gl20 / fill_legendre_rule) matched.

test_that("gl_get_rule matches statmod::gauss.quad", {
  skip_if_not_installed("statmod")
  for (n in c(20L, 31L, 64L)) {
    got <- EMC2:::gl_rule_nodes_weights(n)
    ref <- statmod::gauss.quad(n, "legendre")
    # absolute comparison: odd rules contain a node at exactly 0, where
    # relative tolerance is meaningless
    expect_lt(max(abs(got$nodes - ref$nodes)), 1e-14, label = paste0("nodes n=", n))
    expect_lt(max(abs(got$weights - ref$weights)), 1e-14, label = paste0("weights n=", n))
  }
})

test_that("gl_get_rule rules are symmetric and integrate exactly", {
  for (n in c(20L, 31L)) {
    r <- EMC2:::gl_rule_nodes_weights(n)
    expect_equal(sum(r$weights), 2, tolerance = 1e-14)          # \int_{-1}^{1} 1 dx
    expect_equal(r$nodes, -rev(r$nodes), tolerance = 1e-15)     # antisymmetric nodes
    expect_equal(r$weights, rev(r$weights), tolerance = 1e-15)  # symmetric weights
    # n-point GL is exact for polynomials up to degree 2n-1
    expect_equal(sum(r$weights * r$nodes^4), 2 / 5, tolerance = 1e-14)
    expect_equal(sum(r$weights * r$nodes^(2 * n - 2)), 2 / (2 * n - 1), tolerance = 1e-13)
  }
})
