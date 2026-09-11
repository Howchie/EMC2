# C15: one factorisation per group draw, not one per subject per component.
#
# The prior density under blocking is evaluated against the same full
# non-marginal covariance for every subject and every component of every
# iteration, and each call factorised it again. The cached factor has to
# reproduce the original call exactly -- not to a tolerance -- or the posterior
# moves.

test_that("the hoisted factor reproduces fast_dmvnorm bit for bit", {
  set.seed(21)
  for (p in c(1L, 2L, 8L, 40L)) {
    A <- matrix(rnorm(p * p), p, p)
    sigma <- crossprod(A) + diag(p) * 0.5
    mu <- rnorm(p)
    x <- matrix(rnorm(17L * p), 17L, p)
    want <- EMC2:::fast_dmvnorm(x, mu, sigma)
    f <- EMC2:::fast_dmvnorm_factor(sigma)
    expect_true(f$ok, info = paste("p =", p))
    got <- EMC2:::fast_dmvnorm_rooti(x, mu, f$rooti, f$log_const)
    expect_bit_identical(got, want, paste("p =", p))
  }
})

test_that("a singular covariance keeps fast_dmvnorm's answer", {
  # It returns -Inf for every row rather than regularising. The proposal
  # regularisation ladder in .chol_factor is a different thing and must not
  # leak into the prior.
  p <- 4L
  sigma <- matrix(1, p, p)   # rank 1: Cholesky fails
  f <- EMC2:::fast_dmvnorm_factor(sigma)
  expect_false(f$ok)
  x <- matrix(rnorm(5L * p), 5L, p)
  want <- EMC2:::fast_dmvnorm(x, rep(0, p), sigma)
  expect_true(all(is.infinite(want) & want < 0))
  # And the proposal factoriser does NOT agree, which is the point of keeping
  # them apart: it regularises and returns a finite density.
  reg <- EMC2:::.chol_factor(sigma, p)
  expect_false(is.null(reg))
  expect_true(all(is.finite(reg$rooti)))
})

test_that("the cache carries the full factor and is keyed on the group draw", {
  set.seed(22)
  p <- 6L
  A <- matrix(rnorm(p * p), p, p)
  gv <- crossprod(A) + diag(p)
  idx_list <- list(rep(TRUE, p))
  marginal <- c(rep(FALSE, p - 1L), TRUE)
  cache <- EMC2:::build_group_chol_cache(gv, idx_list, marginal)
  expect_false(is.null(cache$full))
  expect_true(cache$full$ok)
  expect_identical(cache$full$keep, !marginal)
  expect_identical(dim(cache$full$rooti), c(p - 1L, p - 1L))
  # It reproduces the uncached call on the submatrix it was built for.
  keep <- !marginal
  x <- matrix(rnorm(9L * sum(keep)), 9L, sum(keep))
  mu <- rnorm(sum(keep))
  expect_bit_identical(
    EMC2:::fast_dmvnorm_rooti(x, mu, cache$full$rooti, cache$full$log_const),
    EMC2:::fast_dmvnorm(x, mu, gv[keep, keep, drop = FALSE]),
    "cached full factor")
  # Without a marginal index there is nothing to cache, and the component
  # factors are unaffected.
  plain <- EMC2:::build_group_chol_cache(gv, idx_list)
  expect_null(plain$full)
  expect_identical(plain$f, cache$f)
})

test_that("a singular group draw caches the refusal rather than a bad factor", {
  p <- 3L
  gv <- matrix(1, p, p)
  cache <- EMC2:::build_group_chol_cache(gv, list(rep(TRUE, p)),
                                         rep(FALSE, p))
  expect_false(cache$full$ok)
  # The component factor still regularises, because that one is for proposals.
  expect_false(is.null(cache$f[[1L]]))
})

test_that("the proposal factors are untouched by the addition", {
  # C15 adds to the cache; it must not change what was already in it, or every
  # proposal draw moves.
  set.seed(23)
  p <- 5L
  A <- matrix(rnorm(p * p), p, p)
  gv <- crossprod(A) + diag(p)
  idx_list <- list(c(TRUE, TRUE, FALSE, FALSE, FALSE),
                   c(FALSE, FALSE, TRUE, TRUE, TRUE))
  before <- lapply(idx_list, function(idx) EMC2:::.chol_factor(gv[idx, idx, drop = FALSE], sum(idx)))
  cache <- EMC2:::build_group_chol_cache(gv, idx_list, rep(FALSE, p))
  for (k in seq_along(idx_list)) {
    expect_bit_identical(cache$f[[k]]$rooti, before[[k]]$rooti,
                         paste("component", k))
    expect_bit_identical(cache$f[[k]]$sum_log_diag, before[[k]]$sum_log_diag,
                         paste("component", k, "log det"))
  }
})
