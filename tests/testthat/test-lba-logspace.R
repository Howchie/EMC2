call_norm <- function(fun, t, A, b, v, sv, k, posdrift = TRUE,
                      log_out = FALSE) {
  vapply(seq_along(t), function(i) {
    fun(t[i], A[i], b[i], v[i], sv[i], k[i], posdrift, log_out)
  }, numeric(1))
}

test_that("LBA log outputs remain finite past natural-scale underflow", {
  t <- c(0.02, 0.01)
  A <- c(0.3, 0.3)
  b <- c(1.5, 1.5)
  v <- c(1, 1)
  sv <- c(1, 1)

  density <- call_norm(EMC2:::dleakyba_norm, t, A, b, v, sv, rep(0, 2))
  log_density <- call_norm(EMC2:::dleakyba_norm, t, A, b, v, sv, rep(0, 2),
                           log_out = TRUE)
  cdf <- call_norm(EMC2:::pleakyba_norm, t, A, b, v, sv, rep(0, 2))
  log_cdf <- call_norm(EMC2:::pleakyba_norm, t, A, b, v, sv, rep(0, 2),
                       log_out = TRUE)

  expect_true(all(density == 0))
  expect_true(all(cdf == 0))
  expect_true(all(is.finite(log_density)))
  expect_true(all(is.finite(log_cdf)))
  expect_true(all(log_density < log(.Machine$double.xmin)))
  expect_true(all(log_cdf < log(.Machine$double.xmin)))
})

test_that("LBA natural and log outputs agree before underflow", {
  t <- c(0.2, 0.5, 1, 2)
  density <- call_norm(EMC2:::dleakyba_norm, t, rep(0.3, 4), rep(1.5, 4),
                       rep(0.3, 4), rep(1, 4), rep(0, 4))
  log_density <- call_norm(EMC2:::dleakyba_norm, t, rep(0.3, 4), rep(1.5, 4),
                           rep(0.3, 4), rep(1, 4), rep(0, 4), log_out = TRUE)
  cdf <- call_norm(EMC2:::pleakyba_norm, t, rep(0.3, 4), rep(1.5, 4),
                   rep(0.3, 4), rep(1, 4), rep(0, 4))
  log_cdf <- call_norm(EMC2:::pleakyba_norm, t, rep(0.3, 4), rep(1.5, 4),
                       rep(0.3, 4), rep(1, 4), rep(0, 4), log_out = TRUE)

  expect_equal(log_density, log(density), tolerance = 1e-12)
  expect_equal(log_cdf, log(cdf), tolerance = 1e-12)
})

test_that("LBA log density retains the sv term in the tail", {
  t <- 1.494864
  A <- b <- 0.6059518
  v <- -0.2953623
  sv <- 0.05147902
  z_hi <- (b - t * v) / (t * sv)
  z_lo <- (b - A - t * v) / (t * sv)
  denom <- pnorm(v / sv)
  reference <- (v * (pnorm(z_hi) - pnorm(z_lo)) +
                sv * (dnorm(z_lo) - dnorm(z_hi))) / (A * denom)

  expect_equal(
    EMC2:::dleakyba_norm(t, A, b, v, sv, 0, TRUE, TRUE),
    log(reference),
    tolerance = 1e-6
  )
})

test_that("BAwL log-space kernels extend to k > 0", {
  t <- c(0.15, 0.5, 1.5)
  args <- list(A = rep(0.3, 3), b = rep(1.5, 3), v = rep(1, 3),
               sv = rep(1, 3), k = rep(0.4, 3))
  density <- do.call(call_norm, c(list(fun = EMC2:::dleakyba_norm, t = t), args))
  log_density <- do.call(call_norm,
                         c(list(fun = EMC2:::dleakyba_norm, t = t,
                                log_out = TRUE), args))
  cdf <- do.call(call_norm, c(list(fun = EMC2:::pleakyba_norm, t = t), args))
  log_cdf <- do.call(call_norm,
                     c(list(fun = EMC2:::pleakyba_norm, t = t,
                            log_out = TRUE), args))

  expect_equal(log_density, log(density), tolerance = 1e-12)
  expect_equal(log_cdf, log(cdf), tolerance = 1e-12)
  expect_true(all(is.finite(log_density)))
  expect_true(all(is.finite(log_cdf)))
})

test_that("BAwL log-space kernels remain finite in the rare tail", {
  t <- c(0.02, 0.01)
  log_density <- vapply(t, function(tt) {
    EMC2:::dleakyba_norm(tt, 0.3, 1.5, 1, 1, 0.4, TRUE, TRUE)
  }, numeric(1))
  log_cdf <- vapply(t, function(tt) {
    EMC2:::pleakyba_norm(tt, 0.3, 1.5, 1, 1, 0.4, TRUE, TRUE)
  }, numeric(1))

  expect_true(all(is.finite(log_density)))
  expect_true(all(is.finite(log_cdf)))
  expect_true(all(log_density < log(.Machine$double.xmin)))
  expect_true(all(log_cdf < log(.Machine$double.xmin)))
})

test_that("LBA log CDF stays defined when both z endpoints are far in the upper tail", {
  # Large positive z = (b - t*v) / (t*sv): early times or strongly negative
  # drifts.  The Q-antiderivative's log-ratio rounds to >= 0 here and must
  # take the asymptotic branch instead of returning NA via log1m_exp.
  cases <- list(
    list(t = 1e-4, A = 0.3, b = 1.5, v = 3,    sv = 1),     # early time
    list(t = 0.5,  A = 0.3, b = 1.5, v = -1e6, sv = 1),     # extreme drift
    list(t = 0.5,  A = 0.3, b = 1.5, v = -1e6, sv = 1e-3),  # extreme drift, small sv
    list(t = 0.5,  A = 0.3, b = 1.5, v = -1e8, sv = 1)
  )
  for (a in cases) {
    log_cdf <- EMC2:::pleakyba_norm(a$t, a$A, a$b, a$v, a$sv, 0, TRUE, TRUE)
    expect_false(is.na(log_cdf))
    expect_lt(log_cdf, 0)
    cdf <- EMC2:::pleakyba_norm(a$t, a$A, a$b, a$v, a$sv, 0, TRUE, FALSE)
    expect_false(is.na(cdf))
  }
})

test_that("LBA log kernels survive a collapsed start range in the far tail", {
  # span = A / (t*sv) below ~1e-7 with far-tail z endpoints defeats the
  # antiderivative log-difference; the midpoint guard must take over, the
  # log CDF must stay defined and non-positive.
  cases <- list(
    list(t = 0.445, A = 2.32e-10, b = 0.0293, v = -26.6, sv = 0.0359),
    list(t = 7.52,  A = 3.14e-10, b = 0.0189, v = -24.1, sv = 0.037),
    list(t = 84.6,  A = 1.09e-6,  b = 42.9,   v = 5.8e6, sv = 0.775)
  )
  for (a in cases) {
    log_cdf <- EMC2:::pleakyba_norm(a$t, a$A, a$b, a$v, a$sv, 0, TRUE, TRUE)
    log_pdf <- EMC2:::dleakyba_norm(a$t, a$A, a$b, a$v, a$sv, 0, TRUE, TRUE)
    expect_false(is.na(log_cdf))
    expect_false(is.na(log_pdf))
    expect_lte(log_cdf, 0)
  }
})

test_that("BAwL retains its defective upper tail when k > 0", {
  p <- EMC2:::pleakyba_norm(Inf, A = 0.3, b = 1.5, v = 1, sv = 1,
                            k = 0.4, posdrift = TRUE)
  expect_gt(p, 0)
  expect_lt(p, 1)
})

test_that("scalar natural-scale kernels stay exact at truncation-bound regimes", {
  # Truncation normalisers evaluate the scalar CDF at bounds where it is
  # legitimately ~0 or ~1 and immediately clamp; the scalar path must return
  # those natural values (fast) rather than reject them, and must still agree
  # with the strict evaluators away from the boundaries.
  p_scalar <- function(t, v) EMC2:::pkilledleakyba(t, v = v, b = 3.4, A = 1.9,
    sv = 1, t0 = 0, k = 0, lambda_g = 0, lambda_k = 0)
  # Upper bound, strong drift: CDF saturates to ~1 (true survivor ~1e-15).
  expect_gte(p_scalar(20, 8), 1 - 1e-8)
  expect_lte(p_scalar(20, 8), 1)
  # Lower bound, strong negative drift: underflows to 0.
  expect_equal(p_scalar(0.05, -30), 0)
  # Mid regime agrees with the strict log evaluator.
  for (v in c(-2, 0.5, 3)) {
    expect_equal(p_scalar(0.8, v),
                 exp(EMC2:::pleakyba_norm(0.8, 1.9, 3.4, v, 1, 0, TRUE, TRUE)),
                 tolerance = 1e-10)
  }
})

test_that("extreme drifts stay cheap and finite in log tails", {
  # z beyond 37 used to fall back to R::pnorm; the asymptotic tail must give
  # finite, monotone log CDFs without it.
  lp <- vapply(c(-50, -100, -1000, -1e4), function(v) {
    EMC2:::pleakyba_norm(0.8, 1.9, 3.4, v, 1, 0, TRUE, TRUE)
  }, numeric(1))
  expect_true(all(is.finite(lp)))
  expect_true(all(diff(lp) < 0))
})

test_that("LBA exposes no BAwL leak or clock parameters", {
  expect_identical(LBA()$c_name, "LBA")
  expect_false(any(c("k", "mG", "mK") %in% names(LBA()$p_types)))
  expect_false(any(c("k", "mG", "mK") %in% LBA()$p_types_canonical))
})
