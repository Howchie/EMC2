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

test_that("BAwL retains its defective upper tail when k > 0", {
  p <- EMC2:::pleakyba_norm(Inf, A = 0.3, b = 1.5, v = 1, sv = 1,
                            k = 0.4, posdrift = TRUE)
  expect_gt(p, 0)
  expect_lt(p, 1)
})

test_that("LBA exposes no BAwL leak or clock parameters", {
  expect_identical(LBA()$c_name, "LBA")
  expect_false(any(c("k", "mG", "mK") %in% names(LBA()$p_types)))
  expect_false(any(c("k", "mG", "mK") %in% LBA()$p_types_canonical))
})
