# BTAwL: ballistic Smith transient drive with state leak, and its shared-
# strength sustained/transient extension.

btawl_rows <- function(n, launch = 0L, A = 0.3, k = 0.5) {
  if (launch == 1L) {
    data.frame(mu = rep(0.2, n), sigma = rep(0.5, n),
               b = rep(1.3, n), A = rep(A, n), t0 = rep(0, n),
               k = rep(k, n), tau = rep(1, n))
  } else {
    data.frame(v = rep(1, n), sv = rep(1, n),
               b = rep(1.3, n), A = rep(A, n), t0 = rep(0, n),
               k = rep(k, n), tau = rep(1, n))
  }
}

btawl_mix_rows <- function(n, launch = 0L, pi = 0.5, k = 0.5) {
  p <- btawl_rows(n, launch = launch, k = k)
  p$tau_s <- 2
  p$tau_t <- 1
  p$pi <- pi
  p
}

test_that("the shared-strength mixture is exactly transient BTAwL at pi = 0", {
  tt <- c(0.5, 1, 1.5, 2, 2.5, 3, Inf)
  for (launch in 0:1) {
    p0 <- btawl_rows(length(tt), launch = launch)
    pm <- btawl_mix_rows(length(tt), launch = launch, pi = 0)
    expect_equal(
      EMC2:::pBTAwLMix(tt, pm, launch = launch),
      EMC2:::pBTAwL(tt, p0, launch = launch), tolerance = 1e-12
    )
    expect_equal(
      EMC2:::dBTAwLMix(tt, pm, launch = launch),
      EMC2:::dBTAwL(tt, p0, launch = launch), tolerance = 1e-12
    )
  }
})

test_that("transient BTAwL has the analytic hard endpoint for k > 0", {
  tm <- EMC2:::btawl_tmax_vec(0.5, 1)
  expect_true(is.finite(tm))
  tt <- c(tm * 0.8, tm * 0.99, tm, tm * 1.2, Inf)
  p <- btawl_rows(length(tt))
  cdf <- EMC2:::pBTAwL(tt, p)
  expect_true(all(diff(cdf) >= -1e-12))
  expect_equal(cdf[4], cdf[5], tolerance = 1e-12)
  expect_equal(EMC2:::dBTAwL(tt[4], p[4, ]), 0, tolerance = 1e-14)
})

test_that("pi = 1 handles the sustained-only member", {
  tt <- c(1, 2, 3, 5, 10, 20, Inf)
  p <- btawl_mix_rows(length(tt), pi = 1)
  cdf <- EMC2:::pBTAwLMix(tt, p)
  expect_true(all(diff(cdf) >= -1e-10))
  # With k > 0 the sustained channel tends to a finite asymptote, but it does
  # not acquire the transient hard endpoint.
  expect_gt(cdf[6], cdf[5])
  expect_gt(cdf[7], cdf[6])
  expect_gt(EMC2:::dBTAwLMix(10, p[5, ]), 0)
})

test_that("mixed CDF keeps the transient running minimum on a later decline", {
  # A short transient can overshoot the lower sustained asymptote.  Once the
  # transient opportunity has passed, the CDF must stay at its attained value
  # until the slower sustained channel adds new first passages.
  tt <- c(0.2, 0.3, 0.5, 1)
  p <- btawl_mix_rows(length(tt), pi = 0.01, k = 0.5)
  p[, "b"] <- 0.01
  p[, "tau_s"] <- 1
  p[, "tau_t"] <- 0.05
  cdf <- EMC2:::pBTAwLMix(tt, p)
  expect_true(all(diff(cdf) >= -1e-12))
  expect_equal(cdf[3], cdf[4], tolerance = 1e-10)
})

test_that("the no-leak sustained limit is retained", {
  tt <- c(1, 2, 5, Inf)
  p <- btawl_mix_rows(length(tt), pi = 1, k = 0)
  cdf <- EMC2:::pBTAwLMix(tt, p)
  expect_true(all(diff(cdf) >= -1e-12))
  expect_equal(cdf[length(cdf)], 1, tolerance = 1e-12)
})

test_that("constructors expose the transient and shared-strength charts", {
  old <- BTAwL()
  mix <- BTAwL_mixed()
  expect_equal(old$c_name, "BTAwL")
  expect_equal(mix$c_name, "BTAwL_MIX")
  expect_true(all(c("tau_s", "tau_t", "pi") %in% mix$p_types_canonical))
  expect_equal(BTAwL(mixture = TRUE)$c_name, "BTAwL_MIX")
})

test_that("BTAwL simulators select a winner per trial", {
  set.seed(42)
  lR <- factor(rep(c("left", "right"), 4))
  pars <- data.frame(v = rep(3, 8), sv = rep(0.05, 8),
                     B = rep(1, 8), A = rep(0.2, 8), t0 = rep(0.1, 8),
                     k = rep(0.5, 8), tau = rep(1, 8), b = rep(1.2, 8))
  out <- EMC2:::rBTAwL(lR, pars)
  expect_equal(nrow(out), 4)
  expect_true(all(is.finite(out$rt)))
})
