test_that("FRQ exposes h and lambda coordinates", {
  m <- EMC2::FRQ()
  expect_equal(m$p_types_canonical,
               c("alpha", "beta", "h", "lambda", "t0"))
  expect_equal(names(m$p_types)[seq_len(7)],
               c("alpha", "beta", "h", "lambda", "t0", "delta", "cv_u"))
  expect_equal(unname(m$bound$exception[["h"]]), 1)
})

test_that("FRQ direct coordinates agree with the closed form", {
  alpha <- 2.5
  shape_b <- 1.7
  p <- 0.83
  h <- pbeta(p, alpha, shape_b)
  lambda <- 1.4
  x <- c(0.01, 0.1, 0.5, 2)
  q <- p * (-expm1(-lambda * x))
  want_p <- pbeta(q, alpha, shape_b)
  want_d <- p * lambda * exp(-lambda * x) *
    q^(alpha - 1) * (1 - q)^(shape_b - 1) / beta(alpha, shape_b)
  expect_equal(EMC2:::pfrq(x, alpha, shape_b, h, lambda), want_p,
               tolerance = 1e-11)
  expect_equal(EMC2:::dfrq(x, alpha, shape_b, h, lambda), want_d,
               tolerance = 1e-11)
  expect_equal(EMC2:::pfrq(Inf, alpha, shape_b, h, lambda),
               pbeta(p, alpha, shape_b), tolerance = 1e-11)
})

test_that("h = 1 is an exact non-defective FRQ boundary", {
  pars <- cbind(alpha = 2, beta = 3, h = 1, lambda = 0.7, t0 = 0)
  expect_equal(EMC2:::pFRQ(Inf, pars), 1, tolerance = 1e-12)
  expect_equal(EMC2:::sFRQ(Inf, pars), 0, tolerance = 1e-12)
  mapped <- EMC2::FRQ()$Ttransform(pars, NULL)
  expect_equal(unname(mapped[, "p"]), 1, tolerance = 1e-12)
})

test_that("the proper FRQ boundary stays exact with threshold variability", {
  mapped <- EMC2:::frq_rate(
    alpha = c(2, 5), beta = c(3, 2), h = 1, lambda = 1.3,
    delta = c(0.5, 2)
  )
  expect_equal(unname(mapped[, "p"]), c(1, 1), tolerance = 0)

  near <- EMC2:::frq_rate(
    alpha = c(2, 5), beta = c(3, 2), h = 1 - 1e-9, lambda = 1.3,
    delta = c(0.5, 2)
  )
  expect_true(all(near[, "p"] < 1))
})
