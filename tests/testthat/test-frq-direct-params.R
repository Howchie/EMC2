test_that("FRQ exposes the counter coordinates", {
  m <- EMC2::FRQ()
  expect_equal(m$p_types_canonical, c("K", "v", "sv", "h", "t0"))
  expect_equal(names(m$p_types)[seq_len(7)],
               c("K", "v", "sv", "h", "t0", "delta", "clock"))
  expect_equal(unname(m$bound$exception[["h"]]), 1)
})

test_that("FRQ counter coordinates agree with the generative closed form", {
  K <- 2.5; shape_b <- 1.7; p <- 0.83; lambda <- 1.4
  h <- pbeta(p, K, shape_b)
  v <- shape_b * p * lambda
  sv <- v / sqrt(shape_b)
  x <- c(0.01, 0.1, 0.5, 2)
  q <- p * (-expm1(-lambda * x))
  want_p <- pbeta(q, K, shape_b)
  want_d <- p * lambda * exp(-lambda * x) *
    q^(K - 1) * (1 - q)^(shape_b - 1) / beta(K, shape_b)
  expect_equal(EMC2:::pfrq(x, K, v, sv, h, clock = 0), want_p, tolerance = 1e-11)
  expect_equal(EMC2:::dfrq(x, K, v, sv, h, clock = 0), want_d, tolerance = 1e-11)
  expect_equal(EMC2:::pfrq(Inf, K, v, sv, h, clock = 0), h, tolerance = 1e-11)
})

test_that("h = 1 is an exact non-defective FRQ boundary", {
  pars <- cbind(K = 2, v = 2.1, sv = 1.218, h = 1, t0 = 0)
  expect_equal(EMC2:::pFRQ(Inf, pars), 1, tolerance = 1e-12)
  expect_equal(EMC2:::sFRQ(Inf, pars), 0, tolerance = 1e-12)
  mapped <- EMC2::FRQ()$Ttransform(pars, NULL)
  expect_equal(unname(mapped[, "p"]), 1, tolerance = 1e-12)
})

test_that("the proper FRQ boundary stays exact with threshold variability", {
  mapped <- EMC2:::frq_rate(K = c(2, 5), v = 1.3, sv = c(0.754, 0.91), h = 1,
                            delta = c(0.5, 2))
  expect_equal(unname(mapped[, "p"]), c(1, 1), tolerance = 0)
  near <- EMC2:::frq_rate(K = c(2, 5), v = 1.3, sv = c(0.754, 0.91),
                          h = 1 - 1e-9, delta = c(0.5, 2))
  expect_true(all(near[, "p"] < 1))
})
