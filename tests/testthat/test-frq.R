skip_model_validation()

# FRQ, the finite reservoir quorum process -- src/model_FRQ.h, R/model_FRQ.R.
#
# The model is fully closed form, so strong oracles are cheap. This file keeps
# the ones that are fast enough for R CMD check:
#
#   1. A small simulation of the LITERAL generative process (N Bernoulli(p)
#      availability indicators, Exponential(lambda) latencies, K-th order
#      statistic). It shares no formula with the kernel, so it establishes
#      that the kernel's closed form is the right closed form.
#   2. An independent R transcription of the closed-form density and CDF.
#   3. The defective-distribution contract: F(Inf) = h, S(Inf) = 1 - h,
#      integrate(f) = h.
#   4. That the counter coordinates mean what they claim: h is the completion
#      probability, v the mean initial input rate, sv = v/sqrt(beta) the
#      between-trial rate SD, and (K, v, sv, h) inverts to the generating
#      (beta, p, lambda); sv = 0 is the exact Poisson-counter limit.
#
# The references below are written in the generative coordinates
# (alpha, beta, h, lambda); to_chart() maps them to the fitted ones.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

ref_pl <- function(a, b, h, lambda, delta = 0, clock = 0) {
  zh <- if (delta > 0) EMC2:::frq_h_inv_r(h, delta) else h
  p <- qbeta(zh, a, b)
  c(p = p, lambda = lambda)
}

ref_H <- function(z, d) {
  if (d == 0) return(z)
  lz <- log(z) - log1p(-z)
  integrate(function(e) 1 / (1 + exp(-(lz + e))), -d, d,
            rel.tol = 1e-12)$value / (2 * d)
}

ref_pH <- function(x, a, b, h, lambda, d = 0, clock = 0) {
  pl <- ref_pl(a, b, h, lambda, d, clock)
  p <- pl[["p"]]
  c2 <- clock^2
  q <- if (c2 == 0) {
    p * (-expm1(-lambda * x))
  } else {
    p * (1 - (1 + c2 * lambda * x)^(-1 / c2))
  }
  z <- pbeta(q, a, b)
  if (d > 0) vapply(z, ref_H, numeric(1), d = d) else z
}

ref_d <- function(x, a, b, h, lambda, delta = 0, clock = 0) {
  pl <- ref_pl(a, b, h, lambda, delta, clock)
  p <- pl[["p"]]
  c2 <- clock^2
  if (c2 == 0) {
    q <- p * (-expm1(-lambda * x))
    omq <- (1 - p) + p * exp(-lambda * x)
    dens <- p * lambda * exp(-lambda * x) *
      q^(a - 1) * omq^(b - 1) / beta(a, b)
  } else {
    sr <- (1 + c2 * lambda * x)^(-1 / c2)
    g <- lambda * (1 + c2 * lambda * x)^(-1 / c2 - 1)
    q <- p * (1 - sr)
    omq <- (1 - p) + p * sr
    dens <- p * g * q^(a - 1) * omq^(b - 1) / beta(a, b)
  }
  if (delta > 0) {
    z <- pbeta(q, a, b)
    c2_d <- 4 * sinh(delta / 2)^2
    sinhc <- if (delta < 0.1) 1 + delta^2 / 6 + delta^4 / 120 else sinh(delta) / delta
    Hprime <- sinhc / (1 + 2 * z * (1 - z) * (cosh(delta) - 1))
    dens <- dens * Hprime
  }
  dens
}

ref_p <- function(x, a, b, h, lambda, delta = 0, clock = 0) {
  ref_pH(x, a, b, h, lambda, d = delta, clock = clock)
}

ref_s <- function(x, a, b, h, lambda, delta = 0, clock = 0) {
  1 - ref_p(x, a, b, h, lambda, delta = delta, clock = clock)
}

to_chart <- function(a, b, h, lambda, delta = 0) {
  p <- ref_pl(a, b, h, lambda, delta)[["p"]]
  c(K = a, v = b * p * lambda, sv = p * lambda * sqrt(b))
}

cpp_d <- function(x, a, b, h, lambda, delta = 0, clock = 0, log_out = FALSE) {
  ch <- to_chart(a, b, h, lambda, delta)
  EMC2:::dfrq(x, ch[["K"]], ch[["v"]], ch[["sv"]], h, delta = delta,
              log_out = log_out, clock = clock)
}
cpp_p <- function(x, a, b, h, lambda, delta = 0, clock = 0, log_out = FALSE) {
  ch <- to_chart(a, b, h, lambda, delta)
  EMC2:::pfrq(x, ch[["K"]], ch[["v"]], ch[["sv"]], h, delta = delta,
              log_out = log_out, clock = clock)
}
cpp_s <- function(x, a, b, h, lambda, delta = 0, clock = 0, log_out = FALSE) {
  ch <- to_chart(a, b, h, lambda, delta)
  EMC2:::pfrq(x, ch[["K"]], ch[["v"]], ch[["sv"]], h, delta = delta,
              lower_tail = FALSE, log_out = log_out, clock = clock)
}
cpp_rate <- function(a, b, h, lambda, delta = 0, clock = 0) {
  ch <- to_chart(a, b, h, lambda, delta)
  EMC2:::frq_rate(ch[["K"]], ch[["v"]], ch[["sv"]], h, delta, clock)
}

frq_grid <- data.frame(
  a      = c(1,    1,    2,    1.7,  3,     8,    2.5,  1),
  b      = c(1,    3,    1,    2.3,  6,     2.3,  1.4,  6),
  h      = c(0.3,  0.95, 0.75, 0.95, 0.999, 0.6,  0.85, 0.999),
  lambda = c(2.0,  1.5,  1.2,  0.8,  3.0,   1.0,  1.4,  2.5))
frq_x <- c(1e-4, 0.01, 0.05, 0.2, 0.5, 1, 2, 5, 20)

# ---------------------------------------------------------------------------
# 1. The literal finite-reservoir process
# ---------------------------------------------------------------------------

test_that("the closed form is the K-th order statistic of a finite reservoir", {
  skip_on_cran()
  set.seed(20260816)
  nsim <- 100
  N <- 7; K <- 3; p <- 0.8; lam <- 2.5
  a <- K; b <- N - K + 1
  avail <- matrix(runif(nsim * N) < p, nrow = nsim)
  lat <- matrix(rexp(nsim * N, lam), nrow = nsim)
  lat[!avail] <- Inf
  Tk <- apply(lat, 1, function(r) sort(r)[K])

  h <- pbeta(p, a, b)

  expect_true(all(Tk[!is.infinite(Tk)] > 0))
  expect_equal(cpp_p(Inf, a, b, h, lam), h, tolerance = 1e-12)

  # ... and the kernel must recover exactly the (p, lambda) that generated it.
  got <- EMC2:::frq_rate(K, b * p * lam, p * lam * sqrt(b), h, clock = 0)
  expect_equal(unname(got[1, "p"]), p, tolerance = 1e-9)
  expect_equal(unname(got[1, "lambda"]), lam, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 2. Against an independent transcription of the boxed formulas
# ---------------------------------------------------------------------------

test_that("pdf, cdf and survivor match an independent R reference", {
  worst <- c(d = 0, p = 0, s = 0)
  rel <- function(got, want) max(abs(got - want) / pmax(want, 1e-300))
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    worst["d"] <- max(worst["d"], rel(cpp_d(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_d(frq_x, g$a, g$b, g$h, g$lambda)))
    worst["p"] <- max(worst["p"], rel(cpp_p(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_p(frq_x, g$a, g$b, g$h, g$lambda)))
    worst["s"] <- max(worst["s"], rel(cpp_s(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_s(frq_x, g$a, g$b, g$h, g$lambda)))
  }
  expect_lt(worst[["d"]], 1e-11)
  expect_lt(worst[["p"]], 1e-11)
  expect_lt(worst[["s"]], 1e-10)
})

test_that("the log output equals the log of the natural output", {
  for (i in c(1, 4, 5)) {
    g <- frq_grid[i, ]
    expect_equal(cpp_d(frq_x, g$a, g$b, g$h, g$lambda, log_out = TRUE),
                 log(cpp_d(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
    expect_equal(cpp_p(frq_x, g$a, g$b, g$h, g$lambda, log_out = TRUE),
                 log(cpp_p(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
    expect_equal(cpp_s(frq_x, g$a, g$b, g$h, g$lambda, log_out = TRUE),
                 log(cpp_s(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 3. The defective-distribution contract
# ---------------------------------------------------------------------------

test_that("the asymptotic and integral checks of the defective tail hold", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$lambda), g$h, tolerance = 1e-12)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$lambda), 1 - g$h, tolerance = 1e-9)
    I <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$lambda), 0, Inf,
                   rel.tol = 1e-10)$value
    expect_equal(I, g$h, tolerance = 1e-7)
    for (x in c(0.05, 0.3, 1.5)) {
      Ix <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$lambda), 0, x,
                      rel.tol = 1e-10)$value
      expect_equal(Ix, cpp_p(x, g$a, g$b, g$h, g$lambda), tolerance = 1e-7)
    }
    expect_equal(cpp_p(frq_x, g$a, g$b, g$h, g$lambda) +
                   cpp_s(frq_x, g$a, g$b, g$h, g$lambda),
                 rep(1, length(frq_x)), tolerance = 1e-12)
  }
})

test_that("the CDF is monotone and the density non-negative and finite", {
  xs <- exp(seq(log(1e-6), log(50), length.out = 31))
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    d <- cpp_d(xs, g$a, g$b, g$h, g$lambda)
    p <- cpp_p(xs, g$a, g$b, g$h, g$lambda)
    expect_false(any(is.na(d)) || any(is.na(p)))
    expect_true(all(d >= 0) && all(is.finite(d)))
    expect_true(all(diff(p) >= -1e-14))
    expect_true(all(p >= 0 & p <= g$h + 1e-12))
  }
})

# ---------------------------------------------------------------------------
# 4. The (K, v, sv, h) counter coordinates
# ---------------------------------------------------------------------------

test_that("h is the completion probability and v the initial input rate", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$lambda), g$h, tolerance = 1e-12)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$lambda), 1 - g$h, tolerance = 1e-9)
  }
  # At K = 1 the density at the onset is the expected initial event rate, for
  # any rate variability, completion probability or clock shape.
  for (cv in c(0, 0.3, 1)) for (h in c(1, 0.6)) for (cu in c(0, 1, 2))
    expect_equal(EMC2:::dfrq(1e-10, 1, 7.5, cv * 7.5, h, clock = cu), 7.5,
                 tolerance = 1e-8, info = paste(cv, h, cu))
})

test_that("frq_rate is the exact inverse of the counter coordinates", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    got <- cpp_rate(g$a, g$b, g$h, g$lambda)
    want <- ref_pl(g$a, g$b, g$h, g$lambda)
    expect_equal(unname(got[1, "beta"]), g$b, tolerance = 1e-12)
    expect_equal(unname(got[1, "p"]), want[["p"]], tolerance = 1e-12)
    expect_equal(unname(got[1, "lambda"]), g$lambda, tolerance = 1e-10)
    expect_equal(unname(got[1, "Etot"]), g$b * want[["p"]] / (1 - want[["p"]]),
                 tolerance = 1e-10)
    expect_equal(pbeta(unname(got[1, "p"]), g$a, g$b), g$h, tolerance = 1e-10)
  }
  v <- EMC2:::frq_rate(c(1, 2, 3), 4, 0.5, 0.9)
  expect_equal(dim(v), c(3L, 5L))
  expect_equal(colnames(v), c("beta", "p", "lambda", "Etot", "log1mp"))
  expect_equal(unname(v[, "log1mp"]), log1p(-unname(v[, "p"])), tolerance = 1e-12)
  expect_equal(v[2, ], EMC2:::frq_rate(2, 4, 0.5, 0.9)[1, ])
})

test_that("sv = 0 is the exact Poisson-counter limit", {
  x <- c(0.01, 0.1, 0.3, 1, 3)
  K <- 2.7; v <- 6
  for (cu in c(0, 1, 2)) {
    # Proper: a Gamma(K, v) decision time, and the clock drops out.
    expect_equal(EMC2:::pfrq(x, K, v, 0, 1, clock = cu), pgamma(v * x, K),
                 tolerance = 1e-14)
    expect_equal(EMC2:::dfrq(x, K, v, 0, 1, clock = cu), dgamma(x, K, v),
                 tolerance = 1e-12)
    # Defective: finite expected evidence E = P_K^{-1}(h) on the clock v / E.
    E <- qgamma(0.8, K); lam <- v / E
    G <- if (cu == 0) -expm1(-lam * x) else 1 - (1 + cu^2 * lam * x)^(-1 / cu^2)
    expect_equal(EMC2:::pfrq(x, K, v, 0, 0.8, clock = cu), pgamma(E * G, K),
                 tolerance = 1e-14)
    expect_equal(EMC2:::pfrq(Inf, K, v, 0, 0.8, clock = cu), 0.8, tolerance = 1e-14)
    expect_equal(EMC2:::pfrq(Inf, K, v, 0, 0.8, clock = cu, lower_tail = FALSE),
                 0.2, tolerance = 1e-14)
    expect_equal(integrate(function(t) EMC2:::dfrq(t, K, v, 0, 0.8, clock = cu),
                           0, Inf, rel.tol = 1e-10)$value, 0.8, tolerance = 1e-8)
    # The sv > 0 family converges to it.
    prev <- Inf
    for (cv in c(0.1, 0.01, 0.001)) {
      dev <- max(abs(EMC2:::pfrq(x, K, v, cv * v, 0.8, clock = cu) -
                     EMC2:::pfrq(x, K, v, 0, 0.8, clock = cu)))
      expect_lt(dev, prev); prev <- dev
    }
    expect_lt(prev, 1e-5)
  }
  # Below a rate CV of 1e-4 the kernel takes the limit branch itself, and the
  # step there is the O(cv^2) gap, not a loss of accuracy.
  for (h in c(1, 0.8)) {
    expect_identical(EMC2:::pfrq(x, K, v, 0.99e-4 * v, h),
                     EMC2:::pfrq(x, K, v, 0, h))
    expect_lt(max(abs(EMC2:::pfrq(x, K, v, 1.01e-4 * v, h) -
                      EMC2:::pfrq(x, K, v, 0, h))), 1e-8)
  }
  r <- EMC2:::frq_rate(K, v, 0, c(1, 0.8))
  expect_equal(unname(r[, "beta"]), c(Inf, Inf))
  expect_equal(unname(r[, "Etot"]), c(Inf, qgamma(0.8, K)), tolerance = 1e-12)
})

test_that("a constant-rate clock is a counter with a Gamma-mixed input rate", {
  # Lambda ~ Gamma(mean v, SD sv) and Gamma(K) epochs: the decision time is a
  # scaled beta-prime, F(x) = I_{x / (x + b)}(K, a), a = (v/sv)^2, b = a / v.
  x <- c(0.02, 0.1, 0.4, 1.5, 6)
  for (K in c(1, 2.4, 5)) for (cv in c(0.2, 0.6, 1)) {
    a <- 1 / cv^2; b <- a / 9
    expect_equal(EMC2:::pfrq(x, K, 9, cv * 9, 1, clock = 1),
                 pbeta(x / (x + b), K, a), tolerance = 1e-12)
  }
})

test_that("K = N removes p from the conditional RT distribution", {
  a <- 4
  cond <- function(x, b, h, lam = 2) {
    cpp_p(x, a, b, h, lam) / h
  }
  xs <- c(0.1, 0.3, 0.9)
  expect_equal(cond(xs, 1, 0.4), cond(xs, 1, 0.95), tolerance = 1e-10)
  expect_gt(max(abs(cond(xs, 3, 0.4) - cond(xs, 3, 0.95))), 1e-3)
})

test_that("alpha controls the order of the leading edge", {
  xs <- c(1e-7, 2e-7)
  for (a in c(1, 1.5, 3, 5.5)) {
    slope <- diff(log(cpp_d(xs, a, 2.5, 0.9, 1.4))) / diff(log(xs))
    expect_equal(slope, a - 1, tolerance = 1e-4)
  }
})

test_that("degenerate parameters are rejected rather than returning garbage", {
  # Columns are (K, v, sv, h).
  bad <- list(c(0, 2, 0.5, 0.9), c(2, 0, 0.5, 0.9), c(2, -1, 0.5, 0.9),
              c(2, 2, -0.1, 0.9), c(2, 2, 0.5, 0), c(2, 2, 0.5, 1.5),
              c(Inf, 2, 0.5, 0.9), c(NA, 2, 0.5, 0.9), c(2, 2, Inf, 0.9))
  for (v in bad) {
    expect_equal(EMC2:::dfrq(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_equal(EMC2:::pfrq(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_true(is.na(EMC2:::frq_rate(v[1], v[2], v[3], v[4])[1, "p"]))
  }
  expect_equal(cpp_d(c(-1, 0), 2, 3, 0.9, 1.3), c(0, 0))
  expect_equal(cpp_p(c(-1, 0), 2, 3, 0.9, 1.3), c(0, 0))
  expect_equal(cpp_s(c(-1, 0), 2, 3, 0.9, 1.3), c(1, 1))
})

test_that("h at its upper bound stays finite and reaches the non-defective limit", {
  h <- 1 - 1e-9
  expect_equal(cpp_p(Inf, 2, 3, h, 1.3), h, tolerance = 1e-15)
  expect_equal(cpp_s(Inf, 2, 3, h, 1.3), 1 - h, tolerance = 1e-5)
  expect_true(all(is.finite(cpp_d(frq_x, 2, 3, h, 1.3))))
  expect_equal(unname(cpp_rate(2, 3, 1, 1.3)[1, "p"]), 1)
  expect_equal(cpp_p(Inf, 2, 3, 1, 1.3), 1, tolerance = 1e-14)
  expect_equal(cpp_s(Inf, 2, 3, 1, 1.3), 0, tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# 5. The R model interface
# ---------------------------------------------------------------------------

test_that("the constructor exposes the documented contract", {
  m <- FRQ()
  expect_equal(m$type, "RACE")
  expect_equal(m$c_name, "FRQ")
  expect_equal(names(m$p_types),
               c("K", "v", "sv", "h", "t0", "delta", "clock",
                 "pContaminant", "pGuess"))
  expect_equal(m$p_types_canonical, c("K", "v", "sv", "h", "t0"))
  expect_equal(unname(m$transform$func[c("K", "v", "sv", "h", "t0",
                                         "delta", "clock")]),
               c("exp", "exp", "exp", "pnorm", "exp", "exp", "exp"))
  expect_equal(unname(m$transform$lower[["K"]]), 1)
  expect_equal(exp(m$p_types[["sv"]]), 0)
  expect_equal(exp(m$p_types[["delta"]]), 0)
  expect_equal(exp(m$p_types[["clock"]]), 1)
  expect_equal(pnorm(m$p_types[["h"]]), 1, tolerance = 1e-12)
  expect_equal(unname(m$bound$minmax[, "K"]), c(1, Inf))
  expect_equal(unname(m$bound$minmax[, "v"]), c(1e-4, Inf))
  expect_equal(unname(m$bound$minmax[, "sv"]), c(1e-4, Inf))
  expect_equal(unname(m$bound$minmax[, "delta"]), c(1e-4, 6))
  expect_equal(unname(m$bound$minmax[, "clock"]), c(1e-4, Inf))
  expect_equal(unname(m$bound$exception[c("sv", "delta", "clock", "h", "K")]),
               c(0, 0, 0, 1, 1))
  expect_length(formals(FRQ), 0L)
})

test_that("the counter coordinates are representable on the whole bounded box", {
  grid <- data.frame(
    K = c(1, 1.5, 2, 20, 200),
    v = c(0.1, 2, 10, 40, 500),
    cv = c(5, 1, 0.3, 0.05, 1e-3),
    h = c(1e-3, 0.05, 0.5, 0.9, 1 - 1e-9)
  )
  g <- EMC2:::frq_rate(grid$K, grid$v, grid$cv * grid$v, grid$h)
  expect_false(anyNA(g))
  expect_true(all(g[, "p"] > 0 & g[, "p"] < 1))
  expect_true(all(is.finite(g[, "lambda"]) & g[, "lambda"] > 0))
  expect_equal(pbeta(g[, "p"], grid$K, 1 / grid$cv^2), grid$h, tolerance = 1e-8)
  expect_equal(unname(g[, "beta"] * g[, "p"] * g[, "lambda"]), grid$v,
               tolerance = 1e-12)
})

test_that("dfun/pfun apply t0 and call the same kernel", {
  pars <- cbind(K = 2, v = 3.1, sv = 1.86, h = 0.9, t0 = 0.15, delta = 0,
                clock = 0.5)
  pars <- pars[rep(1, 5), ]
  rt <- c(0.05, 0.15, 0.3, 0.8, Inf)
  x <- c(0.15, 0.65)
  expect_equal(EMC2:::dFRQ(rt, pars),
               c(0, 0, EMC2:::dfrq(x, 2, 3.1, 1.86, 0.9, clock = 0.5), 0))
  expect_equal(EMC2:::pFRQ(rt, pars),
               c(0, 0, EMC2:::pfrq(x, 2, 3.1, 1.86, 0.9, clock = 0.5), 0.9))
  expect_equal(EMC2:::sFRQ(rt, pars),
               c(1, 1, EMC2:::pfrq(x, 2, 3.1, 1.86, 0.9, clock = 0.5,
                                   lower_tail = FALSE), 0.1))
  # A design without clock reads the model default, a constant-rate clock.
  expect_equal(EMC2:::dFRQ(rt, pars[, colnames(pars) != "clock"]),
               c(0, 0, EMC2:::dfrq(x, 2, 3.1, 1.86, 0.9, clock = 1), 0))
})

test_that("Ttransform reports the generative parameters", {
  pars <- cbind(K = 2, v = 3.1, sv = 1.86, h = 0.9, t0 = 0.15, delta = 0.3,
                clock = 1)
  out <- FRQ()$Ttransform(pars, NULL)
  b <- 1 / 0.6^2
  p <- qbeta(EMC2:::frq_h_inv_r(0.9, 0.3), 2, b)
  expect_equal(unname(out[1, "beta"]), b, tolerance = 1e-12)
  expect_equal(unname(out[1, "p"]), p, tolerance = 1e-12)
  expect_equal(unname(out[1, "lambda"]), 3.1 / (b * p), tolerance = 1e-10)
  expect_equal(unname(out[1, "Etot"]), b * p / (1 - p), tolerance = 1e-10)
  expect_equal(unname(out[1, "cv"]), 0.6)
  expect_equal(unname(out[1, "sQ"]), 0.3 / sqrt(3))
  expect_equal(out[, colnames(pars), drop = FALSE], pars)
  expect_null(rownames(out))
})

# ---------------------------------------------------------------------------
# 6. The compiled race likelihood
# ---------------------------------------------------------------------------

frq_ll <- function(emc, p, p_types_override = NULL) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  if (!is.null(p_types_override)) p_types <- p_types_override
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform, p_types = p_types,
                    min_ll = log(1e-10), trend = model$trend)
}

ref_race_ll <- function(dadm, pars) {
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  d <- EMC2:::dFRQ(dadm$rt, pars)
  s <- EMC2:::sFRQ(dadm$rt, pars)
  win <- as.logical(dadm$winner)
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    if (!any(win[ix])) return(sum(log(s[ix])))
    log(d[ix][win[ix]]) + sum(log(s[ix][!win[ix]]))
  }, numeric(1))
  sum(pmax(ll, log(1e-10)))
}

frq_truth <- list(hT = pnorm(qnorm(0.93) + 0.45), hF = pnorm(qnorm(0.93) - 0.45))

frq_fixture <- function(n = 60, seed = 20260816, extensions = FALSE) {
  set.seed(seed)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  formula <- list(K ~ 1, v ~ lM, sv ~ 1, h ~ lM, t0 ~ 1)
  if (extensions)
    formula[c("delta", "clock")] <- list(delta ~ 1, clock ~ 1)
  des <- design(data = dat, model = FRQ, matchfun = matchfun,
                formula = formula,
                contrasts = list(h = list(lM = ADmat), v = list(lM = ADmat)))
  p <- sampled_pars(des)
  p["K"] <- log(1.2); p["sv"] <- log(1.65)
  p["h"] <- qnorm(0.93); p["h_lMd"] <- 0.9
  p["v"] <- log(3); p["v_lMd"] <- 0.6
  p["t0"] <- log(0.15)
  if (extensions) {
    p["delta"] <- log(1.2)
    p["clock"] <- log(0.5)
  }
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, p = p, emc = emc)
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- frq_fixture(n = 60)
  got <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll(dadm, pars), tolerance = 1e-6)
})

test_that("the compiled likelihood propagates threshold variability and frailty", {
  skip_on_cran()
  fx <- frq_fixture(n = 40, extensions = TRUE)
  got <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll(dadm, pars), tolerance = 1e-6)
})

test_that("a p_types reordering is caught by the column contract", {
  skip_on_cran()
  fx <- frq_fixture(n = 30)
  swapped <- names(FRQ()$p_types)
  swapped[3:4] <- swapped[4:3]
  expect_error(frq_ll(fx$emc, fx$p, p_types_override = swapped),
               "FRQ kernels expect parameter column")
})

test_that("an omission scores exactly sum_j log(1 - h_j)", {
  skip_on_cran()
  fx <- frq_fixture(n = 30)
  d_om <- fx$dat
  d_om$rt[1:6] <- Inf
  d_om$R[1:6] <- NA
  mk <- function(d) suppressMessages(
    make_emc(d, fx$des, type = "single", n_chains = 1, compress = FALSE,
             rt_resolution = NULL))
  ll_om <- frq_ll(mk(d_om), fx$p)
  expect_true(is.finite(ll_om))
  per_omission <- log(1 - frq_truth$hT) + log(1 - frq_truth$hF)
  expect_equal(ll_om - frq_ll(mk(d_om[-(1:6), ]), fx$p), 6 * per_omission,
               tolerance = 1e-6)
})

test_that("a finite UT normaliser keeps the retained never-finish atom", {
  skip_on_cran()
  fx <- frq_fixture(n = 40)
  LT <- 0.2; UT <- 0.9
  d2 <- fx$dat
  d2$rt[1:6] <- Inf
  d2$R[1:6] <- NA
  d2 <- d2[d2$rt >= LT | is.infinite(d2$rt), ]
  d2 <- d2[d2$rt <= UT | is.infinite(d2$rt), ]
  d2$LT <- LT; d2$UT <- UT
  e2 <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                  compress = FALSE, rt_resolution = NULL))
  dadm <- e2[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, e2[[1]]$model())

  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  ref_trunc_ll <- function(p) {
    pars <- EMC2:::get_pars_matrix_oo(p, dadm, e2[[1]]$model())
    s_at <- function(t) EMC2:::sFRQ(rep(t, nrow(dadm)), pars)
    s_LT <- s_at(LT); s_UT <- s_at(UT); s_Inf <- s_at(Inf)
    d <- EMC2:::dFRQ(dadm$rt, pars)
    s <- EMC2:::sFRQ(dadm$rt, pars)
    win <- as.logical(dadm$winner)
    ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
      logZ <- log(prod(s_LT[ix]) - prod(s_UT[ix]) + prod(s_Inf[ix]))
      num <- if (!is.finite(dadm$rt[ix][1])) sum(log(s[ix]))
             else log(d[ix][win[ix]]) + sum(log(s[ix][!win[ix]]))
      num - logZ
    }, numeric(1))
    sum(pmax(ll, log(1e-10)))
  }
  expect_equal(frq_ll(e2, fx$p), ref_trunc_ll(fx$p), tolerance = 1e-6)

  p_low <- fx$p
  p_low[["h"]] <- qnorm(0.35); p_low[["h_lMd"]] <- 0
  expect_equal(frq_ll(e2, p_low), ref_trunc_ll(p_low), tolerance = 1e-6)

  pars_low <- EMC2:::get_pars_matrix_oo(p_low, dadm, e2[[1]]$model())
  s_low <- function(t) prod(EMC2:::sFRQ(rep(t, nrow(dadm)), pars_low)[seq_len(n_lR)])
  Z_fixed <- s_low(LT) - s_low(UT) + s_low(Inf)
  Z_finite_only <- s_low(LT) - s_low(UT)
  expect_gt(s_low(Inf), 0.8 * Z_finite_only)
  expect_gt(40 * (log(Z_fixed) - log(Z_finite_only)), 10)
})

test_that("truncation, censoring and pContaminant stay well posed", {
  skip_on_cran()
  fx <- frq_fixture(n = 40)
  ll_base <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(ll_base))
  for (bounds in list(list(LT = 0.05, UT = 3), list(LC = 0.05, UC = 2))) {
    d2 <- fx$dat
    for (nm in names(bounds)) d2[[nm]] <- bounds[[nm]]
    e2 <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                    compress = FALSE, rt_resolution = NULL))
    expect_true(is.finite(frq_ll(e2, fx$p)))
  }

  matchfun <- function(d) d$S == d$lR
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des_pc <- suppressMessages(design(
    data = fx$dat, model = FRQ, matchfun = matchfun,
    formula = list(K ~ 1, v ~ lM, sv ~ 1, h ~ lM, t0 ~ 1,
                   pContaminant ~ 1),
    contrasts = list(h = list(lM = ADmat), v = list(lM = ADmat))))
  p_pc <- c(fx$p, pContaminant = qnorm(0.05))[names(sampled_pars(des_pc))]
  e_pc <- suppressMessages(make_emc(fx$dat, des_pc, type = "single",
                                    n_chains = 1, compress = FALSE,
                                    rt_resolution = NULL))
  # The contaminant's link carries its bound, so read its natural value back.
  pc <- EMC2:::get_pars_matrix_oo(p_pc, e_pc[[1]]$data[[1]],
                                  e_pc[[1]]$model())[[1, "pContaminant"]]
  expect_equal(frq_ll(e_pc, p_pc), ll_base + nrow(fx$dat) * log(1 - pc),
               tolerance = 1e-6)
})

test_that("a response below t0 floors cleanly", {
  skip_on_cran()
  fx <- frq_fixture(n = 20)
  p <- fx$p
  p[["t0"]] <- log(max(fx$dat$rt) + 1)
  ll <- frq_ll(fx$emc, p)
  expect_false(is.na(ll))
  expect_equal(ll, nrow(fx$dat) * log(1e-10), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 7. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree distributionally with the CDF", {
  skip_on_cran()
  local_test_rng(kind = "Mersenne-Twister")
  # n = 4000 puts the 0.03 CDF and 0.006 omission-rate tolerances at ~3.8
  # standard errors; at n = 1000 they were ~1.9 and failed for many seeds.
  n <- 4000
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  for (cv in c(0.6, 0)) {
  pars <- cbind(K = 2, v = 3, sv = 3 * cv, h = 0.9, t0 = 0.1)
  pars <- pars[rep(1, length(lR)), ]
  set.seed(4)
  a <- EMC2:::rFRQ(lR, pars)
  set.seed(4)
  b <- withr::with_options(list(emc2.cpp_rfun = TRUE),
                           EMC2:::.rfun_FRQ(lR, pars))
  probe <- c(0.2, 0.35, 0.5, 0.8)
  th <- 1 - EMC2:::pfrq(probe - 0.1, 2, 3, 3 * cv, 0.9, lower_tail = FALSE)^2
  for (dat in list(a, b)) {
    fin <- is.finite(dat$rt)
    emp <- vapply(probe, function(x) mean(fin & dat$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.03)
    expect_lt(abs(mean(!fin) - 0.01), 0.006)
    expect_true(all(is.na(dat$R[!fin])))
    expect_true(all(dat$rt[fin] > 0.1))
  }
  }
  # The proper Poisson counter (sv = 0, h = 1) is a Gamma(K, v) race.
  pars <- cbind(K = 2.5, v = 6, sv = 0, h = 1, t0 = 0)[rep(1, 2 * n), ]
  for (sim in list(EMC2:::rFRQ(lR, pars), EMC2:::rfrq_cpp(pars, levels(lR), rep(TRUE, 2 * n)))) {
    expect_true(all(is.finite(sim$rt)))
    th <- 1 - pgamma(6 * probe, 2.5, lower.tail = FALSE)^2
    emp <- vapply(probe, function(x) mean(sim$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.03)
  }
})

test_that("the simulator reproduces the literal reservoir at integer shapes", {
  skip_on_cran()
  set.seed(99)
  N <- 6; K <- 2; p <- 0.7; lam <- 3
  a <- K; b <- N - K + 1
  h <- pbeta(p, a, b)
  n <- 1e4
  lR <- factor(rep("go", n), levels = "go")
  pars <- cbind(K = a, v = b * p * lam, sv = p * lam * sqrt(b), h = h, t0 = 0,
                clock = 0)
  sim <- EMC2:::.rfun_FRQ(lR, pars[rep(1, n), ])

  avail <- matrix(runif(n * N) < p, nrow = n)
  lat <- matrix(rexp(n * N, lam), nrow = n)
  lat[!avail] <- Inf
  ref <- apply(lat, 1, function(r) sort(r)[K])

  probe <- c(0.05, 0.15, 0.3, 0.6, 1.2)
  emp <- vapply(probe, function(x) mean(is.finite(sim$rt) & sim$rt <= x),
                numeric(1))
  ref_e <- vapply(probe, function(x) mean(ref <= x), numeric(1))
  expect_lt(max(abs(emp - ref_e)), 0.02)
  expect_lt(abs(mean(!is.finite(sim$rt)) - mean(is.infinite(ref))), 0.015)
})

test_that("make_data produces omissions the design can be fit back through", {
  skip_on_cran()
  set.seed(20260816)
  fx <- frq_fixture(n = 30)
  p_omit <- fx$p
  p_omit[["h"]] <- qnorm(0.8)
  p_omit[["h_lMd"]] <- 0.4
  sim <- make_data(p_omit, design = fx$des, n_trials = 150)
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] > 0.15))
  expect_gt(mean(sim$S[fin] == sim$R[fin]), 0.6)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(frq_ll(e, p_omit)))
})

# ---------------------------------------------------------------------------
# 8. Recovery
# ---------------------------------------------------------------------------

test_that("the likelihood is maximised at the generating parameters", {
  skip_on_cran()
  set.seed(20260816)
  fx <- frq_fixture(n = 30)
  p_rec <- fx$p
  p_rec["h"] <- qnorm(0.7)
  p_rec["h_lMd"] <- 0.3
  sim <- make_data(p_rec, design = fx$des, n_trials = 3000)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  ll0 <- frq_ll(e, p_rec)
  expect_true(is.finite(ll0))
  for (nm in names(p_rec)) {
    for (delta in c(-0.5, 0.5)) {
      p2 <- p_rec
      p2[[nm]] <- p2[[nm]] + delta
      expect_lt(frq_ll(e, p2), ll0)
    }
  }
})

# ---------------------------------------------------------------------------
# 9. Threshold variability (delta) and unit-rate frailty (clock)
# ---------------------------------------------------------------------------

test_that("delta defaults to zero and is exactly the untransformed model", {
  x <- c(0.05, 0.2, 0.6, 2, Inf)
  expect_identical(EMC2:::pfrq(x, 2, 3, 0.6, 0.9, 0),
                   EMC2:::pfrq(x, 2, 3, 0.6, 0.9))
  expect_identical(EMC2:::dfrq(x, 2, 3, 0.6, 0.9, 0),
                   EMC2:::dfrq(x, 2, 3, 0.6, 0.9))
  expect_identical(EMC2:::frq_rate(2, 3, 0.6, 0.9, 0),
                   EMC2:::frq_rate(2, 3, 0.6, 0.9))
  expect_equal(cpp_p(x[1:4], 2, 3, 0.9, 1.4, 0), ref_p(x[1:4], 2, 3, 0.9, 1.4))
})

test_that("unit-rate frailty matches the closed-form registration mixture", {
  a <- 2.5; b <- 1.4; h <- 0.85; lambda <- 1.4; cv <- 0.7
  pl <- cpp_rate(a, b, h, lambda, 0, cv)
  p <- qbeta(h, a, b)
  expect_equal(unname(pl[1, "p"]), p, tolerance = 1e-12)
  expect_equal(unname(pl[1, "lambda"]), lambda, tolerance = 1e-12)

  x <- c(0.05, 0.2, 0.6, 1.5, 5)
  c2 <- cv^2
  sr <- (1 + c2 * lambda * x)^(-1 / c2)
  g <- lambda * (1 + c2 * lambda * x)^(-1 / c2 - 1)
  q <- p * (1 - sr)
  want_p <- pbeta(q, a, b)
  want_d <- p * g * q^(a - 1) * (1 - q)^(b - 1) / beta(a, b)
  expect_equal(cpp_p(x, a, b, h, lambda, 0, clock = cv), want_p,
               tolerance = 1e-11)
  expect_equal(cpp_d(x, a, b, h, lambda, 0, clock = cv), want_d,
               tolerance = 1e-11)
  expect_equal(cpp_p(Inf, a, b, h, lambda, 0, clock = cv), h,
               tolerance = 1e-11)
  expect_equal(cpp_s(Inf, a, b, h, lambda, 0, clock = cv), 1 - h,
               tolerance = 1e-11)
  expect_equal(integrate(function(z)
    cpp_d(z, a, b, h, lambda, 0, clock = cv), 0, Inf,
    rel.tol = 1e-9)$value, h, tolerance = 1e-8)
})

test_that("unit-rate frailty nests the exponential member exactly", {
  x <- c(0.05, 0.2, 0.6, 2, Inf)
  # The exponential clock is the clock = 0 member of the Lomax family.
  prev <- Inf
  for (cu in c(1e-2, 1e-3, 1e-4)) {
    dev <- max(abs(cpp_p(x, 2, 3, 0.9, 1.4, 0, clock = cu) -
                   cpp_p(x, 2, 3, 0.9, 1.4, 0, clock = 0)))
    expect_lt(dev, prev); prev <- dev
  }
  expect_lt(prev, 1e-4)
  expect_equal(cpp_p(x[1:4], 2, 3, 0.9, 1.4, 0, clock = 0),
               ref_p(x[1:4], 2, 3, 0.9, 1.4), tolerance = 1e-12)
})

test_that("threshold variability and unit frailty compose", {
  a <- 2.5; b <- 1.4; h <- 0.85; lambda <- 1.4
  d <- 1.25; cv <- 0.55
  x <- c(0.05, 0.2, 0.6, 1.5, 5)
  pl <- cpp_rate(a, b, h, lambda, d, cv)
  p <- unname(pl[1, "p"])
  c2 <- cv^2
  q <- p * (1 - (1 + c2 * lambda * x)^(-1 / c2))
  z <- pbeta(q, a, b)
  want <- vapply(z, ref_H, numeric(1), d = d)
  expect_equal(cpp_p(x, a, b, h, lambda, d, clock = cv), want,
               tolerance = 1e-10)
  expect_equal(cpp_p(Inf, a, b, h, lambda, d, clock = cv), h,
               tolerance = 1e-10)
  expect_equal(integrate(function(t)
    cpp_d(t, a, b, h, lambda, d, clock = cv), 0, Inf,
    rel.tol = 1e-9)$value, h, tolerance = 1e-8)
})

test_that("the proper h = 1 boundary remains proper with unit frailty", {
  pars <- cbind(K = 2, v = 4.2, sv = 2.436, h = 1, t0 = 0, delta = 0,
                clock = 0.7)
  pl <- EMC2:::frq_rate(2, 4.2, 2.436, 1, 0, 0.7)
  expect_equal(unname(pl[1, "p"]), 1)
  expect_equal(unname(pl[1, "Etot"]), Inf)
  expect_equal(EMC2:::pfrq(Inf, 2, 4.2, 2.436, 1, 0, clock = 0.7), 1)
  expect_equal(EMC2:::pfrq(Inf, 2, 4.2, 2.436, 1, 0, clock = 0.7,
                            lower_tail = FALSE), 0)
  set.seed(20260920)
  sim_cpp <- EMC2:::rfrq_cpp(pars[rep(1, 2000), ], "1", rep(TRUE, 2000))
  expect_true(all(is.finite(sim_cpp$rt)))
  sim_r <- EMC2:::rFRQ(factor(rep("1", 2000)),
                       pars[rep(1, 2000), ])
  expect_true(all(is.finite(sim_r$rt)))
  out <- FRQ()$Ttransform(pars, NULL)
  expect_equal(unname(out[1, "lambda"]), 4.2 / (1 / 0.58^2), tolerance = 1e-12)
})

test_that("the delta > 0 CDF matches the marginalised criterion state", {
  x <- c(0.02, 0.1, 0.35, 0.9, 4)
  for (d in c(0.25, 1, 3, 6)) {
    for (i in c(1, 2, 5)) {
      g <- frq_grid[i, ]
      expect_equal(cpp_p(x, g$a, g$b, g$h, g$lambda, d),
                   ref_pH(x, g$a, g$b, g$h, g$lambda, d),
                   tolerance = 1e-9,
                   info = paste("delta", d, "row", i))
    }
  }
})

test_that("delta preserves the defective contract and the meaning of h and lambda", {
  for (d in c(0.25, 1, 3, 6)) for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    lab <- paste("delta", d, "row", i)
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$lambda, d), g$h,
                 tolerance = 1e-10, info = lab)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$lambda, d),
                 1 - g$h, tolerance = 1e-10, info = lab)
    expect_equal(integrate(function(t) cpp_d(t, g$a, g$b, g$h, g$lambda, d),
                           0, Inf, rel.tol = 1e-10,
                           subdivisions = 2000L)$value,
                 g$h, tolerance = 1e-8, info = lab)
  }
})

test_that("the density is the derivative of the CDF under delta", {
  eps <- 1e-5
  for (d in c(0.5, 2, 5)) {
    x <- c(0.1, 0.3, 0.7, 1.5)
    fd <- (cpp_p(x + eps, 2.5, 1.4, 0.85, 1.4, d) -
           cpp_p(x - eps, 2.5, 1.4, 0.85, 1.4, d)) / (2 * eps)
    expect_equal(cpp_d(x, 2.5, 1.4, 0.85, 1.4, d), fd,
                 tolerance = 1e-6, info = paste("delta", d))
  }
})

test_that("H is invertible and the family nests the base model as delta -> 0", {
  y <- c(1e-9, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9)
  for (d in c(0.3, 1, 3, 6))
    expect_equal(vapply(EMC2:::frq_h_inv_r(y, d), ref_H, numeric(1), d = d), y,
                 tolerance = 1e-12, info = paste("delta", d))
  x <- c(0.05, 0.3, 1)
  base <- cpp_p(x, 2, 3, 0.9, 1.4, 0)
  prev <- Inf
  for (d in c(1e-3, 1e-4, 1e-5)) {
    dev <- max(abs(cpp_p(x, 2, 3, 0.9, 1.4, d) - base))
    expect_lt(dev, prev)
    prev <- dev
  }
  expect_lt(prev, 1e-10)
})

test_that("the simulator draws from the delta > 0 likelihood", {
  skip_on_cran()
  set.seed(4)
  n <- 4e4
  d <- 2
  pars <- cbind(K = 2, v = 3, sv = 1.5, h = 0.8, t0 = 0.15,
                delta = d)[rep(1, n), ]
  sim <- EMC2:::rfrq_cpp(pars, "1", rep(TRUE, n))
  fin <- is.finite(sim$rt)
  expect_equal(mean(fin), 0.8, tolerance = 0.01)
  probs <- seq(0.1, 0.9, 0.1)
  qs <- quantile(sim$rt[fin], probs)
  expect_equal(unname(EMC2:::pfrq(qs - 0.15, 2, 3, 1.5, 0.8, d) / 0.8),
               probs, tolerance = 0.02)
  expect_false(isTRUE(all.equal(
    unname(EMC2:::pfrq(qs - 0.15, 2, 3, 1.5, 0.8, 0) / 0.8), probs,
    tolerance = 0.02)))
})

test_that("FRQ keeps the defect when sv > v and rejects extreme delta", {
  # beta = 0.05 puts 1 - p far below double-precision spacing near one; the
  # complementary quantile keeps it, and the defect is exact.
  g <- EMC2:::frq_rate(2, 1.3, 1.3 / sqrt(0.05), 0.99)
  expect_false(anyNA(g))
  expect_equal(unname(g[1, "beta"]), 0.05, tolerance = 1e-12)
  expect_equal(EMC2:::pfrq(Inf, 2, 1.3, 1.3 / sqrt(0.05), 0.99,
                           lower_tail = FALSE), 0.01, tolerance = 1e-12)
  sv <- 1.3 / sqrt(0.05); x <- c(0.1, 1, 10, 1e3)
  omp <- qbeta(0.01, 0.05, 2); lam <- 1.3 / (0.05 * (1 - omp))
  expect_equal(EMC2:::pfrq(x, 2, 1.3, sv, 0.99, lower_tail = FALSE),
               pbeta(omp + (1 - omp) / (1 + lam * x), 0.05, 2), tolerance = 1e-10)
  expect_equal(EMC2:::pfrq(x, 2, 1.3, sv, 0.99, lower_tail = FALSE, clock = 0),
               pbeta(omp + (1 - omp) * exp(-lam * x), 0.05, 2), tolerance = 1e-10)
  expect_true(is.na(EMC2:::frq_rate(2, 1.3, 0.5, 0.9, 1000)[1, "p"]))
  expect_equal(EMC2:::pfrq(0.4, 2, 1.3, 0.5, 0.9, 1000), 0)
})

test_that("the FRQ R simulator keeps malformed rows as omissions", {
  lR <- factor(rep(c("a", "b"), 2), levels = c("a", "b"))
  pars <- cbind(K = 2, v = 3, sv = 1.5, h = 0.9,
                t0 = c(NA, NA, 0.1, 0.1))
  sim <- EMC2:::rFRQ(lR, pars)
  expect_true(is.infinite(sim$rt[1]) || is.na(sim$R[1]))
  expect_error(EMC2:::rFRQ(lR, pars, ok = TRUE), "one value per")
})

test_that("the FRQ R simulator permits a valid accumulator to win", {
  lR <- factor(c("a", "b"), levels = c("a", "b"))
  pars <- cbind(K = c(2, 2), v = c(3, 3), sv = c(1.5, 1.5), h = c(1, 1),
                t0 = c(NA, 0.1))
  sim <- EMC2:::rFRQ(lR, pars)
  expect_equal(sim$R, factor("b", levels = levels(lR)))
  expect_true(is.finite(sim$rt))
})

test_that("sv far above v stays exact where 1 - q underflows", {
  # K = 1 on the exponential clock at h = 1 is exactly Exp(v) for every sv:
  # the first of the units to register.  sv / v = 50 puts 1 - q below double
  # range by u = 0.5.
  u <- c(0.1, 0.5, 1.5, 4)
  for (sv in c(3, 30, 150)) {
    P <- cbind(K = 1, v = 3, sv = sv, h = 1, t0 = 0, delta = 0,
               clock = 0)[rep(1, length(u)), ]
    expect_equal(EMC2:::pFRQ(u, P), pexp(u, 3), tolerance = 1e-12)
    expect_equal(EMC2:::sFRQ(u, P), pexp(u, 3, lower.tail = FALSE),
                 tolerance = 1e-12)
    expect_equal(EMC2:::dFRQ(u, P), dexp(u, 3), tolerance = 1e-10)
  }
  # The simulator holds W = 1 - U as a log for the same reason.
  withr::local_options(emc2.cpp_rfun = TRUE)
  set.seed(5)
  n <- 20000
  P <- cbind(K = 1, v = 3, sv = 150, h = 1, t0 = 0, delta = 0,
             clock = 0)[rep(1, n), ]
  sim <- EMC2:::.rfun_FRQ(factor(rep("a", n)), P)
  expect_gt(ks.test(sim$rt, "pexp", 3)$p.value, 1e-3)
})

test_that("a per-row K matches a fresh derive in every state", {
  # Rows mix the proper sv = 0 state (where only alpha and lgamma(K) change
  # with K) with sv > 0, h < 1 and delta > 0 rows and an invalid K.
  set.seed(9)
  n <- 60
  P <- cbind(K = 1 + rexp(n, 0.5), v = 4, sv = 0, h = 1, t0 = 0.1,
             delta = 0, clock = 1)
  P[c(7, 8, 30), "sv"] <- 1.5
  P[c(12, 31), "h"] <- 0.8
  P[c(20, 21), "delta"] <- 0.7
  P[c(25, 40), "K"] <- 0.5 * 1e-8
  rt <- 0.1 + rexp(n, 2)
  one <- function(f) vapply(seq_len(n), function(i)
    f(rt[i], P[i, , drop = FALSE]), numeric(1))
  expect_equal(EMC2:::dFRQ(rt, P), one(EMC2:::dFRQ), tolerance = 1e-14)
  expect_equal(EMC2:::pFRQ(rt, P), one(EMC2:::pFRQ), tolerance = 1e-14)
  expect_equal(EMC2:::sFRQ(rt, P), one(EMC2:::sFRQ), tolerance = 1e-14)
  expect_equal(EMC2:::pFRQ(rt[1:3], P[1:3, ]),
               pgamma(4 * (rt[1:3] - 0.1), P[1:3, "K"]), tolerance = 1e-12)
})
