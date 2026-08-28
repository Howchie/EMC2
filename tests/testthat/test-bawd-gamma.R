skip_model_validation()

# Regression coverage for fixed clearance fading in BAwD.
# The oracles in this file work in (mu, sigma).  The compiled lognormal launch
# (launch = 1L) is sampled as the natural-scale mean and CV, so the kernel
# wrappers below convert at the boundary.
ln_meancv <- function(p1, p2, launch) {
  if (as.integer(launch) != 1L) return(c(unname(p1), unname(p2)))
  c(unname(exp(p1 + p2^2 / 2)), unname(sqrt(expm1(p2^2))))
}
# All gamma_ref_* helpers below are independent trajectory/quadrature oracles.

gamma_ref_q <- function(u, k) if (k < 1e-12) u else -expm1(-k * u) / k
gamma_ref_c <- function(u, k, gamma) {
  if (k < 1e-12 || gamma < 1e-12) u
  else if (gamma > 1 - 1e-12) gamma_ref_q(u, k)
  else -expm1(-gamma * k * u) / (gamma * k)
}
gamma_ref_psi <- function(s, gamma) {
  if (gamma < 1e-12) expm1(s) - s
  else if (gamma > 1 - 1e-12) 0
  else expm1((1 - gamma) * s) + (1 - gamma) / gamma * expm1(-gamma * s)
}
gamma_ref_sroot <- function(c, gamma) {
  if (!(c > 0) || gamma > 1 - 1e-12) return(0)
  hi <- max(1, log1p(c) + 3)
  while (gamma_ref_psi(hi, gamma) < c) hi <- hi * 2
  uniroot(function(s) gamma_ref_psi(s, gamma) - c, c(1e-12, hi), tol = 1e-13)$root
}
gamma_ref_fp <- function(V, z, b, k, ell, gamma) {
  d <- b - z
  if (d <= 0) return(0)
  if (k < 1e-12) return(if (V > ell) d / (V - ell) else Inf)
  if (gamma > 1 - 1e-12)
    return(if (V > ell && 1 - k * d / (V - ell) > 0)
      -log1p(-k * d / (V - ell)) / k else Inf)
  if (V <= ell) return(Inf)
  up <- log(V / ell) / ((1 - gamma) * k)
  if (V * gamma_ref_q(up, k) - ell * gamma_ref_c(up, k, gamma) < d) return(Inf)
  uniroot(function(u) V * gamma_ref_q(u, k) - ell * gamma_ref_c(u, k, gamma) - d,
          c(0, up), tol = 1e-12)$root
}

gamma_ref_vstar <- function(u, z, b, k, ell, gamma) {
  d <- b - z
  if (d <= 0) return(0)
  if (u <= 0) return(Inf)
  fun <- function(x) (d + ell * gamma_ref_c(x, k, gamma)) / gamma_ref_q(x, k)
  # The minimum over [0,u] is the launch threshold for a first passage by u.
  lo <- 1e-9
  opt <- optimize(fun, c(lo, u), tol = 1e-11)
  min(fun(u), opt$objective)
}
gamma_ref_surv <- function(w, launch, p1, p2, posdrift = TRUE) {
  if (!is.finite(w)) return(0)
  if (launch == 0L) {
    den <- if (posdrift) pnorm(p1 / p2) else 1
    pnorm((p1 - w) / p2) / den
  } else pnorm((p1 - log(w)) / p2)
}
gamma_ref_density <- function(w, launch, p1, p2, posdrift = TRUE) {
  if (launch == 0L) {
    den <- if (posdrift) pnorm(p1 / p2) else 1
    dnorm((w - p1) / p2) / (p2 * den)
  } else dlnorm(w, p1, p2)
}
gamma_ref_F <- function(u, p1, p2, b, A, k, ell, launch, gamma,
                        posdrift = TRUE) {
  if (is.infinite(u)) {
    critical <- function(z) {
      d <- b - z
      if (gamma > 1 - 1e-12) return(ell + k * d)
      ell * exp((1 - gamma) * gamma_ref_sroot(k * d / ell, gamma))
    }
    if (A <= 0) return(gamma_ref_surv(critical(0), launch, p1, p2, posdrift))
    return(integrate(function(z) vapply(z, function(zz)
      gamma_ref_surv(critical(zz), launch, p1, p2, posdrift), numeric(1)),
      0, A, rel.tol = 2e-9, subdivisions = 120L)$value / A)
  }
  if (A <= 0) return(gamma_ref_surv(gamma_ref_vstar(u, 0, b, k, ell, gamma),
                                    launch, p1, p2, posdrift))
  integrate(function(z) vapply(z, function(zz) gamma_ref_surv(
    gamma_ref_vstar(u, zz, b, k, ell, gamma), launch, p1, p2, posdrift),
    numeric(1)), 0, A, rel.tol = 2e-9, subdivisions = 120L)$value / A
}
gamma_ref_f <- function(u, p1, p2, b, A, k, ell, launch, gamma,
                        posdrift = TRUE) {
  q <- gamma_ref_q(u, k)
  if (!(q > 0)) return(0)
  eg <- exp(-gamma * k * u)
  er <- exp(-(1 - gamma) * k * u)
  one <- function(z) {
    d <- b - z
    if (d <= 0) return(0)
    w <- (d + ell * gamma_ref_c(u, k, gamma)) / q
    if (!(w > 0)) return(0)
    zstar <- if (gamma > 1 - 1e-12 || k < 1e-12 || ell <= 1e-12) b else
      b - ell / k * gamma_ref_psi(k * u, gamma)
    if (!(z < zstar)) return(0)
    gamma_ref_density(w, launch, p1, p2, posdrift) *
      eg * (w * er - ell) / q
  }
  if (A <= 0) return(one(0))
  integrate(function(z) vapply(z, one, numeric(1)), 0, A,
            rel.tol = 2e-8, subdivisions = 120L)$value / A
}
gamma_cpp_p <- function(u, p1, p2, b, A, k, ell, launch, gamma) {
  q <- ln_meancv(p1, p2, launch)
  EMC2:::pbawd(t = u, A = A, b = b, p1 = q[1], p2 = q[2], k = k, ell = ell,
               launch = as.integer(launch), posdrift = TRUE, gamma = gamma)
}
gamma_cpp_d <- function(u, p1, p2, b, A, k, ell, launch, gamma) {
  q <- ln_meancv(p1, p2, launch)
  EMC2:::dbawd(t = u, A = A, b = b, p1 = q[1], p2 = q[2], k = k, ell = ell,
               launch = as.integer(launch), posdrift = TRUE, gamma = gamma)
}

test_that("gamma half wall has an independent quadratic oracle", {
  pars <- expand.grid(b = c(.4, 1.1, 2), k = c(.5, .8, 1.7), ell = c(.7, 1, 1.8))
  got <- mapply(function(b, k, ell) EMC2:::bawd_tmax(.3, b, k, ell, .5),
                pars$b, pars$k, pars$ell)
  c0 <- pars$k * pars$b / pars$ell
  r <- (2 + c0 + sqrt(c0^2 + 4 * c0)) / 2
  expect_equal(got, 2 * log(r) / pars$k, tolerance = 2e-9)
  expect_equal(EMC2:::bawd_tmax(.3, 1.1, .8, 1, .5), 2.26676, tolerance = 2e-5)
})

test_that("fading increases omission mass at a common wall", {
  solve_b <- function(g) uniroot(function(b) EMC2:::bawd_tmax(0, b, 1, 1, g) - b,
                                  c(.05, 10))$root
  bs <- vapply(c(0, .5, 2 / 3, .75), solve_b, numeric(1))
  vc <- vapply(seq_along(bs), function(i) {
    q <- ln_meancv(0, 1, 1L)   # the anchors are quoted for mu = 0, sigma = 1
    f <- EMC2:::pbawd(Inf, A = 0, b = bs[i], p1 = q[1], p2 = q[2], k = 1, ell = 1,
                       launch = 1L, posdrift = TRUE, gamma = c(0, .5, 2 / 3, .75)[i])
    exp(-qnorm(f))
  }, numeric(1))
  expect_equal(vc, c(3.513, 5.031, vc[3], vc[4]), tolerance = 1e-2) # first two are anchors
  expect_true(all(diff(vc) > 0))
})

test_that("independent generalized references cover both launches and all regimes", {
  gs <- c(0, .5, 2 / 3, .75, 1)
  for (launch in 0:1) for (g in gs) {
    b <- 1.1; A <- .3; k <- .8; ell <- 1; p1 <- if (launch) 1.2 else .1
    p2 <- if (launch) .9 else .55
    tm <- EMC2:::bawd_tmax(A, b, k, ell, g)
    tt <- c(.2, if (is.finite(tm)) .95 * tm else .8, if (is.finite(tm)) tm else 3, Inf)
    for (u in tt) {
      expect_equal(gamma_cpp_p(u, p1, p2, b, A, k, ell, launch, g),
                   gamma_ref_F(u, p1, p2, b, A, k, ell, launch, g), tolerance = 2e-6)
      if (is.finite(u) && u > .05)
        expect_equal(gamma_cpp_d(u, p1, p2, b, A, k, ell, launch, g),
                     gamma_ref_f(u, p1, p2, b, A, k, ell, launch, g), tolerance = 3e-5)
    }
  }
})

test_that("density integrates to CDF through the fading seams", {
  for (launch in 0:1) for (g in c(0, .5, 2 / 3, .75, 1)) {
    b <- 1.1; A <- .3; k <- .8; ell <- 1; p1 <- if (launch) 1.2 else .1
    p2 <- if (launch) .9 else .55
    for (u in c(.4, 1.1, 2.2)) {
      z <- integrate(function(x) gamma_cpp_d(x, p1, p2, b, A, k, ell, launch, g),
                     1e-5, u, rel.tol = 2e-6)$value
      expect_equal(z, gamma_cpp_p(u, p1, p2, b, A, k, ell, launch, g), tolerance = 1e-4)
    }
  }
})

test_that("frozen fading cases agree with the independent reference", {
  for (g in c(.5, 2 / 3, .75)) for (A in c(1e-6, .2)) {
    args <- list(p1 = -1, p2 = .7, b = 1.1, A = A, k = .8, ell = 1, launch = 1L, gamma = g)
    u <- EMC2:::bawd_tmax(A, args$b, args$k, args$ell, g)
    expect_equal(do.call(gamma_cpp_p, c(list(u = .8 * u), args)),
                 do.call(gamma_ref_F, c(list(u = .8 * u), args)), tolerance = 3e-6)
  }
  args <- list(p1 = .55, p2 = .05, b = 1.1, A = 1e-5, k = .8,
               ell = 1, launch = 1L, gamma = .75)
  tmax <- EMC2:::bawd_tmax(args$A, args$b, args$k, args$ell, args$gamma)
  tsat <- gamma_ref_sroot(args$k * (args$b - args$A) / args$ell,
                          args$gamma) / args$k
  u <- (tsat + tmax) / 2
  got <- do.call(gamma_cpp_p, c(list(u = u), args))
  ref <- do.call(gamma_ref_F, c(list(u = u), args))
  expect_gt(ref, 1e-12)
  expect_equal(got, ref, tolerance = 1e-4)
})

test_that("gamma one has co-decay structure and positive density", {
  expect_true(is.infinite(EMC2:::bawd_tmax(.3, 1.1, .8, 1, 1)))
  for (u in c(.2, 1, 3, 10)) {
    q <- gamma_ref_q(u, .8)
    want <- pnorm((.9 - log(1.1 / q + 1)) / .6)
    qg <- ln_meancv(.9, .6, 1L)
    expect_equal(EMC2:::pbawd(t = u, A = 0, b = 1.1, p1 = qg[1], p2 = qg[2],
                              k = .8, ell = 1, launch = 1L, gamma = 1),
                 want, tolerance = 1e-12)
    expect_gt(EMC2:::dbawd(t = u, A = 0, b = 1.1, p1 = qg[1], p2 = qg[2],
                           k = .8, ell = 1, launch = 1L, gamma = 1), 0)
  }
})

test_that("k-zero and deep co-decay density limits stay on their exact paths", {
  common <- list(t = 1e9, A = .3, b = 1.1, p1 = .9, p2 = .6,
                 k = 1e-10, ell = 1, launch = 1L, posdrift = TRUE)
  p0 <- do.call(EMC2:::pbawd, c(common, list(gamma = 0)))
  d0 <- do.call(EMC2:::dbawd, c(common, list(gamma = 0, log_out = TRUE)))
  for (g in c(.5, 2 / 3, .75, 1)) {
    expect_identical(do.call(EMC2:::pbawd, c(common, list(gamma = g))), p0)
    expect_identical(do.call(EMC2:::dbawd,
                             c(common, list(gamma = g, log_out = TRUE))), d0)
  }

  deep <- function(u, log_out) EMC2:::dbawd(
    t = u, A = .3, b = 1.1, p1 = .9, p2 = .6, k = 1, ell = 1,
    launch = 1L, log_out = log_out, gamma = 1)
  expect_identical(deep(800, FALSE), 0)
  expect_true(is.finite(deep(800, TRUE)))
  expect_equal(deep(801, TRUE) - deep(800, TRUE), -1, tolerance = 1e-12)
})

test_that("lognormal negative-moment survivor integral has quadrature and tail oracles", {
  for (m in c(1, 3)) for (v in c(.4, 1, 3)) {
    got <- exp(EMC2:::lognormal_power_stoploss_log(v, .2, .8, m))
    want <- integrate(function(w) w^(-(m + 1)) * pnorm((.2 - log(w)) / .8),
                      v, Inf, rel.tol = 1e-10)$value
    expect_equal(got, want, tolerance = 2e-8)
  }
  x <- 1.7; v <- exp(.2 + .8 * x)
  expect_equal(exp(EMC2:::lognormal_power_stoploss_log(v, .2, .8, 0)),
               .8 * (dnorm(x) - x * pnorm(-x)), tolerance = 2e-10)
  m <- 3; sg <- 1e-13; vv <- exp(1); a <- m * sg
  xx <- log(vv) / sg
  yy <- xx + a
  series <- a / (xx * yy) -
    a * (xx^2 + xx * yy + yy^2) / (xx^3 * yy^3)
  log_asym <- -m * log(vv) + dnorm(xx, log = TRUE) + log(series) - log(m)
  log_got <- EMC2:::lognormal_power_stoploss_log(vv, 0, sg, m)
  expect_true(is.finite(log_got))
  expect_equal(log_got, log_asym, tolerance = 1e-10)
})
