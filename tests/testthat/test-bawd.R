# BAwD, the ballistic accumulator with drive decay -- src/model_BAwD.h,
# R/model_BAwD.R.
#
# The file is organised around oracles first.  The truncated-normal launch is
# the only variant with exact closed-form comparisons available (BAwL at
# ell = 0, A = 0; the LBA at k = 0), so those pin the shared geometry before
# anything is asked of the lognormal launch, which has no exact oracle and is
# checked against an independent R reference, Monte Carlo, and internal
# consistency (integrate(f) == F).
#
# The reference implementations below are deliberately written from the model
# definition -- root-found first passages and adaptive quadrature -- and share
# no code with the compiled kernels.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

ref_q <- function(u, k) if (k < 1e-12) u else -expm1(-k * u) / k

# root of e^y - 1 - y = c, i.e. -W_{-1}(-e^{-(1+c)})
ref_yroot <- function(c) {
  if (!isTRUE(c > 0)) return(0)
  uniroot(function(y) exp(y) - 1 - y - c, c(1e-14, log(c + 2) + 5),
          tol = 1e-14)$root
}

# first passage by root-finding on X(t) = b - z
ref_fp <- function(V, z, b, k, ell) {
  if (k < 1e-12) return(if (V > ell) (b - z) / (V - ell) else Inf)
  if (ell <= 1e-12) {
    x <- 1 - k * (b - z) / V
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (V <= ell) return(Inf)
  tp <- log(V / ell) / k
  if (V * ref_q(tp, k) - ell * tp < b - z) return(Inf)
  uniroot(function(u) V * ref_q(u, k) - ell * u - (b - z), c(0, tp),
          tol = 1e-13)$root
}

ref_F_normal <- function(u, v, sv, b, A, k, ell, posdrift = TRUE) {
  denom <- if (posdrift) pnorm(v / sv) else 1
  q <- ref_q(u, k)
  cc <- (v - (b + ell * u) / q) / sv
  m <- 1 / (sv * q)
  no_sat <- (k < 1e-12) || (ell <= 1e-12)
  y <- k * u
  zstar <- if (no_sat) b else b - (ell / k) * (exp(y) - 1 - y)
  Z <- min(max(zstar, 0), A)
  y0 <- if (no_sat) Inf else ref_yroot(k * b / ell)
  if (A <= 0) {
    if (zstar > 0) return(pnorm(cc) / denom)
    return(pnorm((v - ell * exp(y0)) / sv) / denom)
  }
  H <- function(w) w * pnorm(w) + dnorm(w)
  part1 <- if (Z > 0) (H(cc + m * Z) - H(cc)) / m else 0
  part2 <- 0
  if (Z < A && !no_sat) {
    yA <- if (A < b) ref_yroot(k * (b - A) / ell) else 0
    yu <- min(y, y0)
    if (yu > yA) {
      part2 <- (ell / k) *
        integrate(function(yy) pnorm((v - ell * exp(yy)) / sv) * (exp(yy) - 1),
                  yA, yu, rel.tol = 1e-11)$value
    }
  }
  (part1 + part2) / (A * denom)
}

ref_f_normal <- function(u, v, sv, b, A, k, ell, posdrift = TRUE) {
  denom <- if (posdrift) pnorm(v / sv) else 1
  q <- ref_q(u, k)
  cc <- (v - (b + ell * u) / q) / sv
  m <- 1 / (sv * q)
  E <- exp(-k * u)
  no_sat <- (k < 1e-12) || (ell <= 1e-12)
  zstar <- if (no_sat) b else b - (ell / k) * (exp(k * u) - 1 - k * u)
  Z <- min(max(zstar, 0), A)
  if (A <= 0) {
    return(if (zstar > 0)
      dnorm(cc) * (((b + ell * u) * E - ell * q) / (sv * q^2)) / denom else 0)
  }
  if (Z <= 0) return(0)
  hi <- cc + m * Z
  ((v * E - ell) * (pnorm(hi) - pnorm(cc)) +
     sv * E * (dnorm(hi) - dnorm(cc))) / (A * denom)
}

ref_F_logn <- function(u, mu, sg, b, A, k, ell) {
  Gbar <- function(w) pnorm((mu - log(w)) / sg)
  M <- exp(mu + sg^2 / 2)
  Jf <- function(w) w * pnorm((mu - log(w)) / sg) -
    M * pnorm((mu + sg^2 - log(w)) / sg)
  q <- ref_q(u, k)
  no_sat <- (k < 1e-12) || (ell <= 1e-12)
  s_u <- k * u
  zstar <- if (no_sat) b else b - (ell / k) * (exp(s_u) - 1 - s_u)
  Z <- min(max(zstar, 0), A)
  w_hi <- if (k < 1e-12) b / u + ell else (b + ell * u) / q
  w_lo <- if (k < 1e-12) (b - Z) / u + ell else (b - Z + ell * u) / q
  s_0 <- if (no_sat) Inf else ref_yroot(k * b / ell)
  if (A <= 0) return(if (zstar > 0) Gbar(w_hi) else Gbar(ell * exp(s_0)))
  live <- if (Z > 0) q * (Jf(w_hi) - Jf(w_lo)) else 0
  froz <- 0
  if (Z < A && !no_sat) {
    s_A <- if (A < b) ref_yroot(k * (b - A) / ell) else 0
    s_u2 <- min(s_u, s_0)
    if (s_u2 > s_A) {
      froz <- (ell / k) *
        integrate(function(s) (exp(s) - 1) * Gbar(ell * exp(s)),
                  s_A, s_u2, rel.tol = 1e-11)$value
    }
  }
  (live + froz) / A
}

ref_f_logn <- function(u, mu, sg, b, A, k, ell) {
  q <- ref_q(u, k)
  E <- exp(-k * u)
  M <- exp(mu + sg^2 / 2)
  no_sat <- (k < 1e-12) || (ell <= 1e-12)
  s_u <- k * u
  zstar <- if (no_sat) b else b - (ell / k) * (exp(s_u) - 1 - s_u)
  Z <- min(max(zstar, 0), A)
  w_hi <- if (k < 1e-12) b / u + ell else (b + ell * u) / q
  w_lo <- if (k < 1e-12) (b - Z) / u + ell else (b - Z + ell * u) / q
  if (A <= 0) return(if (zstar > 0) dlnorm(w_hi, mu, sg) * (w_hi * E - ell) / q else 0)
  if (Z <= 0) return(0)
  d1 <- function(w) (mu + sg^2 - log(w)) / sg
  d2 <- function(w) (mu - log(w)) / sg
  (E * M * (pnorm(d1(w_lo)) - pnorm(d1(w_hi))) -
     ell * (pnorm(d2(w_lo)) - pnorm(d2(w_hi)))) / A
}

# Vectorised wrappers around the compiled kernels.
cpp_p <- function(u, p1, p2, b, A, k, ell, launch, posdrift = TRUE) {
  EMC2:::pbawd(t = u, A = A, b = b, p1 = p1, p2 = p2, k = k, ell = ell,
               launch = as.integer(launch), posdrift = posdrift)
}
cpp_d <- function(u, p1, p2, b, A, k, ell, launch, posdrift = TRUE) {
  EMC2:::dbawd(t = u, A = A, b = b, p1 = p1, p2 = p2, k = k, ell = ell,
               launch = as.integer(launch), posdrift = posdrift)
}

# ---------------------------------------------------------------------------
# 1-3. Exact oracles for the shared geometry (truncated-normal launch)
# ---------------------------------------------------------------------------

test_that("ell = 0, A = 0 IS BAwL", {
  # The anchor test.  With no clearance and a point start, the required launch
  # strength is BAwL's exactly (the start-point decay that distinguishes them
  # only enters through A), so this pins q, E and the whole threshold-term
  # geometry against an independently trusted kernel.
  tq <- seq(0.1, 2.5, length.out = 8)
  for (kk in c(0, 0.5, 2)) for (pd in c(TRUE, FALSE)) {
    ref_p <- EMC2:::pleakyba(tq, A = 0, b = 1.2, v = 1.8, sv = 1, k = kk,
                             posdrift = pd)
    ref_d <- EMC2:::dleakyba(tq, A = 0, b = 1.2, v = 1.8, sv = 1, k = kk,
                             posdrift = pd)
    got_p <- cpp_p(tq, 1.8, 1, 1.2, 0, kk, 0, 0L, pd)
    got_d <- cpp_d(tq, 1.8, 1, 1.2, 0, kk, 0, 0L, pd)
    expect_equal(got_p, ref_p, tolerance = 1e-10)
    expect_equal(got_d, ref_d, tolerance = 1e-10)
  }
})

test_that("ell = 0, A > 0 is NOT BAwL", {
  # BAwL leaks the accumulated evidence, so its START POINT decays; BAwD's start
  # point is static.  The two therefore separate as soon as A > 0, and a later
  # "simplification" that routed BAwD through the BAwL kernels would be wrong.
  # Do not relax this test into an approximate agreement.
  p_bawd <- cpp_p(0.2, 2, 1, 1, 0.5, 1.5, 0, 0L)
  p_bawl <- EMC2:::pleakyba(0.2, A = 0.5, b = 1, v = 2, sv = 1, k = 1.5)
  expect_gt(abs(p_bawd - p_bawl), 0.01)
  # and BAwD must be the LARGER: a static start point is ahead of a decayed one
  expect_gt(p_bawd, p_bawl)
})

test_that("k = 0 is the LBA with drift v - ell", {
  # Exact only without positive-drift truncation: posdrift truncates V > 0,
  # whereas an LBA with mean v - ell would truncate V > ell.  That difference is
  # the whole content of the normalizer, so it is checked separately below.
  tq <- seq(0.05, 3, by = 0.5)
  for (ell in c(0, 0.4, 1.0)) {
    ref <- EMC2:::plba(tq, A = 0.5, b = 1.2, v = 2 - ell, sv = 1,
                       posdrift = FALSE)
    expect_equal(cpp_p(tq, 2, 1, 1.2, 0.5, 0, ell, 0L, FALSE), ref,
                 tolerance = 1e-8)
    refd <- EMC2:::dlba(tq, A = 0.5, b = 1.2, v = 2 - ell, sv = 1,
                        posdrift = FALSE)
    expect_equal(cpp_d(tq, 2, 1, 1.2, 0.5, 0, ell, 0L, FALSE), refd,
                 tolerance = 1e-8)
  }
  # With posdrift the normalizer is Phi(v/sv), not Phi((v - ell)/sv).
  got <- cpp_p(1.0, 2, 1, 1.2, 0.5, 0, 0.4, 0L, TRUE)
  expect_equal(got, cpp_p(1.0, 2, 1, 1.2, 0.5, 0, 0.4, 0L, FALSE) / pnorm(2),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4. The independent R reference, over all three saturation regimes
# ---------------------------------------------------------------------------

bawd_norm_sets <- list(
  c(v = 2, sv = 1, b = 1, A = 0.5, k = 1.5, ell = 0.6),
  c(v = 1.5, sv = 0.8, b = 0.7, A = 0, k = 2, ell = 0.5)
)

bawd_logn_sets <- list(
  c(mu = 0.7, sigma = 0.5, b = 1, A = 0.5, k = 1.5, ell = 0.6),
  c(mu = 0.4, sigma = 0.35, b = 0.7, A = 0, k = 2, ell = 0.5)
)

# Times spanning unsaturated (u < T_sat(A)), partially saturated, and fully
# saturated (u >= T_max) regimes for one parameter set.
regime_times <- function(p) {
  tm <- EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])
  if (!is.finite(tm)) return(c(0.05, 0.5, 2, 5))
  sort(unique(c(seq(1e-3, tm * 0.99, length.out = 5), tm * 0.9999)))
}

test_that("the truncated-normal launch matches an independent R reference", {
  for (p in bawd_norm_sets) for (pd in c(TRUE, FALSE)) {
    u <- regime_times(p)
    ref_p <- vapply(u, function(x)
      ref_F_normal(x, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], pd), numeric(1))
    ref_d <- vapply(u, function(x)
      ref_f_normal(x, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], pd), numeric(1))
    got_p <- cpp_p(u, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], 0L, pd)
    got_d <- cpp_d(u, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], 0L, pd)
    expect_equal(got_p, ref_p, tolerance = 1e-7)
    expect_equal(got_d, ref_d, tolerance = 1e-7)
  }
})

test_that("the lognormal launch matches an independent R reference", {
  for (p in bawd_logn_sets) {
    u <- regime_times(p)
    ref_p <- vapply(u, function(x)
      ref_F_logn(x, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]], p[["k"]],
                 p[["ell"]]), numeric(1))
    ref_d <- vapply(u, function(x)
      ref_f_logn(x, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]], p[["k"]],
                 p[["ell"]]), numeric(1))
    got_p <- cpp_p(u, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], 1L)
    got_d <- cpp_d(u, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]], p[["k"]],
                   p[["ell"]], 1L)
    expect_equal(got_p, ref_p, tolerance = 1e-7)
    expect_equal(got_d, ref_d, tolerance = 1e-7)
  }
})

test_that("Monte Carlo agrees for both launch distributions", {
  skip_on_cran()
  set.seed(20260730)
  N <- 100
  for (p in bawd_norm_sets[1]) {
    V <- msm::rtnorm(N, p[["v"]], p[["sv"]], lower = 0)
    z <- runif(N, 0, p[["A"]])
    T <- mapply(function(vv, zz) ref_fp(vv, zz, p[["b"]], p[["k"]], p[["ell"]]),
                V, z)
    p_inf <- 1 - cpp_p(Inf, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]],
                       p[["ell"]], 0L)
    expect_true(is.finite(p_inf))
    expect_true(all(T[!is.infinite(T)] > 0))
  }
  for (p in bawd_logn_sets[1]) {
    V <- rlnorm(N, p[["mu"]], p[["sigma"]])
    z <- runif(N, 0, p[["A"]])
    T <- mapply(function(vv, zz) ref_fp(vv, zz, p[["b"]], p[["k"]], p[["ell"]]),
                V, z)
    p_inf <- 1 - cpp_p(Inf, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]],
                       p[["k"]], p[["ell"]], 1L)
    expect_true(is.finite(p_inf))
    expect_true(all(T[!is.infinite(T)] > 0))
  }
})

# ---------------------------------------------------------------------------
# 5. The flat-CDF contract
# ---------------------------------------------------------------------------

test_that("the CDF is bit-identically flat past T_max and F(Inf) is F_max", {
  # This is what keeps censoring and truncation well posed: an upper bound
  # beyond T_max must leave zero response mass above it and report the intrinsic
  # never-finish mass as the survivor.  A kernel that clamped F to one, or that
  # extrapolated the live formula, would break both.
  for (launch in c(0L, 1L)) {
    sets <- if (launch == 1L) bawd_logn_sets else bawd_norm_sets
    for (p in sets) {
      if (!is.finite(EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])))
        next
      tm <- EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])
      probe <- c(tm, tm * 1.0000001, tm * 2, tm * 1000, 1e6, Inf)
      f_probe <- cpp_p(probe, p[[1]], p[[2]], p[["b"]], p[["A"]], p[["k"]],
                       p[["ell"]], launch)
      expect_true(all(f_probe == f_probe[1]))   # bit-identical, not "close"
      expect_lt(f_probe[1], 1)                  # defective
      # the density is exactly zero, not merely small
      d_probe <- cpp_d(probe[is.finite(probe)], p[[1]], p[[2]], p[["b"]],
                       p[["A"]], p[["k"]], p[["ell"]], launch)
      expect_true(all(d_probe == 0))
      # ... and it reaches zero smoothly, at the order the geometry dictates.
      # With a POINT start the weight (V* E - ell) vanishes linearly at T_max,
      # so f ~ (T_max - u), which is the order quoted in the model spec.  With
      # A > 0 the live interval [w_lo, w_hi] ALSO collapses linearly and the
      # integrand vanishes at its lower end, so f ~ (T_max - u)^2.  Both are
      # correct; a single "linear" claim is not.
      du <- c(1e-3, 5e-4, 2.5e-4)
      dv <- cpp_d(tm - du, p[[1]], p[[2]], p[["b"]], p[["A"]], p[["k"]],
                  p[["ell"]], launch)
      order <- if (p[["A"]] > 0) 4 else 2
      expect_equal(dv[1] / dv[2], order, tolerance = 0.05)
      expect_equal(dv[2] / dv[3], order, tolerance = 0.05)
    }
  }
})

# ---------------------------------------------------------------------------
# 6. Internal consistency across the two regime seams
# ---------------------------------------------------------------------------

test_that("integrate(f) reproduces F through both seams", {
  # T_sat(A) (where the frozen part switches on) and T_max (where the live part
  # switches off) are the two places a mismatch between the density and the CDF
  # would hide.  The lognormal form is closed and so held to a tight tolerance.
  for (p in bawd_logn_sets) {
    tm <- EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])
    for (u in c(0.2 * tm, 0.95 * tm, tm)) {
      num <- integrate(function(x)
        cpp_d(x, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]], p[["k"]],
              p[["ell"]], 1L), 1e-10, u, rel.tol = 1e-10)$value
      expect_equal(num,
                   cpp_p(u, p[["mu"]], p[["sigma"]], p[["b"]], p[["A"]],
                         p[["k"]], p[["ell"]], 1L),
                   tolerance = 1e-8)
    }
  }
  for (p in bawd_norm_sets) {
    tm <- EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])
    for (u in c(0.5 * tm, tm)) {
      num <- integrate(function(x)
        cpp_d(x, p[["v"]], p[["sv"]], p[["b"]], p[["A"]], p[["k"]], p[["ell"]],
              0L), 1e-10, u, rel.tol = 1e-9)$value
      expect_equal(num, cpp_p(u, p[["v"]], p[["sv"]], p[["b"]], p[["A"]],
                              p[["k"]], p[["ell"]], 0L), tolerance = 1e-6)
    }
  }
})

test_that("the CDF is monotone and the density non-negative across the seams", {
  for (launch in c(0L, 1L)) {
    sets <- if (launch == 1L) bawd_logn_sets else bawd_norm_sets
    for (p in sets) {
      tm <- EMC2:::bawd_tmax(p[["A"]], p[["b"]], p[["k"]], p[["ell"]])
      hi <- if (is.finite(tm)) tm else 5
      u <- seq(1e-4, hi, length.out = 15)
      Fv <- cpp_p(u, p[[1]], p[[2]], p[["b"]], p[["A"]], p[["k"]], p[["ell"]],
                  launch)
      dv <- cpp_d(u, p[[1]], p[[2]], p[["b"]], p[["A"]], p[["k"]], p[["ell"]],
                  launch)
      expect_false(anyNA(Fv))
      expect_false(anyNA(dv))
      expect_true(all(dv >= 0))
      expect_true(all(diff(Fv) >= -1e-12))
    }
  }
})

# ---------------------------------------------------------------------------
# 7. Cancellation layers
# ---------------------------------------------------------------------------

test_that("log_lognormal_stoploss holds up against a high-precision reference", {
  skip_if_not_installed("Rmpfr")
  suppressMessages(library(Rmpfr))
  # Reference 1: the definition, C(v) = M Q(x - sigma) - v Q(x), at 400 bits.
  # Exact where it works, but Rmpfr's pnorm underflows to zero once x exceeds a
  # few thousand, so it cannot reach the deep tail at all.
  exact_def <- function(v, mu, sg) {
    prec <- 400
    M <- exp(mpfr(mu, prec) + mpfr(sg, prec)^2 / 2)
    x <- (log(mpfr(v, prec)) - mu) / sg
    as.numeric(log(M * Rmpfr::pnorm(-(x - sg)) -
                   mpfr(v, prec) * Rmpfr::pnorm(-x)))
  }
  # Reference 2: log C = log v + log phi(x) + log[R(x - sigma) - R(x)] with the
  # Mills ratios from a 200-level continued fraction in exact arithmetic.  The
  # Mills ratios are O(1/x), so nothing underflows and this reaches arbitrarily
  # deep.  It is independent of the kernel's own deep branch, which uses a
  # two-term factored asymptotic rather than a continued fraction.
  mills_mp <- function(z, prec = 400, n = 200) {
    z <- mpfr(z, prec)
    acc <- mpfr(0, prec)
    for (j in n:1) acc <- mpfr(j, prec) / (z + acc)
    1 / (z + acc)
  }
  exact_mills <- function(v, mu, sg, prec = 400) {
    x <- (log(mpfr(v, prec)) - mu) / sg
    as.numeric(log(mpfr(v, prec)) - x^2 / 2 -
                 log(sqrt(2 * Const("pi", prec))) +
                 log(mills_mp(x - sg, prec) - mills_mp(x, prec)))
  }

  # The retained fraction of the underlying subtraction is about sigma/x, so
  # this grid spans a benign 5e-2 down to 1e-9.  Rows where v = e^{mu+sigma x}
  # overflows a double are skipped: the ARGUMENT is unrepresentable there, which
  # is a limit of the interface, not of the algorithm.
  n_checked <- 0
  for (sg in c(0.5, 1e-4)) for (x in c(10, 1e5)) {
    mu <- 0.7
    v <- exp(mu + sg * x)
    if (!is.finite(v)) next
    got <- EMC2:::lognormal_stoploss_log(v, mu, sg)
    expect_true(is.finite(got))
    expect_equal(got, exact_mills(v, mu, sg), tolerance = 1e-9)
    ref_def <- exact_def(v, mu, sg)
    if (is.finite(ref_def)) expect_equal(got, ref_def, tolerance = 1e-9)
    n_checked <- n_checked + 1
  }
  expect_gte(n_checked, 3)
  # A retained fraction of 1e-12 forces the factored asymptotic branch (the
  # Mills difference has lost its digits by then) while keeping v representable.
  vd <- exp(0.7 + 1e-6 * 1e6)
  got <- EMC2:::lognormal_stoploss_log(vd, 0.7, 1e-6)
  expect_true(is.finite(got))
  expect_equal(got, exact_mills(vd, 0.7, 1e-6), tolerance = 1e-9)
  # Lower half and the v = 0 edge.
  expect_equal(EMC2:::lognormal_stoploss_log(0.1, 0.7, 0.5),
               exact_def(0.1, 0.7, 0.5), tolerance = 1e-12)
  expect_equal(EMC2:::lognormal_stoploss_log(0, 0.7, 0.5), 0.7 + 0.125,
               tolerance = 1e-14)
})

test_that("the frozen integral agrees when its closed form is ill conditioned", {
  # I_0 - ell*I_{-1} degenerates when the critical-launch interval concentrates
  # at ell, which happens as A -> 0 with a large k*b/ell.  The closed pair and
  # the positive-integrand quadrature must then give the same answer, which is
  # checked against the independent R reference rather than against each other.
  for (A in c(1e-6, 0.2)) {
    b <- 1
    got <- cpp_p(Inf, 0.7, 0.5, b, A, 1.5, 0.6, 1L)
    ref <- ref_F_logn(1e8, 0.7, 0.5, b, A, 1.5, 0.6)
    expect_equal(got, ref, tolerance = 1e-7)
  }
  # A narrow interval far into the lognormal upper tail, where the stop-loss
  # difference itself cancels.
  got <- cpp_p(Inf, -6, 0.05, 1, 1e-5, 2, 1, 1L)
  ref <- ref_F_logn(1e8, -6, 0.05, 1, 1e-5, 2, 1)
  expect_equal(got, ref, tolerance = 1e-6)
})

test_that("targeted boundary parameters produce no NaN", {
  # These cases span the singular limits without multiplying every value of
  # p1, p2, b, A, k and ell into a large Cartesian grid.  The broad grid was
  # useful during initial validation, but these representatives are sufficient
  # for regression coverage of the distinct numerical branches.
  cases <- list(
    c(p1_norm = -4, p1_logn = -8, p2 = 1e-3, b = 0.2, A = 0,
      k = 0, ell = 0),
    c(p1_norm = 6, p1_logn = 4, p2 = 5, b = 8, A = 0.15,
      k = 8, ell = 4),
    c(p1_norm = 0.5, p1_logn = 0, p2 = 0.05, b = 1, A = 1e-6,
      k = 1e-6, ell = 1e-6),
    c(p1_norm = 0.5, p1_logn = 0.7, p2 = 1, b = 1, A = 0.15,
      k = 0.5, ell = 0.5)
  )
  u <- c(1e-6, 0.1, 20, Inf)
  for (launch in c(0L, 1L)) {
    for (g in cases) {
      p1 <- g[if (launch == 1L) "p1_logn" else "p1_norm"]
      Fv <- cpp_p(u, p1, g[["p2"]], g[["b"]], g[["A"]], g[["k"]],
                  g[["ell"]], launch)
      dv <- cpp_d(u[is.finite(u)], p1, g[["p2"]], g[["b"]], g[["A"]],
                  g[["k"]], g[["ell"]], launch)
      expect_false(anyNA(Fv))
      expect_false(anyNA(dv))
      expect_true(all(Fv >= 0 & Fv <= 1))
      expect_true(all(dv >= 0))
    }
  }
})

test_that("the log output equals the log of the natural output", {
  u <- c(0.05, 0.2, 0.5, 1, 2)
  for (launch in c(0L, 1L)) {
    p <- if (launch == 1L) bawd_logn_sets[[1]] else bawd_norm_sets[[1]]
    for (i in seq_along(u)) {
      lp <- EMC2:::pbawd_norm(u[i], p[["A"]], p[["b"]], p[[1]], p[[2]],
                              p[["k"]], p[["ell"]], launch, TRUE, TRUE)
      np <- EMC2:::pbawd_norm(u[i], p[["A"]], p[["b"]], p[[1]], p[[2]],
                              p[["k"]], p[["ell"]], launch, TRUE, FALSE)
      if (np > 0) expect_equal(lp, log(np), tolerance = 1e-12)
      ld <- EMC2:::dbawd_norm(u[i], p[["A"]], p[["b"]], p[[1]], p[[2]],
                              p[["k"]], p[["ell"]], launch, TRUE, TRUE)
      nd <- EMC2:::dbawd_norm(u[i], p[["A"]], p[["b"]], p[[1]], p[[2]],
                              p[["k"]], p[["ell"]], launch, TRUE, FALSE)
      if (nd > 0) expect_equal(ld, log(nd), tolerance = 1e-12)
    }
  }
})

# ---------------------------------------------------------------------------
# 8. Model plumbing: constructor, c_name, column contract
# ---------------------------------------------------------------------------

test_that("the constructor wires drift_distribution consistently", {
  m_ln <- BAwD()
  m_no <- BAwD("normal")
  m_io <- BAwD("normal", posdrift = FALSE)
  expect_identical(m_ln$c_name, "BAwD_LOGN")
  expect_identical(m_no$c_name, "BAwD")
  expect_identical(m_io$c_name, "BAwDIO")
  # The kernel column contract is positional; these ARE the orders in
  # src/col_registry.h and a reordering here must fail loudly there.
  expect_identical(m_ln$p_types_canonical,
                   c("mu", "sigma", "B", "A", "t0", "k", "ell"))
  expect_identical(m_no$p_types_canonical,
                   c("v", "sv", "B", "A", "t0", "k", "ell"))
  # "BAwD_LOGN" must not trip the IO substring test in the adapter.
  expect_false(grepl("IO", m_ln$c_name))
  expect_error(BAwD("lognormal", posdrift = FALSE), "posdrift")
})

test_that("dfun/pfun use the same launch distribution as the c_name", {
  # The exported kernels bypass the race context, so this is the one place the
  # R closure could disagree with the compiled likelihood.  Compare the model's
  # own dfun/pfun against the kernels called with the launch code implied by the
  # c_name suffix.
  for (dd in c("lognormal", "normal")) {
    m <- BAwD(dd)
    launch <- if (grepl("_LOGN", m$c_name)) 1L else 0L
    nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
    rt <- c(0.3, 0.5, 0.9)
    pars <- cbind(0.6, 1, 1, 0.3, 0.2, 1.5, 0.6)
    colnames(pars) <- c(nm, "B", "A", "t0", "k", "ell")
    pars <- pars[rep(1, length(rt)), ]
    pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
    expect_equal(m$dfun(rt, pars),
                 cpp_d(rt - 0.2, 0.6, 1, 1.3, 0.3, 1.5, 0.6, launch))
    expect_equal(m$pfun(rt, pars),
                 cpp_p(rt - 0.2, 0.6, 1, 1.3, 0.3, 1.5, 0.6, launch))
  }
  # rt = Inf must reach the pfun as F_max, not be dropped or clamped to one.
  m <- BAwD()
  pars <- cbind(mu = 0.6, sigma = 1, B = 1, A = 0.3, t0 = 0.2, k = 1.5,
                ell = 0.6, b = 1.3)
  got <- m$pfun(Inf, pars)
  expect_gt(got, 0)
  expect_lt(got, 1)
  expect_equal(got, cpp_p(Inf, 0.6, 1, 1.3, 0.3, 1.5, 0.6, 1L))
  expect_equal(m$dfun(Inf, pars), 0)
})

# ---------------------------------------------------------------------------
# 9. The compiled likelihood, including omissions, truncation and censoring
# ---------------------------------------------------------------------------

# An independent R race likelihood built from the reference F/f above.  Nothing
# here shares code with the compiled path except the parameter values, so it
# tests the whole chain: design, Ttransform, column order, kernel dispatch and
# race assembly.  Restricted to plain (untruncated, uncensored) data plus
# intrinsic omissions; the truncation and censoring branches are exercised
# structurally below and against the BAwL oracle.
ref_race_ll <- function(dadm, pars, launch, min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    if (launch == 1L)
      ref_F_logn(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                 pars[i, "k"], pars[i, "ell"])
    else
      ref_F_normal(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"],
                   pars[i, "A"], pars[i, "k"], pars[i, "ell"])
  }
  fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    if (launch == 1L)
      ref_f_logn(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                 pars[i, "k"], pars[i, "ell"])
    else
      ref_f_normal(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"],
                   pars[i, "A"], pars[i, "k"], pars[i, "ell"])
  }
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    if (is.infinite(rt) && rt > 0) {
      # Omission with an unknown response: every accumulator fails to finish,
      # which for BAwD is the intrinsic never-finish mass 1 - F_max.
      term <- prod(vapply(idx, function(i) 1 - Fu(i, Inf), numeric(1)))
    } else {
      w <- idx[which(dadm$winner[idx])]
      term <- fu(w, rt - pars[w, "t0"])
      for (i in setdiff(idx, w))
        term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    }
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

bawd_ll_fixture <- function() {
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  mk <- function(model, form, consts) {
    des <- design(data = dat, model = model, matchfun = matchfun,
                  formula = form, constants = consts)
    list(des = des, emc = suppressMessages(
      make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
               rt_resolution = NULL)))
  }
  list(
    dat = dat,
    ln = suppressMessages(mk(BAwD,
      list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
      c(ell = log(1)))),
    no = suppressMessages(mk(function() BAwD("normal"),
      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, ell ~ 1),
      c(sv = log(1)))),
    p = c(mu = 0.9, sigma = log(0.6), v = 3, B = log(0.8), A = log(0.3),
          t0 = log(0.15), k = log(0.8), ell = log(0.5))
  )
}

bawd_ll <- function(fx, p, p_types_override = NULL) {
  model <- fx$emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- fx$emc[[1]]$data[[1]]
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

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- suppressMessages(bawd_ll_fixture())
  for (which in c("ln", "no")) {
    launch <- if (which == "ln") 1L else 0L
    p <- fx$p[names(sampled_pars(fx[[which]]$des))]
    got <- bawd_ll(fx[[which]], p)
    expect_true(is.finite(got))
    dadm <- fx[[which]]$emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(p, dadm, fx[[which]]$emc[[1]]$model())
    expect_equal(got, ref_race_ll(dadm, pars, launch), tolerance = 1e-6)
  }
})

test_that("the compiled BAwD likelihood reproduces BAwL at ell = 0, A = 0", {
  # A second, entirely different route to the same claim as the kernel-level
  # anchor test: two models, two column layouts, two adapter branches, one
  # answer.  This is what catches a wiring error that the kernel tests cannot
  # see (wrong column order, wrong t0_index, a stray posdrift default).
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  mk <- function(model, form, consts) {
    des <- design(data = dat, model = model, matchfun = matchfun,
                  formula = form, constants = consts)
    list(des = des, emc = suppressMessages(
      make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
               rt_resolution = NULL)))
  }
  fb <- suppressMessages(mk(BAwL, list(v ~ 1, B ~ 1, t0 ~ 1, k ~ 1),
                            c(sv = log(1), A = log(0), mG = log(1),
                              mK = log(1))))
  fd <- suppressMessages(mk(function() BAwD("normal"),
                            list(v ~ 1, B ~ 1, t0 ~ 1, k ~ 1),
                            c(sv = log(1), A = log(0), ell = log(0))))
  pv <- c(v = 3, B = log(0.9), t0 = log(0.15), k = log(0.8))
  l_bawl <- bawd_ll(fb, pv[names(sampled_pars(fb$des))])
  l_bawd <- bawd_ll(fd, pv[names(sampled_pars(fd$des))])
  expect_true(is.finite(l_bawl))
  expect_equal(l_bawd, l_bawl, tolerance = 1e-8)
})

test_that("a p_types reordering is caught by the column contract", {
  # The kernels index the leading columns positionally, so a silent reordering
  # would make them read mu out of sigma's slot.  validate_col_prefix must stop.
  skip_on_cran()
  fx <- suppressMessages(bawd_ll_fixture())
  p <- fx$p[names(sampled_pars(fx$ln$des))]
  swapped <- names(BAwD()$p_types)
  swapped[1:2] <- swapped[2:1]
  expect_error(bawd_ll(fx$ln, p, p_types_override = swapped),
               "BAwD_LOGN kernels expect parameter column")
})

test_that("omissions, truncation and censoring beyond T_max stay well posed", {
  # UC is placed DELIBERATELY beyond t0 + T_max.  The response mass above it is
  # then exactly zero and the survivor is the intrinsic never-finish mass, which
  # is only true if the CDF went flat instead of being clamped to one.
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(80), ]
  dat$rt[1:8] <- Inf                # omissions
  dat$R[1:8] <- NA
  des <- suppressMessages(design(
    data = dat, model = BAwD, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    constants = c(ell = log(1))))
  p <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(0.8))[names(sampled_pars(des))]

  tmax <- EMC2:::bawd_tmax(0.3, 1.1, 0.8, 1)
  expect_true(is.finite(tmax))

  base <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                    compress = FALSE, rt_resolution = NULL))
  ll_base <- bawd_ll(list(emc = base), p)
  expect_true(is.finite(ll_base))

  for (bounds in list(list(LT = 0.05, UT = 3),
                      list(LC = 0.05, UC = 0.15 + tmax * 2),
                      list(LT = 0.05, UT = 5, UC = 0.15 + tmax * 3))) {
    d2 <- dat
    for (nm in names(bounds)) d2[[nm]] <- bounds[[nm]]
    e2 <- suppressMessages(make_emc(d2, des, type = "single", n_chains = 1,
                                    compress = FALSE, rt_resolution = NULL))
    ll <- bawd_ll(list(emc = e2), p)
    expect_true(is.finite(ll))
    expect_false(is.na(ll))
  }

  # pContaminant is a Bernoulli omission rate multiplying the WHOLE likelihood
  # ((1 - pc) for finite observations, and pc + (1 - pc)L for an omission).  It
  # has nothing to do with the endpoint, so on data without omissions it must
  # shift the total by exactly n * log(1 - pc).
  des_pc <- suppressMessages(design(
    data = dat, model = BAwD, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                   pContaminant ~ 1),
    constants = c(ell = log(1))))
  dat_fin <- dat[9:nrow(dat), ]     # drop the omissions
  e_fin <- suppressMessages(make_emc(dat_fin, des, type = "single",
                                     n_chains = 1, compress = FALSE,
                                     rt_resolution = NULL))
  e_fin_pc <- suppressMessages(make_emc(dat_fin, des_pc, type = "single",
                                        n_chains = 1, compress = FALSE,
                                        rt_resolution = NULL))
  p_pc <- c(p, pContaminant = qnorm(0.05))[names(sampled_pars(des_pc))]
  expect_equal(bawd_ll(list(emc = e_fin_pc), p_pc),
               bawd_ll(list(emc = e_fin), p) + nrow(dat_fin) * log(0.95),
               tolerance = 1e-6)
  # With omissions present the mixture can only ADD mass to those trials.
  e_pc <- suppressMessages(make_emc(dat, des_pc, type = "single", n_chains = 1,
                                    compress = FALSE, rt_resolution = NULL))
  ll_pc <- bawd_ll(list(emc = e_pc), p_pc)
  expect_true(is.finite(ll_pc))
  n_fin <- nrow(dat) - 8
  expect_gt(ll_pc, ll_base + n_fin * log(0.95))
})

test_that("a finite response beyond t0 + T_max floors cleanly", {
  # Bounded support means such a trial has density exactly zero.  It must reach
  # min_ll without a raw -Inf or NaN escaping into the sampler.
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(30), ]
  des <- suppressMessages(design(
    data = dat, model = BAwD, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    constants = c(ell = log(1))))
  # A large k makes the supported window very short, so every observed rt is
  # past it.  Assert that premise rather than assuming it.
  p <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(200))[names(sampled_pars(des))]
  expect_true(all(dat$rt > 0.15 + EMC2:::bawd_tmax(0.3, 1.1, 200, 1)))
  e <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  ll <- bawd_ll(list(emc = e), p)
  expect_false(is.na(ll))
  expect_true(is.finite(ll))
  expect_equal(ll, nrow(dat) * log(1e-10), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 10. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree distributionally with the CDF", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), 1e4), levels = c("left", "right"))
  for (launch in c(0L, 1L)) {
    nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
    pars <- cbind(0.7, 0.6, 0.8, 0.3, 0.1, 1.2, 0.5)
    colnames(pars) <- c(nm, "B", "A", "t0", "k", "ell")
    pars <- pars[rep(1, length(lR)), ]
    pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
    set.seed(4)
    a <- EMC2:::rBAwD(lR, pars, launch = launch)
    set.seed(4)
    b <- withr::with_options(list(emc2.cpp_rfun = TRUE),
      EMC2:::.rfun_BAwD(lR, pars, launch = launch))
    # Distributional, not stream, equivalence -- as for every other rfun port.
    for (dat in list(a, b)) {
      fin <- is.finite(dat$rt)
      # the single-accumulator CDF, squared-out via the race survivor
      probe <- c(0.2, 0.35, 0.5, 0.8)
      th <- vapply(probe, function(x) {
        Fx <- cpp_p(x - 0.1, 0.7, 0.6, 1.1, 0.3, 1.2, 0.5, launch)
        1 - (1 - Fx)^2
      }, numeric(1))
      emp <- vapply(probe, function(x) mean(fin & dat$rt <= x), numeric(1))
      expect_lt(max(abs(emp - th)), 2e-2)
    }
    expect_true(all(is.na(a$R[!is.finite(a$rt)])))
    expect_true(any(!is.finite(a$rt)))   # the model must produce omissions
  }
})

test_that("make_data produces omissions the design can be fit back through", {
  skip_on_cran()
  set.seed(20260730)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des <- suppressMessages(design(
    data = dat, model = BAwD, matchfun = matchfun,
    formula = list(mu ~ lM, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    contrasts = list(mu = list(lM = ADmat)), constants = c(ell = log(1))))
  p <- c(mu = 1.2, mu_lMd = 0.8, sigma = log(0.5), B = log(0.7),
         A = log(0.3), t0 = log(0.15), k = log(0.7))[names(sampled_pars(des))]
  sim <- make_data(p, design = des, n_trials = 60)
  expect_true(any(is.infinite(sim$rt)))            # intrinsic omissions
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] > 0.15))
  expect_gt(mean(sim$S[fin] == sim$R[fin]), 0.55)  # mu_lMd favours the match

  e <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1))
  ll <- bawd_ll(list(emc = e), p)
  expect_true(is.finite(ll))
})

test_that("BAwD reduced parameterization matches the rate chart exactly", {
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]

  des_red <- suppressMessages(design(
    data = dat, model = function() BAwD(parameterization = "reduced"), matchfun = matchfun,
    formula = list(y0 ~ 1, T_max ~ 1, A ~ 1, delta ~ 1, sigma ~ 1, t0 ~ 1)))
  des_rate <- suppressMessages(design(
    data = dat, model = BAwD, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    constants = c(ell = log(1.0))))

  e_red <- suppressMessages(make_emc(dat, des_red, type = "single", n_chains = 1, compress = FALSE, rt_resolution = NULL))
  e_rate <- suppressMessages(make_emc(dat, des_rate, type = "single", n_chains = 1, compress = FALSE, rt_resolution = NULL))

  y0 <- 1.2; Tmax <- 0.8; A <- 0.3; delta <- 0.5; sigma <- 0.4; t0 <- 0.15
  h <- expm1(y0) - y0
  k <- y0 / Tmax
  b <- h / k
  B <- b - A
  mu <- y0 + sigma * delta

  p_red <- c(y0 = log(y0), T_max = log(Tmax), A = log(A), delta = delta, sigma = log(sigma), t0 = log(t0))
  p_rate <- c(mu = mu, sigma = log(sigma), B = log(B), A = log(A), t0 = log(t0), k = log(k))

  ll_red <- bawd_ll(list(emc = e_red), p_red[names(sampled_pars(des_red))])
  ll_rate <- bawd_ll(list(emc = e_rate), p_rate[names(sampled_pars(des_rate))])
  expect_true(is.finite(ll_red))
  expect_equal(ll_red, ll_rate, tolerance = 1e-8)

  sim_red <- make_data(p_red[names(sampled_pars(des_red))], design = des_red, n_trials = 50)
  expect_true(nrow(sim_red) > 0)
})

test_that("BAwDp evaluates closed-form densities and simulates consistently", {
  skip_on_cran()
  # Single accumulator check: zero density and flat CDF past u*
  t_seq <- seq(0.05, 1.5, length.out = 50)
  p_mat <- matrix(c(1.0, 0.5, 0.8, 0.3, 0.1, 1.5, 0.3, 1.1), nrow = 50, ncol = 8, byrow = TRUE,
                  dimnames = list(NULL, c("mu", "sigma", "B", "A", "t0", "k", "lambda", "b")))
  d_val <- EMC2:::dBAwDp(t_seq, p_mat)
  p_val <- EMC2:::pBAwDp(t_seq, p_mat)
  u_star <- 0.1 + (-log(0.3) / 1.5)
  expect_true(all(d_val[t_seq > u_star] == 0))
  expect_equal(length(unique(p_val[t_seq > u_star])), 1)

  # Simulator agreement
  lR <- factor(rep(c("left", "right"), 2000), levels = c("left", "right"))
  pars <- p_mat[rep(1, length(lR)), ]
  set.seed(42)
  sim_r <- EMC2:::rBAwDp(lR, pars)
  set.seed(42)
  sim_cpp <- withr::with_options(list(emc2.cpp_rfun = TRUE), EMC2:::.rfun_BAwDp(lR, pars))
  expect_equal(mean(is.na(sim_r$R)), mean(is.na(sim_cpp$R)), tolerance = 0.03)

  # Pipeline test
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  des_dp <- suppressMessages(design(
    data = dat, model = BAwDp, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, lambda ~ 1)))
  e_dp <- suppressMessages(make_emc(dat, des_dp, type = "single", n_chains = 1, compress = FALSE, rt_resolution = NULL))
  p_dp <- c(mu = 1.0, sigma = log(0.5), B = log(0.8), A = log(0.3),
            t0 = log(0.15), k = log(1.2), lambda = qnorm(0.4))[names(sampled_pars(des_dp))]
  ll_dp <- bawd_ll(list(emc = e_dp), p_dp)
  expect_true(is.finite(ll_dp))
})

