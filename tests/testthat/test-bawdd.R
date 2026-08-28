test_that("BAwDD exposes alpha = 1/rho with the exponential boundary", {
  m <- EMC2::BAwDD()
  expect_equal(m$c_name, "BAwDD_LOGN")
  expect_true("alpha" %in% m$p_types_canonical)
  expect_false("rho" %in% m$p_types_canonical)
  expect_true("ell" %in% m$p_types_canonical)
  expect_equal(unname(m$bound$minmax[, "alpha"]), c(0, Inf))
  expect_equal(unname(m$bound$exception["alpha"]), 0)
  expect_equal(EMC2::BAwDD("normal", posdrift = FALSE)$c_name, "BAwDDIO")
})

test_that("BAwDD closed forms reuse the BAwD kernel under reciprocal rho", {
  alpha <- c(0, 0.5, 1)
  rho <- c(Inf, 2, 1)
  pars <- data.frame(mean = 1, cv = 0.4, b = 1.2, A = 0.2,
                     t0 = 0, k = 0.5, ell = 0, alpha = alpha)
  rt <- rep(0.8, length(alpha))
  got_d <- EMC2:::dBAwDD(rt, pars)
  got_p <- EMC2:::pBAwDD(rt, pars)
  ref_d <- vapply(rho, function(r)
    EMC2:::dbawd(0.8, A = 0.2, b = 1.2, p1 = 1, p2 = 0.4,
                 k = 0.5, ell = 0, rho = r), numeric(1))
  ref_p <- vapply(rho, function(r)
    EMC2:::pbawd(0.8, A = 0.2, b = 1.2, p1 = 1, p2 = 0.4,
                 k = 0.5, ell = 0, rho = r), numeric(1))
  expect_equal(got_d, ref_d, tolerance = 1e-10)
  expect_equal(got_p, ref_p, tolerance = 1e-10)
  expect_equal(got_d[3], ref_d[3], tolerance = 1e-12)
})

test_that("BAwDD compiled simulator accepts row-wise alpha", {
  pars <- as.matrix(data.frame(mean = c(1, 1), cv = c(0.4, 0.4),
                               B = c(1, 1), A = c(0.2, 0.2),
                               t0 = c(0.05, 0.05), k = c(0.5, 0.5),
                               ell = c(0, 0), alpha = c(0, 1), eta = c(0, 0),
                               pContaminant = c(0, 0), pGuess = c(0, 0),
                               b = c(1.2, 1.2)))
  out <- EMC2:::.rfun_BAwDD(factor(c("a", "b"), levels = c("a", "b")),
                            pars, ok = c(TRUE, TRUE), launch = 1L)
  expect_true(is.data.frame(out))
  expect_named(out, c("R", "rt"))
  expect_length(out$rt, 1)
})

# ---------------------------------------------------------------------------
# alpha > 1 (rho < 1) with ell > 0.  This combination used to be rejected by
# bawd_geometry() and by the simulator, which pinned every fit to alpha <= 1.
# ---------------------------------------------------------------------------

## Reference straight from the model definition:
##   h(t) = (1 + alpha k t)^(-1/alpha),  X(t) = z + V Q(t) - ell t
##   CDF(u) = (1/A) int_0^A P(V >= min_{s<=u} (b - z + ell s)/Q(s)) dz
.bawdd_ref_Q <- function(t, k, alpha) {
  if (alpha == 0) return((1 - exp(-k * t)) / k)
  cc <- alpha * k
  if (abs(alpha - 1) < 1e-12) return(log1p(cc * t) / cc)
  p <- 1 - 1 / alpha
  ((1 + cc * t)^p - 1) / (cc * p)
}

## Every start point has peaked by T_max, and the earliest (z = A) peaks at
## T_sat_A, so only t > T_sat_A exercises the frozen branch -- the half of the
## closed form that carries the sign of rho - 1.  Both edges move by more than
## an order of magnitude across the alpha grid, so the times have to be chosen
## per alpha rather than fixed.
.bawdd_ref_times <- function(A, b, k, ell, alpha) {
  rho <- if (alpha > 0 && is.finite(alpha)) 1 / alpha else Inf
  t_sat <- EMC2:::bawd_tmax(A = 0, b = b - A, k = k, ell = ell, gamma = 0, rho = rho)
  t_max <- EMC2:::bawd_tmax(A = A, b = b, k = k, ell = ell, gamma = 0, rho = rho)
  c(live = 0.5 * t_sat, partial = 0.5 * (t_sat + t_max), frozen = 1.5 * t_max)
}

.bawdd_ref_cdf <- function(u, A, b, k, ell, alpha, Sbar) {
  vreq <- function(z) {
    f <- function(s) (b - z + ell * s) / .bawdd_ref_Q(s, k, alpha)
    s <- exp(seq(log(u) - 30, log(u), length.out = 20001))
    v <- f(s)
    i <- which.min(v)
    lo <- s[max(1, i - 1)]
    hi <- s[min(length(s), i + 1)]
    if (hi > lo) min(v[i], stats::optimize(f, c(lo, hi))$objective) else v[i]
  }
  g <- function(z) vapply(z, function(zz) Sbar(vreq(zz)), numeric(1))
  stats::integrate(g, 0, A, rel.tol = 1e-10, subdivisions = 500L)$value / A
}

test_that("alpha > 1 with ell > 0 matches the reference CDF for every launch", {
  skip_on_cran()
  A <- 0.4; b <- 1.2; k <- 0.9; ell <- 0.25
  cases <- list(
    list(launch = 1L, p1 = 0.3, p2 = 0.5,
         Sbar = function(v) stats::pnorm((log(v) - 0.3) / 0.5, lower.tail = FALSE)),
    list(launch = 3L, p1 = 2, p2 = 1.5,
         Sbar = local({
           scale <- 1.5 / gamma(1 + 1 / 2)
           function(v) exp(-(v / scale)^2)
         })),
    list(launch = 0L, p1 = 1.4, p2 = 0.6,
         Sbar = function(v) stats::pnorm((1.4 - v) / 0.6) / stats::pnorm(1.4 / 0.6))
  )
  for (cs in cases) for (a in c(1.1, 1.5, 3, 10, 40)) {
    for (tt in c(0.3, 0.8, 2.2, .bawdd_ref_times(A, b, k, ell, a))) {
      got <- EMC2:::pbawdd(tt, A = A, b = b, p1 = cs$p1, p2 = cs$p2, k = k, ell = ell,
                    alpha = a, launch = cs$launch)
      ref <- .bawdd_ref_cdf(tt, A, b, k, ell, a, cs$Sbar)
      expect_equal(got, ref, tolerance = 1e-8)
    }
    # Nothing can still be running past T_max, so the CDF is flat from there on
    # and the omission mass reported at t = Inf has to agree with it exactly.
    t_max <- EMC2:::bawd_tmax(A = A, b = b, k = k, ell = ell, gamma = 0, rho = 1 / a)
    expect_equal(EMC2:::pbawdd(Inf, A = A, b = b, p1 = cs$p1, p2 = cs$p2, k = k,
                               ell = ell, alpha = a, launch = cs$launch),
                 EMC2:::pbawdd(3 * t_max, A = A, b = b, p1 = cs$p1, p2 = cs$p2,
                               k = k, ell = ell, alpha = a, launch = cs$launch),
                 tolerance = 1e-10)
  }
})

test_that("alpha < 1 and the rho = 1 seam are unchanged by the rho < 1 branch", {
  skip_on_cran()
  A <- 0.4; b <- 1.2; k <- 0.9; ell <- 0.25
  Sbar <- function(v) stats::pnorm((log(v) - 0.3) / 0.5, lower.tail = FALSE)
  for (a in c(0.25, 0.5, 0.9, 1)) {
    for (tt in c(0.8, .bawdd_ref_times(A, b, k, ell, a))) {
      got <- EMC2:::pbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k, ell = ell, alpha = a)
      expect_equal(got, .bawdd_ref_cdf(tt, A, b, k, ell, a, Sbar), tolerance = 1e-8)
    }
  }
  # The closed form switches branch at rho = 1; it must not leave a step there.
  # t = 0.8 is still pre-saturation, so run the sweep past T_sat_A as well --
  # the frozen branch is the only place the sign of rho - 1 shows up at all.
  eps <- 10^-(4:11)
  P <- function(a, tt) EMC2:::pbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5,
                                     k = k, ell = ell, alpha = a)
  for (tt in c(0.8, 4.1, 12)) {
    lo <- vapply(1 - eps, P, numeric(1), tt = tt)
    hi <- vapply(1 + eps, P, numeric(1), tt = tt)
    at1 <- P(1, tt)
    expect_equal(lo[length(lo)], at1, tolerance = 1e-10)
    expect_equal(hi[length(hi)], at1, tolerance = 1e-10)
    # Both sides approach the seam monotonically from their own direction.
    expect_true(all(diff(lo) > 0) && all(diff(hi) < 0))
  }
})

test_that("the density is the derivative of the CDF for alpha > 1 with ell > 0", {
  skip_on_cran()
  A <- 0.4; b <- 1.2; k <- 0.9; ell <- 0.25
  P <- function(tt, a) EMC2:::pbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5,
                                     k = k, ell = ell, alpha = a)
  D <- function(tt, a) EMC2:::dbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5,
                                     k = k, ell = ell, alpha = a)
  # Fourth-order stencil on a t-relative step.  Late in the tail the density is
  # ~1e-8 against a CDF of ~1, so a central difference on a fixed small h is
  # pure round-off there; the wider step plus the higher order keeps both the
  # truncation and the cancellation below the density being measured.
  deriv <- function(tt, a) {
    h <- 1e-3 * tt
    (-P(tt + 2 * h, a) + 8 * P(tt + h, a) - 8 * P(tt - h, a) + P(tt - 2 * h, a)) / (12 * h)
  }
  for (a in c(1.1, 1.5, 3, 10, 40)) {
    times <- .bawdd_ref_times(A, b, k, ell, a)
    for (tt in c(0.3, 0.8, 2.2, times[["live"]], times[["partial"]]))
      expect_equal(D(tt, a), deriv(tt, a), tolerance = 1e-5)
    # Past T_max every trajectory has peaked, so the density is exactly zero.
    expect_identical(D(times[["frozen"]], a), 0)
  }
})

test_that("alpha is finite and well behaved across its whole range", {
  skip_on_cran()
  A <- 0.4; b <- 1.2; k <- 0.9; ell <- 0.25
  as <- 10^seq(-3, 4, length.out = 22)
  p <- vapply(as, function(a)
    EMC2:::pbawdd(1, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k, ell = ell, alpha = a), numeric(1))
  d <- vapply(as, function(a)
    EMC2:::dbawdd(1, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k, ell = ell, alpha = a), numeric(1))
  expect_true(all(is.finite(p)) && all(p > 0) && all(p < 1))
  expect_true(all(is.finite(d)) && all(d > 0))
  expect_true(all(diff(p) > 0))  # slower decay always helps
  # alpha = Inf is the no-decay limit the sequence converges to (as 1/alpha,
  # so take the limit far enough out to see it).
  expect_equal(EMC2:::pbawdd(1, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k,
                             ell = ell, alpha = 1e6),
               EMC2:::pbawdd(1, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k,
                             ell = ell, alpha = Inf),
               tolerance = 1e-4)
})

test_that("the simulator reaches threshold for alpha > 1 with ell > 0", {
  skip_on_cran()
  set.seed(20260827)
  A <- 0.4; B <- 0.8; b <- B + A; k <- 0.9; ell <- 0.25
  N <- 60000
  for (a in c(1.5, 3)) {
    pars <- cbind(mean = rep(0.3, N), cv = 0.5, B = B, A = A, t0 = 0, k = k,
                  ell = ell, alpha = a, eta = 0, b = b)
    out <- EMC2:::rbawdd_cpp(pars, "a", rep(TRUE, N), 1L, TRUE)
    hit <- is.finite(out$rt)
    p_inf <- EMC2:::pbawdd(Inf, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k, ell = ell, alpha = a)
    expect_gt(mean(hit), 0.5)          # not the all-omission behaviour it had
    # Binomial standard errors, so compare absolute deviations, not relative.
    expect_lt(abs(mean(hit) - p_inf), 4 * sqrt(p_inf * (1 - p_inf) / N))
    for (tt in c(0.6, 1.2)) {
      ana <- EMC2:::pbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k, ell = ell, alpha = a)
      expect_lt(abs(mean(hit & out$rt <= tt) - ana),
                4 * sqrt(ana * (1 - ana) / N))
    }
  }
})

test_that("Ttransform reports the finite endpoint that ell > 0 creates", {
  m <- EMC2::BAwDD()
  pars <- as.matrix(data.frame(
    mean = 1, cv = 0.5, B = 0.8, A = 0.4, t0 = 0.1, k = 0.9,
    ell = c(0.25, 0.25, 0, 0.25), alpha = c(0.5, 3, 3, Inf), eta = 0))
  out <- m$Ttransform(pars, NULL)
  expect_equal(unname(out[, "b"]), rep(1.2, 4))
  # ell > 0 and k > 0 give a finite endpoint on both sides of alpha = 1;
  # ell = 0 and the alpha = Inf (no-decay) member do not.
  expect_true(all(is.finite(out[1:2, "Tmax"])))
  expect_equal(unname(out[3:4, "Tmax"]), c(Inf, Inf))
  expect_equal(unname(out[, "rt_max"]), unname(out[, "t0"] + out[, "Tmax"]))
  # It must agree with the shared BAwD endpoint at the matching fixed rho.
  expect_equal(unname(out[2, "Tmax"]),
               EMC2:::bawd_tmax(A = 0.4, b = 1.2, k = 0.9, ell = 0.25,
                                gamma = 0, rho = 1 / 3))
})

test_that("the saturation-solve memo is transparent to evaluation order", {
  skip_on_cran()
  # bawd_geometry() memoises its Newton solve for the saturation boundary on
  # (c, rho, gamma).  A key that dropped any of those, or a set-associative
  # collision that returned the wrong way, would show up as a row picking up a
  # neighbouring row's geometry -- so the same rows must agree however they are
  # batched and whatever order they arrive in.
  set.seed(99)
  n <- 400L
  A <- 0.4; b <- 1.2; ell <- 0.25
  alpha <- sample(c(0, 0.3, 1, 1 + 1e-12, 2.5, 40, Inf), n, TRUE)
  k <- sample(c(0.2, 0.9, 5), n, TRUE)
  # Past T_sat_A, so the memoised saturation boundary really is in the answer:
  # below it the frozen branch contributes nothing and y_0/y_A go unread.
  tt <- stats::runif(n, 3, 40)
  batched <- EMC2:::dbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k,
                           ell = ell, alpha = alpha)
  one_at_a_time <- vapply(rev(seq_len(n)), function(i)
    EMC2:::dbawdd(tt[i], A = A, b = b, p1 = 0.3, p2 = 0.5, k = k[i],
                  ell = ell, alpha = alpha[i]), numeric(1))
  expect_identical(batched, rev(one_at_a_time))
  expect_identical(EMC2:::pbawdd(tt, A = A, b = b, p1 = 0.3, p2 = 0.5, k = k,
                                 ell = ell, alpha = alpha),
                   vapply(seq_len(n), function(i)
                     EMC2:::pbawdd(tt[i], A = A, b = b, p1 = 0.3, p2 = 0.5,
                                   k = k[i], ell = ell, alpha = alpha[i]),
                     numeric(1)))
})

test_that("more distinct geometries than the memo holds still match the reference", {
  skip_on_cran()
  # Far more keys than the table has slots, so most rows evict; the answers
  # must come from the closed form either way.
  A <- 0.4; b <- 1.2; k <- 0.9; ell <- 0.25
  Sbar <- function(v) stats::pnorm((log(v) - 0.3) / 0.5, lower.tail = FALSE)
  alpha <- seq(0.4, 6, length.out = 300)
  u <- 6  # past T_sat_A for every alpha on this grid
  got <- EMC2:::pbawdd(rep(u, length(alpha)), A = A, b = b, p1 = 0.3,
                       p2 = 0.5, k = k, ell = ell, alpha = alpha)
  expect_true(all(is.finite(got)) && all(diff(got) > 0))
  for (i in c(1L, 87L, 150L, 251L, 300L))
    expect_equal(got[i], .bawdd_ref_cdf(u, A, b, k, ell, alpha[i], Sbar),
                 tolerance = 1e-8)
})
