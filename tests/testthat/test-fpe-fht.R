# Fokker-Planck first-passage solver (src/fpe_solver.h, src/fpe_models.h).
#
# Deterministic checks only -- the Monte Carlo arms of the validation ladder live
# in .Rtmp/fpe_validate.R, since they are slow and non-deterministic.  Tolerances
# are loose relative to the measured errors so that ordinary retuning of the
# defaults does not trip them, but tight enough to catch a real regression.

tg <- seq(0.05, 2.0, by = 0.025)

wald_pdf <- function(t, mu, b) b / sqrt(2 * pi * t^3) * exp(-(b - mu * t)^2 / (2 * t))
wald_cdf <- function(t, mu, b) {
  pnorm((mu * t - b) / sqrt(t)) + exp(2 * mu * b) * pnorm(-(b + mu * t) / sqrt(t))
}

test_that("BM first passage matches the Wald closed form", {
  for (cs in list(c(1, 1, 1), c(2, 1, 1), c(0.5, 1, 1.5))) {
    mu <- cs[1]; sigma <- cs[2]; b0 <- cs[3]
    r <- EMC2:::fpe_bm_fht_pdf_cdf_vec(tg, mu, sigma, 0, b0, b0, 1, 1, 512L, 1024L)

    expect_lt(max(abs(r$pdf - wald_pdf(tg, mu / sigma, b0 / sigma))), 2e-3)
    expect_lt(max(abs(r$cdf - wald_cdf(tg, mu / sigma, b0 / sigma))), 1e-4)
  }
})


test_that("the two routes to the cdf agree once the solve is resolved", {
  # CDF from lost mass (1 - h*sum(q)) versus CDF from the integrated boundary
  # flux.  Agreement is not automatic -- the mass identity alone is a tautology
  # of the discretisation -- so this is a genuine per-solve error estimate.
  expect_lt(EMC2:::fpe_bm_fht_pdf_cdf_vec(tg, 1, 1, 0, 1, 1, 1, 1, 512L, 1024L)$mismatch,
            1e-6)
  expect_lt(EMC2:::fpe_ou_fht_pdf_cdf_vec(tg, 2, 2, 1, 0, 1, 0.5, 1, 1, 512L, 1024L)$mismatch,
            1e-6)
  expect_lt(EMC2:::fpe_gompertz_fht_pdf_cdf_vec(tg, 2, 0.5, 0.5, 2, 2, 1, 1,
                                                512L, 1024L, 1e-3)$mismatch, 1e-6)
})

test_that("the scheme is second order in (h, dt) jointly", {
  ref <- wald_cdf(tg, 1, 1)
  err <- vapply(c(1, 2, 4), function(k) {
    r <- EMC2:::fpe_bm_fht_pdf_cdf_vec(tg, 1, 1, 0, 1, 1, 1, 1,
                                       as.integer(128 * k), as.integer(256 * k))
    max(abs(r$cdf - ref))
  }, numeric(1))

  # halving h and dt should cut the error by ~4
  expect_gt(err[1] / err[2], 3.0)
  expect_gt(err[2] / err[3], 3.0)
})

test_that("uniform start-point variability is uniform in the physical state", {
  # Regression for a solver bug: GBM and Gompertz integrate Y = log X, and
  # seeding the start uniform in Y rather than in X tilts the initial density by
  # (Zhi/Zlo).  That is nearly invisible over a narrow start range and enormous
  # over a wide one, so test Gompertz, whose range here spans 500x.
  #
  # Oracle: the uniform-start solution must equal the mixture of point-start
  # solutions over the physical start distribution.  A point start at z is
  # obtained by setting both z0 and start_floor to z.
  zs <- seq(1e-3, 0.5, length.out = 40)
  zs <- zs[-1] - diff(zs)[1] / 2                      # midpoint rule
  mixture <- rowMeans(vapply(zs, function(z) {
    EMC2:::fpe_gompertz_fht_pdf_cdf_vec(tg, 2, 0.5, z, 2, 2, 1, 1, 256L, 512L, z)$cdf
  }, numeric(length(tg))))

  uniform <- EMC2:::fpe_gompertz_fht_pdf_cdf_vec(tg, 2, 0.5, 0.5, 2, 2, 1, 1,
                                                 256L, 512L, 1e-3)$cdf
  # The residual here is the point-start SEED error carried by the 39 members of
  # the oracle, not quadrature error -- it does not shrink as the mixture is
  # refined.  It fell 4x when the mesh was graded and t_seed capped.
  expect_lt(max(abs(uniform - mixture)), 1.5e-3)
})

test_that("GBM reduces to the Wald form in log space for a point start", {
  # X starts at start_floor and is absorbed at b0, so Y = log X is Brownian with
  # drift mu - sigma^2/2 travelling log(b0) - log(start_floor).
  mu <- 0.5; sigma <- 0.5; b0 <- 2; z <- 1
  r <- EMC2:::fpe_gbm_fht_pdf_cdf_vec(tg, mu, sigma, 0, b0, b0, 1, 1, 512L, 1024L, z)

  drift <- (mu - 0.5 * sigma^2) / sigma
  dist  <- (log(b0) - log(z)) / sigma
  expect_lt(max(abs(r$pdf - wald_pdf(tg, drift, dist))), 3e-4)
  expect_lt(max(abs(r$cdf - wald_cdf(tg, drift, dist))), 1e-4)
})

test_that("the graded mesh is a strict improvement and grade = 1 is uniform", {
  # Cells are uniform in a stretched coordinate and clustered at the absorbing
  # barrier, where the solution actually varies; the far field holds almost no
  # mass.  `grade` is the ratio of far-field to barrier cell width.
  ref <- wald_cdf(tg, 1, 1)
  err <- vapply(c(1, 2, 4, 8), function(gr) {
    # nt is deliberately generous: at nt = 512 the time error masks the gain.
    r <- EMC2:::fpe_bm_fht_pdf_cdf_vec(tg, 1, 1, 0, 1, 1, 1, 1, 256L, 2048L, gr)
    max(abs(r$cdf - ref))
  }, numeric(1))

  expect_true(all(diff(err) < 0))          # monotone in the grading
  expect_gt(err[1] / err[4], 4.0)          # ~5x at the default grade of 8

  # grade = 1 must reproduce the uniform mesh bit for bit, since the whole
  # geometry collapses back to a constant h.
  expect_equal(err[1], 2.352e-04, tolerance = 1e-3)
})

test_that("grading does not disturb the order of accuracy", {
  ref <- wald_cdf(tg, 1, 1)
  err <- vapply(c(1, 2, 4), function(k) {
    r <- EMC2:::fpe_bm_fht_pdf_cdf_vec(tg, 1, 1, 0, 1, 1, 1, 1,
                                       as.integer(256 * k), as.integer(2048 * k), 8)
    max(abs(r$cdf - ref))
  }, numeric(1))
  expect_gt(err[1] / err[2], 3.0)
  expect_gt(err[2] / err[3], 3.0)
})

# ---------------------------------------------------------------------------
# Fixes for the four defects found while evaluating the solver as a likelihood
# backend (see race_ou_integration_plan.md §2.3).
# ---------------------------------------------------------------------------

test_that("graded time stepping resolves the rising flank at realistic t_max", {
  # A uniform dt = t_max/nt spends its steps in proportion to elapsed time, but
  # the density climbs from ~0 to its mode inside the first ~100 ms.  With t_max
  # set by the largest RT in a data set the flank is then badly under-resolved,
  # and the error there is a function of dt alone -- it is identical at
  # t_max = 0.25 and t_max = 2.5.  `tgrade` doubles dt between blocks, so the
  # early steps are small, at no extra cost (same total step budget).
  ou <- function(t, tgrade, nx = 256L, nt = 512L)
    EMC2:::fpe_ou_fht_pdf_cdf_vec(t, 4, 2, 1, 0, 1, 1, 1, 1, nx, nt, 8, tgrade)

  tq  <- c(0.05, 0.08, 0.12, 0.5, 1.0)
  ref <- ou(tq, 1, 2048L, 16384L)$pdf
  rel <- function(tgrade) {
    r <- ou(c(tq, 2.5), tgrade)$pdf[seq_along(tq)]
    abs(r - ref) / ref
  }

  # the flank at t = 0.05 is where it bites: ~26% uniform, ~1% graded
  expect_gt(rel(1)[1], 0.15)
  expect_lt(rel(32)[1], 0.03)
  # and the mode / body must not be traded away to get it
  expect_lt(max(rel(32)[2:4]), 5e-3)

  # tgrade <= 1 must reproduce the uniform grid exactly, so that every test
  # above this line keeps testing what it says it tests.
  expect_identical(ou(tq, 1)$pdf, ou(tq, 0.5)$pdf)
})

test_that("the reported pdf and cdf are safe to take logs of", {
  # The absorbing-face flux is a one-sided difference of two nearly equal cell
  # values, so once the sub-density has decayed it goes slightly negative
  # (-3e-16 measured where the true density is 1.5e-17).  log() of that is NaN,
  # which would silently poison a likelihood rather than fail loudly.  Likewise
  # 1 - sum(dx*q) is not reliably inside [0,1] nor monotone out there.
  t <- c(seq(1e-4, 0.01, length.out = 20), seq(0.05, 3, length.out = 200))
  r <- EMC2:::fpe_ou_fht_pdf_cdf_vec(t, 4, 2, 1, 0, 1, 1, 1, 1, 256L, 512L, 8, 32)

  expect_true(all(r$pdf >= 0))
  expect_true(all(r$pdf_grid >= 0))
  expect_true(all(r$cdf >= 0 & r$cdf <= 1))
  expect_true(all(diff(r$cdf) >= 0))
  expect_false(anyNA(log(pmax(r$pdf, 1e-300))))
})

test_that("the pdf is zero before the seed time, not clamped to it", {
  # A point start seeds the march at t_seed > 0, so the grid does not reach back
  # to 0.  Clamping the lookup to the first grid value hands back the flux AT
  # t_seed for every earlier time -- a positive constant where the density is
  # ~0.  The cdf keeps the clamp: cdf[1] is mass genuinely absorbed by t_seed.
  # t_seed is capped at 0.25 * t_max, so the query has to carry a realistic
  # t_max for the seed to sit where it would in a fit.
  r <- EMC2:::fpe_ou_fht_pdf_cdf_vec(c(1e-6, 1e-4, 2.5), 4, 2, 1, 0, 1, 1, 1, 1,
                                     256L, 512L, 8, 1)
  expect_gt(r$t_grid[1], 1e-4)          # first two query times precede the seed
  expect_identical(r$pdf[1:2], c(0, 0))
  expect_gt(r$pdf[3], 0)
})
