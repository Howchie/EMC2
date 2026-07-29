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

    expect_lt(max(abs(r$pdf - wald_pdf(tg, mu / sigma, b0 / sigma))), 3e-3)
    expect_lt(max(abs(r$cdf - wald_cdf(tg, mu / sigma, b0 / sigma))), 3e-4)
  }
})

test_that("OU first passage matches the closed form when b0 == theta", {
  # NB this case nearly zeroes the Volterra kernel and so tests little there;
  # for the Fokker-Planck solver it does not degenerate and is a real check.
  for (lambda in c(0.5, 2, 5)) {
    r <- EMC2:::fpe_ou_fht_pdf_cdf_vec(tg, lambda, 1, 1, 0, 1, 1, 1, 1, 512L, 1024L)

    expect_lt(max(abs(r$pdf - EMC2:::ou_fht_pdf_vec_closed_form(
      tg, lambda, 1, 1, 0, 1, 1, 1, 1))), 3e-3)
    expect_lt(max(abs(r$cdf - EMC2:::ou_fht_cdf_vec_closed_form(
      tg, lambda, 1, 1, 0, 1, 1, 1, 1))), 2e-4)
  }
})

test_that("OU with b0 != theta agrees with the Volterra solver", {
  # The two disagree by ~0.3% at the pdf peak: the Fokker-Planck solution is
  # self-converged across nx there and the Volterra one is not, so this is a
  # coarse cross-check rather than a tight one.
  for (lambda in c(0.5, 2)) {
    r <- EMC2:::fpe_ou_fht_pdf_cdf_vec(tg, lambda, 2, 1, 0, 1, 1, 1, 1, 512L, 1024L)
    expect_lt(max(abs(r$pdf - EMC2:::ou_fht_pdf_vec(tg, lambda, 2, 1, 0, 1, 1, 1, 1))),
              2e-2)
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
  expect_lt(max(abs(uniform - mixture)), 5e-3)
})

test_that("GBM reduces to the Wald form in log space for a point start", {
  # X starts at start_floor and is absorbed at b0, so Y = log X is Brownian with
  # drift mu - sigma^2/2 travelling log(b0) - log(start_floor).
  mu <- 0.5; sigma <- 0.5; b0 <- 2; z <- 1
  r <- EMC2:::fpe_gbm_fht_pdf_cdf_vec(tg, mu, sigma, 0, b0, b0, 1, 1, 512L, 1024L, z)

  drift <- (mu - 0.5 * sigma^2) / sigma
  dist  <- (log(b0) - log(z)) / sigma
  expect_lt(max(abs(r$pdf - wald_pdf(tg, drift, dist))), 3e-3)
  expect_lt(max(abs(r$cdf - wald_cdf(tg, drift, dist))), 3e-4)
})
