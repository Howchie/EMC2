# Race Lévy Flight nonlocal Fokker-Planck solver and CMS simulator.

tg <- seq(0.1, 2.0, by = 0.05)

wald_pdf <- function(t, mu, b) {
  b / sqrt(2 * pi * t^3) * exp(-(b - mu * t)^2 / (2 * t))
}
wald_cdf <- function(t, mu, b) {
  pnorm((mu * t - b) / sqrt(t)) +
    exp(2 * mu * b) * pnorm(-(b + mu * t) / sqrt(t))
}

test_that("RLF PDE solver reduces continuously to Brownian first passage", {
  for (cs in list(c(1, 1, 1), c(2, 1, 1), c(0.5, 1, 1.5))) {
    v <- cs[1]
    sigma <- cs[2]
    b0 <- cs[3]
    r <- EMC2:::rlf_fht_pdf_cdf_vec(
      tg, v, sigma, 2.0, b0, 0, 250L, 500L
    )

    expect_lt(
      max(abs(r$pdf - wald_pdf(tg, v / sigma, b0 / sigma))),
      0.18
    )
    expect_lt(
      max(abs(r$cdf - wald_cdf(tg, v / sigma, b0 / sigma))),
      0.022
    )
    expect_lt(r$mismatch, 1e-10)
  }

  near <- EMC2:::rlf_fht_pdf_cdf_vec(
    tg, 1, 1, 1.9999, 1, 0, 250L, 500L
  )
  limit <- EMC2:::rlf_fht_pdf_cdf_vec(
    tg, 1, 1, 2, 1, 0, 250L, 500L
  )
  expect_lt(max(abs(near$pdf - limit$pdf)), 5e-4)
  expect_lt(max(abs(near$cdf - limit$cdf)), 1e-4)
})

test_that("RLF flux, mass, and conservative lower closure are consistent", {
  r <- EMC2:::rlf_fht_pdf_cdf_vec(
    tg, 1, 1, 1.3, 1, 0, 250L, 500L
  )

  expect_true(all(diff(r$cdf) >= -1e-10))
  expect_true(all(r$cdf >= 0 & r$cdf <= 1))
  expect_true(all(r$pdf >= 0))
  expect_gte(r$min_density, -1e-10)
  expect_lt(r$mismatch, 1e-10)
  expect_lt(r$operator_conservation_error, 1e-10)
  expect_gt(r$lower_boundary_pressure, 0)
  expect_lt(r$lower_boundary_pressure, 0.01)
  expect_true(r$refinement_checked)
  expect_false(r$refinement_skipped)
})

test_that("RLF uses the cheap path for safely resolved routine cases", {
  r <- EMC2:::rlf_fht_pdf_cdf_vec(
    tg, 1, 1, 1.7, 1, 0, 250L, 500L
  )

  expect_false(r$refinement_checked)
  expect_true(r$refinement_skipped)
  expect_equal(r$domain_refinements, 0L)
  expect_equal(r$spatial_refinements, 0L)
  expect_true(is.nan(r$domain_cdf_error))
  expect_true(is.nan(r$spatial_cdf_error))
  expect_true(all(diff(r$cdf) >= -1e-10))
  expect_true(all(r$pdf >= 0))
  expect_lt(r$lower_boundary_pressure, 0.002)
})

test_that("RLF automatically expands and refines space and time", {
  r <- EMC2:::rlf_fht_pdf_cdf_vec(
    seq(0.01, 2, length.out = 201),
    1, 1, 1.7, 1, 0, 80L, 100L
  )

  expect_gt(r$domain_refinements, 0L)
  expect_gt(r$spatial_refinements, 0L)
  expect_lt(r$x_lo, -12)
  expect_gt(r$nx_used, 80L)
  expect_gt(r$nt_used, 100L)
  expect_lte(r$domain_pdf_error, 0.01)
  expect_lte(r$domain_cdf_error, 0.002)
  expect_lte(r$spatial_pdf_error, 0.05)
  expect_lte(r$spatial_cdf_error, 0.0125)
  expect_true(all(diff(r$cdf) >= -1e-10))
  expect_true(all(r$cdf >= 0 & r$cdf <= 1 + 1e-12))
  expect_true(all(r$pdf >= 0))
  expect_gte(r$min_density, -1e-10)
  expect_lt(r$mismatch, 1e-10)
})

test_that("RLF rejects spatial regimes that exceed its refinement budget", {
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(
      seq(0.001, 0.3, length.out = 101),
      10, 0.05, 1.99, 1, 0, 120L, 120L
    ),
    "automatic spatial refinement failed to converge"
  )
})

test_that("RLF rejects unsupported drift and parameter ranges", {
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, 0, 1, 1.7, 1),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, -1, 1, 1.7, 1),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, NA_real_, 1, 1.7, 1),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, Inf, 1, 1.7, 1),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, 1, 1, 1, 1),
    "alpha must be finite and in"
  )
  expect_error(
    EMC2:::rlf_fht_pdf_cdf_vec(1, 1, 1, 1.7, 1, 1),
    "z0 must be finite and in"
  )

  expect_error(
    EMC2:::simulate_rlf_hit_times_cpp(
      10L, 0, 1, 1.7, 1, 0, 1, 0.01, 1L
    ),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::simulate_rlf_hit_times_cpp(
      10L, NA_real_, 1, 1.7, 1, 0, 1, 0.01, 1L
    ),
    "v must be finite and positive"
  )
  expect_error(
    EMC2:::simulate_rlf_hit_times_cpp(
      10L, 1, 1, 1, 1, 0, 1, 0.01, 1L
    ),
    "alpha must be finite and in"
  )
})

test_that("CMS simulator uses the same Brownian scale as the PDE", {
  hits <- EMC2:::simulate_rlf_hit_times_cpp(
    30000L, 1, 1, 2, 1, 0, 2, 0.01, 42L
  )
  empirical <- mean(!is.na(hits) & hits <= 2)
  expect_equal(empirical, wald_cdf(2, 1, 1), tolerance = 0.008)
})

test_that("RLF simulator stops exactly at a partial terminal step", {
  hits <- EMC2:::simulate_rlf_hit_times_cpp(
    20L, 1, 1e-9, 2, 0.104, 0, 0.105, 0.02, 42L
  )

  expect_true(all(hits <= 0.105))
  expect_equal(hits, rep(0.105, length(hits)), tolerance = 1e-14)
})

test_that("heavy-tailed CMS simulation agrees with the nonlocal PDE", {
  hits <- EMC2:::simulate_rlf_hit_times_cpp(
    20000L, 1, 1, 1.7, 1, 0, 2, 0.001, 42L
  )
  empirical <- mean(!is.na(hits) & hits <= 2)
  pde <- EMC2:::rlf_fht_pdf_cdf_vec(
    2, 1, 1, 1.7, 1, 0, 250L, 500L
  )

  expect_equal(empirical, pde$cdf[1], tolerance = 0.015)
  expect_lt(pde$lower_boundary_pressure, 0.005)
})
