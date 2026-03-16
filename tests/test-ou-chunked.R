test_that("chunked OU solver matches single Volterra at long horizons", {
  skip_on_cran()

  lambda <- 2.81674853375414
  theta <- 0.783897707122378
  sigma <- 0.51536637947429
  z0 <- 0.333658571015112
  tau <- 4
  pow <- 3
  b0 <- theta
  binf <- b0

  steps_fineness <- 5e-4
  min_steps <- 0
  min_t_fineness <- 1e-3
  chunk_ratio <- 2
  chunk_base_panels <- 200
  chunk_max <- 15

  t_grid <- seq(0, 10, by = 0.01)
  target_idx <- which.min(abs(t_grid - 0.75))
  t_eval <- t_grid[target_idx]

  pdf_chunked <- EMC2:::ou_fht_pdf_vec_grid_chunked(
    t_grid, lambda, theta, sigma, z0, b0, binf,
    tau, pow,
    steps_fineness, min_steps, min_t_fineness,
    chunk_ratio, chunk_base_panels, chunk_max
  )

  pdf_uniform <- EMC2:::ou_fht_pdf_vec(
    t_eval, lambda, theta, sigma, z0, b0, binf,
    tau, pow,
    steps_fineness, min_steps
  )

  pdf_analytic <- EMC2:::ou_fht_pdf_vec_closed_form(
    t_eval, lambda, theta, sigma, z0, b0, binf
  )

  expect_equal(pdf_chunked[target_idx], pdf_uniform[1], tolerance = 1e-6)
  expect_equal(pdf_chunked[target_idx], pdf_analytic[1], tolerance = 1e-6)
  expect_true(pdf_chunked[target_idx] > 0)

  cdf_chunked <- EMC2:::ou_fht_cdf_vec_grid_chunked(
    t_grid, lambda, theta, sigma, z0, b0, binf,
    tau, pow,
    steps_fineness, min_steps, min_t_fineness,
    chunk_ratio, chunk_base_panels, chunk_max,
    rt_resolution = 0.01
  )

  cdf_uniform <- EMC2:::ou_fht_cdf_vec(
    t_eval, lambda, theta, sigma, z0, b0, binf,
    tau, pow,
    num_steps = 4000
  )

  expect_equal(cdf_chunked[target_idx], cdf_uniform[1], tolerance = 1e-6)

  trap_cdf <- numeric(length(t_grid))
  for (i in seq_len(length(t_grid))[-1]) {
    dt <- t_grid[i] - t_grid[i - 1]
    trap_cdf[i] <- trap_cdf[i - 1] + 0.5 * dt * (pdf_chunked[i] + pdf_chunked[i - 1])
  }

  expect_equal(cdf_chunked[target_idx], trap_cdf[target_idx], tolerance = 5e-4)
})

test_that("chunked OU solver dataset matches recorded diffs", {
  skip_on_cran()

  cases <- read.csv(test_path("ou_chunked_cases.csv"))
  expected <- read.csv(test_path("ou_chunked_eval_results.csv"))

  steps_fineness <- 1e-3
  min_steps <- 0
  min_t_fineness <- 2e-3
  chunk_ratio <- 1.5
  chunk_base_panels <- 80
  chunk_max <- 12
  tau <- 4
  pow <- 3
  rt_resolution <- 0.02

  t_grid <- seq(0, 10, by = 0.01)

  pdf_chunked <- numeric(nrow(cases))
  pdf_uniform <- numeric(nrow(cases))
  cdf_chunked <- numeric(nrow(cases))
  cdf_uniform <- numeric(nrow(cases))
  idx <- integer(nrow(cases))
  t_eval <- numeric(nrow(cases))

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    idx[i] <- max(1L, min(length(t_grid), as.integer(floor(case$t))))
    t_eval[i] <- t_grid[idx[i]]

    pdf_chunked_vec <- EMC2:::ou_fht_pdf_vec_grid_chunked(
      t_grid,
      case$lambda, case$theta, case$sigma, case$z0,
      case$theta, case$theta,
      tau, pow,
      steps_fineness, min_steps, min_t_fineness,
      chunk_ratio, chunk_base_panels, chunk_max
    )
    pdf_uniform_vec <- EMC2:::ou_fht_pdf_vec(
      t_eval[i],
      case$lambda, case$theta, case$sigma, case$z0,
      case$theta, case$theta,
      tau, pow,
      steps_fineness, min_steps
    )

    cdf_chunked_vec <- EMC2:::ou_fht_cdf_vec_grid_chunked(
      t_grid,
      case$lambda, case$theta, case$sigma, case$z0,
      case$theta, case$theta,
      tau, pow,
      steps_fineness, min_steps, min_t_fineness,
      chunk_ratio, chunk_base_panels, chunk_max,
      rt_resolution
    )
    cdf_uniform_vec <- EMC2:::ou_fht_cdf_vec(
      t_eval[i],
      case$lambda, case$theta, case$sigma, case$z0,
      case$theta, case$theta,
      tau, pow,
      num_steps = 2000
    )

    pdf_chunked[i] <- pdf_chunked_vec[idx[i]]
    pdf_uniform[i] <- pdf_uniform_vec[1]
    cdf_chunked[i] <- cdf_chunked_vec[idx[i]]
    cdf_uniform[i] <- cdf_uniform_vec[1]
  }

  computed <- data.frame(
    row = seq_len(nrow(cases)),
    idx = idx,
    t_eval = t_eval,
    pdf_chunked = pdf_chunked,
    pdf_uniform = pdf_uniform,
    pdf_diff = pdf_chunked - pdf_uniform,
    cdf_chunked = cdf_chunked,
    cdf_uniform = cdf_uniform,
    cdf_diff = cdf_chunked - cdf_uniform
  )

  expect_equal(computed$pdf_diff, expected$pdf_diff, tolerance = 1e-10)
  expect_equal(computed$cdf_diff, expected$cdf_diff, tolerance = 1e-10)

  expect_equal(sum(abs(computed$pdf_diff) > 1e-6), 15L)
  expect_equal(sum(abs(computed$pdf_diff) > 1e-5), 1L)
  expect_equal(sum(abs(computed$pdf_diff) > 1e-4), 1L)
  expect_equal(sum(abs(computed$cdf_diff) > 1e-6), 0L)
})
