test_that("rejected sample iteration repeats previous stored draw", {
  samples <- list(
    theta_mu = matrix(1:6, nrow = 2),
    theta_var = array(1:12, dim = c(2, 2, 3)),
    last_theta_var_inv = diag(2),
    stage = c("init", "burn", "burn"),
    idx = 2
  )

  out <- reject_sample_iteration(samples, 3)

  expect_equal(out$theta_mu[, 3], samples$theta_mu[, 2])
  expect_equal(out$theta_var[, , 3], samples$theta_var[, , 2])
  expect_equal(out$last_theta_var_inv, samples$last_theta_var_inv)
  expect_equal(out$idx, 3)
})

test_that("failed particle proposal is rejected to previous subject state", {
  pm_settings <- list(list(epsilon = 1))
  parameters <- list(alpha = matrix(c(1, 2), ncol = 1))
  out <- safe_new_particle(
    s = 1,
    data = NULL,
    pm_settings = pm_settings,
    prev_ll = -7,
    parameters = parameters,
    model = NULL,
    stage = "preburn",
    type = "standard",
    tune = list(components = c(1, 1), shared_ll_idx = c(1, 1)),
    r_cores = 1
  )

  expect_equal(out$proposal, parameters$alpha[, 1])
  expect_equal(out$ll, -7)
  expect_equal(out$pm_settings, pm_settings)
})

test_that("direct worker state rejects to its local current alpha", {
  pm_settings <- list(list(epsilon = 1))
  out <- safe_new_particle(
    s = 9,
    data = NULL,
    pm_settings = pm_settings,
    prev_ll = -11,
    parameters = NULL,
    current_alpha = c(a = 4, b = 5),
    population_mu = c(a = 0, b = 0),
    population_var = diag(2),
    model = NULL,
    stage = "preburn",
    type = "standard",
    tune = list(components = c(1, 1), shared_ll_idx = c(1, 1)),
    r_cores = 1
  )
  expect_equal(out$proposal, c(a = 4, b = 5))
  expect_equal(out$ll, -11)
})
