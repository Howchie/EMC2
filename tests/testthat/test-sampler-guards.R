test_that("inter-block proposal scaling converts variance to standard deviation", {
  pmwgs <- list(
    samples = structure(
      list(),
      pm_settings = list(list(list(epsilon = c(1, 2))))
    )
  )

  out <- update_epsilon_scale(pmwgs, prop_var_ratio = 4)
  expect_equal(attr(out$samples, "pm_settings")[[1]][[1]]$epsilon,
               c(2, 4))

  bounded <- update_epsilon_scale(pmwgs, prop_var_ratio = 1e100)
  expect_true(all(is.finite(attr(bounded$samples, "pm_settings")[[1]][[1]]$epsilon)))
  expect_lte(max(attr(bounded$samples, "pm_settings")[[1]][[1]]$epsilon), 5)
})

test_that("empirical covariance fallback resets proposal multipliers", {
  expect_true(.emc_empirical_covariance_ok(diag(2)))
  expect_false(.emc_empirical_covariance_ok(diag(c(1, 1e-20)), subject = 7))

  pmwgs <- list(
    samples = structure(
      list(),
      pm_settings = list(list(list(epsilon = c(1e38, 1e-20))))
    )
  )
  out <- update_epsilon_scale(pmwgs, reset = TRUE)
  expect_equal(attr(out$samples, "pm_settings")[[1]][[1]]$epsilon,
               c(1, 1))
})

test_that("each burn block resets particle adaptation counters", {
  samples <- structure(
    list(stage = "burn"),
    pm_settings = list(list(list(
      iter = 99L, proposal_counts = c(10, 20), acc_counts = c(3, 4),
      mix = c(.2, .8), epsilon = c(1, 1)
    )))
  )
  emc <- list(list(samples = samples))
  out <- reset_pm_settings(emc, "burn")
  settings <- attr(out[[1]]$samples, "pm_settings")[[1]][[1]]
  expect_equal(settings$iter, 25L)
  expect_null(settings$mix)
  expect_equal(settings$proposal_counts, c(0, 0))
  expect_equal(settings$acc_counts, c(0, 0))
})

test_that("ill-conditioned inverse-Wishart scales report diagnostics", {
  expect_error(
    riwish(5, diag(c(1, 1e-20)), context = "standard Gibbs block 7"),
    "standard Gibbs block 7.*rcond"
  )
})

test_that("group Gibbs failures are fatal and carry location", {
  cond <- tryCatch(
    .emc_gibbs_abort(
      structure(list(message = "singular"), class = c("error", "condition")),
      stage = "burn", iteration = 123
    ),
    error = identity
  )
  expect_s3_class(cond, "emc_gibbs_failure")
  expect_match(conditionMessage(cond), "stage=burn iteration=123")
})

test_that("long runs of identical sample iterations are detectable", {
  samples <- list(
    theta_mu = matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2),
    alpha = array(c(1, 2, 3, 3), dim = c(2, 1, 2)),
    stage = rep("burn", 2),
    idx = 2
  )
  expect_false(.emc_sample_iteration_equal(samples, 2))
  samples$theta_mu[, 2] <- samples$theta_mu[, 1]
  samples$alpha[, , 2] <- samples$alpha[, , 1]
  expect_true(.emc_sample_iteration_equal(samples, 2))
})
