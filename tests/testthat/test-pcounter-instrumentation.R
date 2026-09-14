test_that("PCOUNTER instrumentation is off by default and resettable", {
  old <- EMC2:::emc_kernel_stats()
  on.exit({
    EMC2:::emc_kernel_stats(old)
    EMC2:::emc_kernel_stats_reset()
  }, add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(FALSE)

  invisible(EMC2:::dpcounter(
    c(.2, .5), c(8, 8), c(0, 1), c(0, .5), c(0, 2), c(0, 1),
    c(0, 0), log_out = TRUE
  ))
  expect_identical(nrow(EMC2:::emc_kernel_stats_read()), 0L)

  EMC2:::emc_kernel_stats(TRUE)
  invisible(EMC2:::dpcounter(
    c(.2, .5), c(8, 8), c(0, 1), c(0, .5), c(0, 2), c(0, 1),
    c(0, 0), log_out = TRUE
  ))
  st <- EMC2:::emc_kernel_stats_read()
  pc <- st[st$model == "PCOUNTER", , drop = FALSE]
  expect_equal(nrow(pc), 1L)
  expect_gt(pc$k_observations, 0)
  numeric_pc <- pc[vapply(pc, is.numeric, logical(1))]
  expect_true(all(vapply(numeric_pc, function(x) all(is.finite(x)), logical(1))))

  EMC2:::emc_kernel_stats_reset()
  expect_identical(nrow(EMC2:::emc_kernel_stats_read()), 0L)
})

test_that("PCOUNTER instrumentation records branches, support exits, and K terms", {
  old <- EMC2:::emc_kernel_stats()
  on.exit({
    EMC2:::emc_kernel_stats(old)
    EMC2:::emc_kernel_stats_reset()
  }, add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(TRUE)

  # The first three rows cover all degenerate switches; the final row forces a
  # large finite threshold and exercises the bounded tail machinery when it is
  # numerically needed.  The two invalid rows distinguish parameter and K
  # support exits without changing the successful rows.
  invisible(EMC2:::dpcounter(
    c(.2, .4, .7, .2, .2, .2),
    c(8, 8, 8, 8, 8, 8),
    c(0, 0, 1, 1, 1, 1),
    c(0, .5, 0, .5, .5, .5),
    c(0, 2, 2, 1023, 0, 1e100),
    c(0, 0, 0, 1, 1, 1),
    rep(0, 6), log_out = TRUE
  ))
  invisible(EMC2:::dpcounter(.2, -1, 0, 0, 0, 0, 0, log_out = TRUE))
  invisible(EMC2:::dpcounter(.2, 8, 0, 0, -1, 0, 0, log_out = TRUE))
  invisible(EMC2:::dpcounter(.2, 8, 0, 0, NaN, 0, 0, log_out = TRUE))
  invisible(EMC2:::dpcounter(.2, 8, 0, 0, Inf, 0, 0, log_out = TRUE))

  st <- EMC2:::emc_kernel_stats_read()
  pc <- st[st$model == "PCOUNTER", , drop = FALSE]
  expect_equal(nrow(pc), 1L)
  expect_gte(pc$sv_zero, 1)
  expect_gte(pc$gamma_zero, 1)
  numeric_pc <- pc[vapply(pc, is.numeric, logical(1))]
  expect_true(all(vapply(numeric_pc, function(x) all(is.finite(x)), logical(1))))
  expect_gte(pc$support_exits, 1)
  expect_gte(pc$invalid_exits, 3)
  expect_gt(pc$stirling_terms + pc$rising_terms, 0)
  expect_gt(pc$fallback_tail, 0)
})

test_that("PCOUNTER wrappers count nonfinite launch times as invalid", {
  old <- EMC2:::emc_kernel_stats()
  on.exit({
    EMC2:::emc_kernel_stats(old)
    EMC2:::emc_kernel_stats_reset()
  }, add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(TRUE)

  invisible(EMC2:::dpcounter(.2, 8, 0, 0, 0, 0, NaN, log_out = TRUE))
  st <- EMC2:::emc_kernel_stats_read()
  pc <- st[st$model == "PCOUNTER", , drop = FALSE]
  expect_equal(pc$invalid_exits, 1)

  invisible(EMC2:::ppcounter(.2, 8, 0, 0, 0, 0, NaN,
                             lower_tail = TRUE, log_out = TRUE))
  st <- EMC2:::emc_kernel_stats_read()
  pc <- st[st$model == "PCOUNTER", , drop = FALSE]
  expect_equal(pc$invalid_exits, 2)
})

test_that("PCOUNTER raw adapters report calls, rows, and timing", {
  old <- EMC2:::emc_kernel_stats()
  on.exit({
    EMC2:::emc_kernel_stats(old)
    EMC2:::emc_kernel_stats_reset()
  }, add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(TRUE)

  fx <- audit_fixture("PCOUNTER", n_trials = 64L, n_particles = 2L,
                      branch = "truncation")
  ll <- audit_ll_direct(fx)
  expect_true(all(is.finite(ll)))
  st <- EMC2:::emc_kernel_stats_read()
  pc <- st[st$model == "PCOUNTER", , drop = FALSE]
  expect_equal(nrow(pc), 1L)
  expect_gt(pc$dpcounter_calls, 0)
  expect_gt(pc$ppcounter_calls, 0)
  expect_gt(pc$pcounter_logS_calls, 0)
  expect_gt(pc$dpcounter_rows + pc$ppcounter_rows + pc$pcounter_logS_rows, 0)
  expect_true(all(c(pc$dpcounter_seconds, pc$ppcounter_seconds,
                    pc$pcounter_logS_seconds, pc$preparation_seconds,
                    pc$rt_sum_seconds) >= 0))
})

test_that("PCOUNTER instrumentation leaves likelihood values bit-identical", {
  old <- EMC2:::emc_kernel_stats()
  on.exit({
    EMC2:::emc_kernel_stats(old)
    EMC2:::emc_kernel_stats_reset()
  }, add = TRUE)
  args <- list(
    t = c(.05, .2, .5, 1, 2),
    nu = c(3, 8, 12, 5, 20),
    sv = c(0, 1, 2, 0, 1),
    gamma = c(0, .5, .5, 0, .2),
    k = c(0, 1, 2, 5, 1023),
    omega = c(0, 1, 0, 1, 2),
    t0 = rep(0, 5),
    log_out = TRUE
  )
  EMC2:::emc_kernel_stats(FALSE)
  off <- do.call(EMC2:::dpcounter, args)
  EMC2:::emc_kernel_stats(TRUE)
  on <- do.call(EMC2:::dpcounter, args)
  expect_identical(on, off)
})
