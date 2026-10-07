gd_of <- function(emc, select) {
  g <- unlist(gd_summary.emc(emc, selection = select, stage = "sample",
                             omit_mpsrf = TRUE, stat = NULL, version = "old"))
  g[is.na(g)] <- Inf
  g
}

run_check <- function(mean_gd = NULL, max_gd = NULL, selection = c("alpha", "mu")) {
  EMC2:::check_gd(samples_LNR, stage = "sample", max_gd = max_gd, mean_gd = mean_gd,
                  omit_mpsrf = TRUE, trys = 0, verbose = FALSE, selection = selection,
                  iter = 50, truncate = FALSE, rhat_version = "old")
}

test_that("a scalar R-hat rule is applied to each selection separately", {
  a <- gd_of(samples_LNR, "alpha"); m <- gd_of(samples_LNR, "mu")
  out <- run_check(mean_gd = 10)
  expect_equal(out$gd_stats$selection, c("alpha", "mu"))
  expect_equal(out$gd_stats$statistic, c("mean", "mean"))
  expect_equal(out$gd_stats$value, c(mean(a), mean(m)))
  expect_true(out$gd_done)
  # Passing on one selection does not carry the other.
  thr <- c(alpha = mean(a) + 1e-8, mu = mean(m) - 1e-8)
  expect_false(run_check(mean_gd = thr)$gd_done)
  expect_true(run_check(mean_gd = thr + c(0, 2e-8))$gd_done)
})

test_that("named thresholds mix statistics and skip unnamed selections", {
  a <- gd_of(samples_LNR, "alpha"); m <- gd_of(samples_LNR, "mu")
  out <- run_check(mean_gd = c(alpha = 10), max_gd = c(mu = 10))
  expect_equal(out$gd_stats$statistic, c("mean", "max"))
  expect_equal(out$gd_stats$value, c(mean(a), max(m)))
  # Only mu is ruled, so alpha's R-hats are not computed.
  out <- run_check(max_gd = c(mu = max(m) + 1e-8))
  expect_equal(out$gd_stats$selection, "mu")
  expect_equal(out$gd, unname(m))
  expect_true(out$gd_done)
  expect_false(run_check(max_gd = c(mu = max(m) - 1e-8))$gd_done)
})

test_that("R-hat thresholds are validated against the selection", {
  gsc <- EMC2:::get_stop_criteria
  expect_error(gsc("burn", list(mean_gd = c(Sigma = 1.1), selection = c("alpha", "mu")), "standard"),
               "not in selection")
  expect_error(gsc("burn", list(mean_gd = c(1.1, 1.2)), "standard"), "single number")
  expect_silent(gsc("burn", list(mean_gd = c(alpha = 1.1), max_gd = c(mu = 1.2)), "standard"))
})
