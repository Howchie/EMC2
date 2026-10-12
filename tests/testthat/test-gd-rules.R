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

test_that("sampled-scale R-hats are not range-transformed", {
  emc <- samples_LNR
  idx <- which(emc[[1]]$samples$stage == "sample")
  # Every draw positive, one chain wandering near zero: coda would log it.
  for (j in seq_along(emc)) {
    x <- emc[[j]]$samples$alpha[1, 1, idx]
    emc[[j]]$samples$alpha[1, 1, idx] <- x - min(x) + if (j == 1) 1e-3 else 0.5
  }
  mcl <- get_pars(emc, selection = "alpha", stage = "sample", by_subject = TRUE)[[1]]
  raw <- EMC2:::gelman_diag_robust(mcl, transform = FALSE)
  logged <- EMC2:::gelman_diag_robust(mcl, transform = TRUE)
  expect_gt(abs(logged[1] - raw[1]), 1e-3)
  got <- gd_summary.emc(emc, selection = "alpha", stage = "sample", stat = NULL,
                        version = "old", digits = 10)
  expect_equal(unname(got[[1]][1]), unname(raw[1]), tolerance = 1e-8)
  # Positive-support selections keep the transform.
  vs <- get_pars(emc, selection = "sigma2", stage = "sample")[[1]]
  expect_equal(unlist(gd_summary.emc(emc, selection = "sigma2", stage = "sample",
                                     stat = NULL, digits = 10)),
               EMC2:::gelman_diag_robust(vs, transform = TRUE), tolerance = 1e-8,
               ignore_attr = TRUE)
})
