test_that("LBA log outputs remain finite past natural-scale underflow", {
  t <- c(0.02, 0.01)
  A <- c(0.3, 0.3)
  b <- c(1.5, 1.5)
  v <- c(1, 1)
  sv <- c(1, 1)

  density <- EMC2:::dlba(t, A, b, v, sv, posdrift = TRUE, log_out = FALSE)
  log_density <- EMC2:::dlba(t, A, b, v, sv, posdrift = TRUE, log_out = TRUE)
  cdf <- EMC2:::plba(t, A, b, v, sv, posdrift = TRUE, log_out = FALSE)
  log_cdf <- EMC2:::plba(t, A, b, v, sv, posdrift = TRUE, log_out = TRUE)

  expect_true(all(density == 0))
  expect_true(all(cdf == 0))
  expect_true(all(is.finite(log_density)))
  expect_true(all(is.finite(log_cdf)))
  expect_true(all(log_density < log(.Machine$double.xmin)))
  expect_true(all(log_cdf < log(.Machine$double.xmin)))
})

test_that("LBA natural and log outputs agree before underflow", {
  t <- c(0.2, 0.5, 1, 2)
  density <- EMC2:::dlba(t, rep(0.3, 4), rep(1.5, 4), rep(0.3, 4),
                         rep(1, 4), TRUE, FALSE)
  log_density <- EMC2:::dlba(t, rep(0.3, 4), rep(1.5, 4), rep(0.3, 4),
                             rep(1, 4), TRUE, TRUE)
  cdf <- EMC2:::plba(t, rep(0.3, 4), rep(1.5, 4), rep(0.3, 4),
                     rep(1, 4), TRUE, FALSE)
  log_cdf <- EMC2:::plba(t, rep(0.3, 4), rep(1.5, 4), rep(0.3, 4),
                         rep(1, 4), TRUE, TRUE)

  expect_equal(log_density, log(density), tolerance = 1e-12)
  expect_equal(log_cdf, log(cdf), tolerance = 1e-12)
})
