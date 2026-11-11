test_that("PSIS warming is once per process and preserves the RNG", {
  state <- EMC2:::.emc_psis_warm_state
  old <- state$pid
  on.exit(state$pid <- old, add = TRUE)
  state$pid <- NULL
  set.seed(812)
  before <- .Random.seed
  EMC2:::.emc_warm_psis()
  expect_identical(.Random.seed, before)
  expect_identical(state$pid, Sys.getpid())
  EMC2:::.emc_warm_psis()
  expect_identical(.Random.seed, before)
  expect_identical(state$pid, Sys.getpid())
})
