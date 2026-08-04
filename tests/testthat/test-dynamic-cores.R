# Dynamic reallocation of cores from finished chains to still-running ones.

test_that("the arena is only built when it can help", {
  withr::local_options(list(emc2.dynamic_cores = TRUE))
  expect_null(EMC2:::.emc_core_ctl(1, 8))                     # nothing to donate
  withr::local_options(list(emc2.dynamic_cores = FALSE))
  expect_null(EMC2:::.emc_core_ctl(3, 8))                     # off by default
})

test_that("a chain claims the cores released by finished siblings", {
  skip_on_os("windows")
  withr::local_options(list(emc2.dynamic_cores = TRUE))
  ctl <- EMC2:::.emc_core_ctl(3, 8)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  expect_equal(ctl$total, 24)

  # Nobody has finished: stay on the static share.
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 8)

  file.create(file.path(ctl$dir, "done_1"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 12)   # 24 / 2 remaining

  file.create(file.path(ctl$dir, "done_2"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 24)   # last chain takes the lot

  # Never drops below the chain's own share, whatever the bookkeeping says.
  file.create(file.path(ctl$dir, "done_3"))
  expect_gte(EMC2:::.emc_cores_now(ctl, 8), 8)
})

test_that("a chain releases its cores even when it fails", {
  skip_on_os("windows")
  withr::local_options(list(emc2.dynamic_cores = TRUE))
  ctl <- EMC2:::.emc_core_ctl(3, 4)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  boom <- function() {
    on.exit(EMC2:::.emc_core_release(ctl), add = TRUE)
    stop("chain failed")
  }
  expect_error(boom(), "chain failed")
  expect_length(list.files(ctl$dir), 1)
})

test_that("disabled control leaves the static budget untouched", {
  expect_equal(EMC2:::.emc_cores_now(NULL, 8), 8)
  expect_silent(EMC2:::.emc_core_release(NULL))
})
