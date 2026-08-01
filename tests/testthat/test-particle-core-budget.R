test_that("single-subject chains spend their core budget on proposals", {
  expect_identical(
    EMC2:::.particle_core_budget(1L, n_cores = 4L, r_cores = 1L),
    list(subject = 1L, likelihood = 4L)
  )
  expect_identical(
    EMC2:::.particle_core_budget(1L, n_cores = 2L, r_cores = 4L),
    list(subject = 1L, likelihood = 4L)
  )
})

test_that("multi-subject chains retain outer parallelism", {
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 1L),
    list(subject = 4L, likelihood = 1L)
  )
})

test_that("log_likelihood_joint forwards r_cores to calc_ll_manager", {
  formals_joint <- formals(EMC2:::log_likelihood_joint)
  expect_true("r_cores" %in% names(formals_joint))
  expect_identical(formals_joint$r_cores, 1)
})

