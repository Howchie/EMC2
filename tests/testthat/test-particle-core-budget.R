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

test_that("nested likelihood workers fit inside the chain budget", {
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 4L,
                                 total_cores = 4L),
    list(subject = 1L, likelihood = 4L)
  )
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 8L, r_cores = 4L,
                                 total_cores = 16L),
    list(subject = 4L, likelihood = 4L)
  )
  # An inner request larger than the chain allocation is capped rather than
  # multiplying past it.
  capped <- EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 8L,
                                         total_cores = 4L)
  expect_equal(capped$subject * capped$likelihood, 4L)
})

test_that("proposal likelihood splits never create empty workers", {
  idx <- EMC2:::.split_work_indices(5L, 4L)
  expect_equal(idx, c(1L, 1L, 2L, 3L, 4L))
  expect_true(all(tabulate(idx, nbins = 4L) > 0L))
})

test_that("log_likelihood_joint forwards r_cores to calc_ll_manager", {
  formals_joint <- formals(EMC2:::log_likelihood_joint)
  expect_true("r_cores" %in% names(formals_joint))
  expect_identical(formals_joint$r_cores, 1)
})

test_that("the inner fan-out absorbs cores the outer split would truncate away", {
  pb <- EMC2:::.particle_core_budget
  # 8 cores at r_cores = 3 used to give 2 subject workers x 3 = 6 of 8.
  b <- pb(112, n_cores = 8, r_cores = 3, total_cores = 8)
  expect_equal(b$subject * b$likelihood, 8)
  expect_gte(b$likelihood, 3)   # r_cores stays a lower bound

  # The r_cores = 1 default must not gain inner workers it never asked for.
  b1 <- pb(3, n_cores = 8, r_cores = 1, total_cores = 8)
  expect_equal(b1$likelihood, 1)
  expect_equal(b1$subject, 3)

  # Widening never breaks the budget, at any combination.
  for (ns in c(1, 3, 5, 112)) for (nc in c(1, 4, 8, 16)) for (rc in c(1, 2, 3, 4, 32)) {
    b <- pb(ns, n_cores = nc, r_cores = rc, total_cores = nc)
    expect_lte(b$subject * b$likelihood, nc)
    expect_gte(b$subject, 1L)
    expect_gte(b$likelihood, 1L)
  }
})
