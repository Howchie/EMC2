# C12: one core means one call.
#
# `nrow(proposals) <= r_cores` sent every multi-proposal call into the
# split/subset/lapply branch even at r_cores = 1, where the split is a no-op:
# `.split_work_indices(n, 1)` is all ones, subsetting by it copies the whole
# matrix, and `auto_mclapply` over a single index builds a list to unlist it.

test_that("the split at one core was a no-op that still copied", {
  # The premise, stated so a change to .split_work_indices cannot quietly
  # invalidate the reasoning above.
  expect_identical(unique(EMC2:::.split_work_indices(100L, 1L)), 1L)
  expect_identical(length(EMC2:::.split_work_indices(100L, 1L)), 100L)
})

test_that("serial and split routes give identical likelihoods", {
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 300L, n_particles = 12L, cells = 8L)
    serial <- EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = 1)
    # More cores than proposals also takes the direct call; between them lies
    # the genuinely parallel route.
    direct <- EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = 20)
    split <- EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = 3)
    expect_bit_identical(direct, serial, paste(model, "- direct vs serial"))
    expect_bit_identical(split, serial, paste(model, "- split vs serial"))
    expect_true(all(is.finite(serial)), info = model)
    expect_identical(length(serial), nrow(fx$prop))
  }
})

test_that("the serial branch holds for one proposal and for many", {
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 1L, cells = 4L)
  expect_identical(length(EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model,
                                                 r_cores = 1)), 1L)
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 50L, cells = 4L)
  one <- EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = 1)
  expect_identical(length(one), 50L)
  # Order is preserved: the split branch reassembles by index, and the direct
  # branch must not differ from it anywhere.
  expect_bit_identical(one, EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model,
                                                   r_cores = 4),
                       "50 proposals, 1 core vs 4")
})
