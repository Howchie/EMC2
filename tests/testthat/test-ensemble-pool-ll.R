# Ensemble pool likelihoods: failures are zero-density candidates counted under
# the failure policy, and computing them uses no random numbers.

local_rng_guard()

test_that("checked evaluation turns likelihood failures into recorded rejections", {
  EMC2:::.emc_reject_reset("ensemble")
  assign("announced", character(), envir = EMC2:::.emc_profile_state)
  props <- matrix(0, 1, 2)
  local_mocked_bindings(calc_ll_manager = function(...) -Inf, .package = "EMC2")
  expect_identical(EMC2:::.emc_ensemble_ll_checked(props, NULL, NULL), -Inf)
  local_mocked_bindings(calc_ll_manager = function(...) stop("non-finite likelihood"),
                        .package = "EMC2")
  expect_silent(expect_identical(EMC2:::.emc_ensemble_ll_candidate(props, NULL, NULL), -Inf))
  expect_equal(EMC2:::.emc_reject_counts("ensemble")[["numerical"]], 1L)
  local_mocked_bindings(calc_ll_manager = function(...) stop("broken evaluator"),
                        .package = "EMC2")
  expect_warning(expect_identical(EMC2:::.emc_ensemble_ll_candidate(props, NULL, NULL), -Inf),
                 "ensemble update failed")
  expect_equal(EMC2:::.emc_reject_counts("ensemble")[["unknown"]], 1L)
  local_mocked_bindings(calc_ll_manager = function(...) NaN, .package = "EMC2")
  expect_silent(expect_identical(EMC2:::.emc_ensemble_ll_candidate(props, NULL, NULL), -Inf))
})

test_that("pool likelihoods fill any subject the subject step missed, without RNG", {
  ll_mock <- function(proposals, dadm, model, ...) -0.5 * rowSums(sweep(proposals, 2, dadm$y)^2)
  local_mocked_bindings(calc_ll_manager = ll_mock, .package = "EMC2")
  s <- list(n_subjects = 2L, data = list(list(y = c(0, 0)), list(y = c(1, 1))),
            model = NULL, marginalise = NULL)
  rows <- list(matrix(c(1, 2, 3, 4), 2), matrix(c(0, 1, 1, 0), 2))
  pools <- list(rows = rows, ll = list(c(-9, -9), NULL))
  set.seed(24); before <- .Random.seed
  ll <- EMC2:::.emc_ensemble_pool_ll(s, pools)
  expect_identical(.Random.seed, before)
  expect_equal(ll[, 1], c(-9, -9))                     # returned by the subject step
  expect_equal(ll[, 2], ll_mock(rows[[2]], s$data[[2]]))  # filled in here
})
