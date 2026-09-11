# C10: an explicit update mask for blocked proposals.
#
# `make_pt_mapper` decided a design was invariant from column NAMES -- none of
# its coefficients appear in the proposal matrix. Under a blocked proposal most
# of those named columns hold the current value repeated, so every design that
# read only repeated coordinates was re-mapped per particle for nothing.
#
# The two things that must hold: the mask changes no number, and a wrong mask is
# refused rather than acted on.

mask_args <- function(fx, varying = NULL) {
  c(audit_ll_args(fx), list(type = fx$model()$c_name, min_ll = log(1e-10),
                            varying = varying))
}

test_that("a mask changes no likelihood", {
  for (model in c("RDM", "LBA", "DDM")) {
    for (cells in c(4L, 64L)) {
      fx <- audit_fixture(model, n_trials = 400L, n_particles = 6L, cells = cells)
      # Freeze every coordinate but the first: the rest of the proposal matrix
      # holds one repeated value, which is what a blocked proposal looks like.
      prop <- fx$prop
      keep <- rep(FALSE, ncol(prop))
      keep[1L] <- TRUE
      prop[, !keep] <- rep(prop[1L, !keep], each = nrow(prop))
      fx$prop <- prop
      lbl <- paste(model, "cells =", cells)
      unmasked <- do.call(EMC2:::calc_ll_oo, mask_args(fx))
      masked <- do.call(EMC2:::calc_ll_oo, mask_args(fx, keep))
      expect_bit_identical(masked, unmasked, lbl)
      expect_true(all(is.finite(masked)), info = lbl)
    }
  }
})

test_that("an all-TRUE mask is the behaviour that existed before", {
  fx <- audit_fixture("RDM", n_trials = 400L, n_particles = 5L, cells = 8L)
  none <- do.call(EMC2:::calc_ll_oo, mask_args(fx))
  all_true <- do.call(EMC2:::calc_ll_oo,
                      mask_args(fx, rep(TRUE, ncol(fx$prop))))
  expect_bit_identical(all_true, none, "all-TRUE mask vs no mask")
})

test_that("a mask that lies is refused, not acted on", {
  # The danger of an explicit mask is that a wrong one silently freezes a
  # coordinate the sampler meant to move. The check costs O(particles) against
  # mapping's O(particles x trials), so it is always on.
  fx <- audit_fixture("RDM", n_trials = 200L, n_particles = 4L, cells = 4L)
  lying <- rep(FALSE, ncol(fx$prop))   # nothing varies, but everything does
  expect_error(do.call(EMC2:::calc_ll_oo, mask_args(fx, lying)),
               "marked unchanged")
  # And the message names the column and the particle, so it can be chased.
  err <- tryCatch(do.call(EMC2:::calc_ll_oo, mask_args(fx, lying)),
                  error = conditionMessage)
  expect_match(err, colnames(fx$prop)[1L], fixed = TRUE)
})

test_that("a mask of the wrong length is refused", {
  fx <- audit_fixture("RDM", n_trials = 200L, n_particles = 4L, cells = 4L)
  expect_error(do.call(EMC2:::calc_ll_oo, mask_args(fx, c(TRUE, FALSE))),
               "varying mask has")
})

test_that("NA in the mask is read as varying", {
  # An unknown coordinate must be mapped, not frozen: the safe reading of "I do
  # not know" is "assume it moves".
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 5L, cells = 8L)
  m <- rep(NA, ncol(fx$prop))
  expect_bit_identical(do.call(EMC2:::calc_ll_oo, mask_args(fx, m)),
                       do.call(EMC2:::calc_ll_oo, mask_args(fx)),
                       "all-NA mask")
})

test_that("a one-particle call is unaffected by the mask", {
  # The invariant lane only engages for more than one particle; a single
  # particle has nothing to reuse across.
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 1L, cells = 8L)
  keep <- rep(FALSE, ncol(fx$prop))
  keep[1L] <- TRUE
  expect_bit_identical(do.call(EMC2:::calc_ll_oo, mask_args(fx, keep)),
                       do.call(EMC2:::calc_ll_oo, mask_args(fx)),
                       "single particle")
})

test_that("the mask reaches the mapper through calc_ll_manager", {
  fx <- audit_fixture("RDM", n_trials = 400L, n_particles = 6L, cells = 8L)
  prop <- fx$prop
  keep <- rep(FALSE, ncol(prop))
  keep[1L] <- TRUE
  prop[, !keep] <- rep(prop[1L, !keep], each = nrow(prop))
  a <- EMC2:::calc_ll_manager(prop, fx$dadm, fx$model, r_cores = 1)
  b <- EMC2:::calc_ll_manager(prop, fx$dadm, fx$model, r_cores = 1,
                              varying = keep)
  expect_bit_identical(b, a, "via calc_ll_manager")
  # And a lying mask still fails at that level rather than being dropped.
  expect_error(EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = 1,
                                      varying = rep(FALSE, ncol(fx$prop))),
               "marked unchanged")
})
