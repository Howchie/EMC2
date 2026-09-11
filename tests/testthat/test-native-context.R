# C11: give the native per-call state an explicit owner and lifetime.
#
# Two halves. The mixed (censored/truncated) path held its trial partitions in
# std::vectors that were copy-assigned from state the call already owned, under
# a comment claiming "no copy"; they are const views now. And the compiled plan
# is separated from the per-call state in the type, with an explicit key for
# what would have to invalidate a cache of it.

key_args <- function(fx) {
  m <- fx$model()
  list(data = fx$dadm, constants = if (identical(fx$constants, NA)) numeric(0)
                                   else fx$constants,
       designs = audit_designs(fx), bounds = m$bound, transforms = m$transform,
       pretransforms = m$pre_transform, trend = m$trend, p_types = fx$p_types)
}

test_that("the key is equal for identical inputs", {
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 3L, cells = 6L)
  a <- do.call(EMC2:::emc_pt_invalidation_key, key_args(fx))
  b <- do.call(EMC2:::emc_pt_invalidation_key, key_args(fx))
  expect_identical(a, b)
  # It reports values, not a hash: a cache must compare what actually differs
  # rather than trust that two plans could not collide.
  expect_true(all(c("n_trials", "columns", "levels", "designs", "constants",
                    "transforms", "pretransforms", "bounds", "trend", "p_types",
                    "guess_window") %in% names(a)))
  expect_equal(a$n_trials, nrow(fx$dadm))
})

test_that("the key changes for everything the audit says invalidates a plan", {
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 3L, cells = 6L)
  base <- do.call(EMC2:::emc_pt_invalidation_key, key_args(fx))

  # Data shape.
  smaller <- key_args(fx); smaller$data <- fx$dadm[seq_len(nrow(fx$dadm) - 2L), ]
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, smaller), base))

  # Levels behind the designs -- the same column, a different factor.
  lv <- key_args(fx)
  d2 <- fx$dadm
  d2$lR <- factor(as.character(d2$lR), levels = c(levels(d2$lR), "extra"))
  lv$data <- d2
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, lv), base))

  # Designs.
  des <- key_args(fx)
  des$designs$v <- des$designs$v[, 1L, drop = FALSE]
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, des), base))

  # Constants, transforms, bounds, p_types.
  cst <- key_args(fx); cst$constants <- c(cst$constants, extra = 1)
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, cst), base))
  tf <- key_args(fx); tf$transforms <- list(func = c(v = "exp"))
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, tf), base))
  bd <- key_args(fx); bd$bounds$minmax[1L, 1L] <- bd$bounds$minmax[1L, 1L] - 1
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, bd), base))
  pt <- key_args(fx); pt$p_types <- rev(pt$p_types)
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, pt), base))

  # A trend where there was none.
  tr <- key_args(fx); tr$trend <- list(kind = "placeholder")
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, tr), base))

  # And the guess window, which is data-attached rather than an argument.
  gw <- key_args(fx)
  d3 <- fx$dadm
  attr(d3, "guess_window") <- c(0.1, 0.9)
  gw$data <- d3
  expect_false(identical(do.call(EMC2:::emc_pt_invalidation_key, gw), base))
})

test_that("the mixed path still gives the same likelihoods", {
  # The partitions that became const views are the censored/truncated path's
  # finite/other/nogo trial splits, read by every particle. Exercise a fixture
  # that actually takes that path.
  skip_if_not(exists("audit_fixture"))
  set.seed(4)
  fx <- audit_fixture("RDM", n_trials = 400L, n_particles = 6L, cells = 6L)
  dat <- fx$data
  cens <- sample(nrow(dat), floor(0.2 * nrow(dat)))
  dat$rt[cens] <- Inf
  dat$R[cens] <- NA
  mixed <- audit_fixture("RDM", n_trials = 400L, n_particles = 6L, cells = 6L,
                         data = dat)
  ll <- audit_ll_direct(mixed)
  expect_identical(length(as.numeric(ll)), nrow(mixed$prop))
  expect_true(all(is.finite(as.numeric(ll))))
  # Repeated calls agree exactly: a view whose owner had gone would not.
  expect_bit_identical(audit_ll_direct(mixed), ll, "mixed path, repeated")
  # And the two mapping routes agree on it, which is the C2 gate applied to the
  # path this commit touched.
  expect_mapper_agrees(mixed, label = "mixed path")
})

test_that("a censored fixture exercises the partitioned branch", {
  # Guard against the test above silently drifting onto the all-finite fast
  # path, where none of this code runs.
  set.seed(5)
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 2L, cells = 4L)
  dat <- fx$data
  dat$rt[seq(1L, nrow(dat), by = 5L)] <- Inf
  dat$R[seq(1L, nrow(dat), by = 5L)] <- NA
  mixed <- audit_fixture("RDM", n_trials = 300L, n_particles = 2L, cells = 4L,
                         data = dat)
  expect_false(isTRUE(attr(mixed$dadm, "emc2_all_finite_trials")))
})
