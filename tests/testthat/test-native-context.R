# C11: give the native per-call state an explicit owner and lifetime.
#
# The mixed (censored/truncated) path holds its trial partitions as const views
# into state owned by the call, avoiding copies without extending their lifetime.

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
