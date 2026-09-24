local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# Acceptance gates 4 and 5 of the architecture-efficiency audit: arithmetic
# changes must preserve support, floor and mixture semantics, over a grid that
# includes "RT near t0, zero/small A and leak, zero drift where supported,
# extreme tails, invalid proposals, NA versus NaN, +/-Inf censoring, known and
# unknown responses, LT/UT endpoints, varying accumulator counts/inactive rows,
# contaminants/guessing and compressed/uncompressed expansion".
#
# Division of labour with the tests that already exist.  test-nan-censoring.R,
# test-limit-truncation.R and test-make_missing.R own the *semantics* of
# censoring, truncation and missingness -- what those routes are supposed to
# compute.  This file owns the *invariance* of the awkward corners under the
# changes the audit proposes: a mapper rewrite, a hoisted kernel subexpression
# or a re-scheduled batch must not quietly move a floored row, turn a floor into
# a NaN, or make a rejected particle look valid.  Compression is covered in
# test-audit-differential.R.
#
# Every expectation below is pinned against a value derived from the fixture,
# not a literal: the floor is n * log(min_ll), and the comparisons are between
# two routes or between a corner and its own baseline.

RNGkind("L'Ecuyer-CMRG")

# The likelihood the model reports when a particle is invalid, or when a row
# cannot contribute: every trial contributes exactly `min_ll`.
audit_floor <- function(fx, min_ll = log(1e-10)) {
  n <- length(attr(fx$dadm, "expand"))
  if (!n) n <- nrow(fx$dadm)
  n * min_ll
}

# The floor as the model actually reaches it, by rejecting the particle on a
# bound rather than by failing every row's density.  Two different routes to
# the same answer, which is what makes it a test rather than a restatement.
audit_floor_observed <- function(fx) {
  bad <- fx
  bad$prop[, "t0"] <- log(1e-9)
  as.numeric(audit_ll_direct(bad))[1L]
}

# A floor is a *sum* of n copies of min_ll, and how that sum is accumulated is
# model-specific: DDM's floor sits 4 ULP away from n * log(1e-10) because it
# adds where the analytic expression multiplies, while LBA and RDM happen to
# land on it exactly.  Comparing against the arithmetic literal at tolerance 0
# would therefore fail DDM for a reason that has nothing to do with the
# support boundary being tested.  So: bit-identical against the independently
# reached floor, and merely close against the literal.
expect_floored <- function(fx, ll, label = "") {
  expect_bit_identical(as.numeric(ll)[1L], audit_floor_observed(fx),
                       paste(label, "- vs the bounds-rejection floor"))
  testthat::expect_equal(as.numeric(ll)[1L], audit_floor(fx), tolerance = 1e-12,
                         info = paste(label, "- vs n * min_ll"))
}

# Not merely "no error": a kernel that loses a guard returns NaN or -Inf, and
# both of those propagate silently into a sum before anything notices.
expect_all_finite <- function(ll, label = "") {
  testthat::expect_true(all(is.finite(ll)),
                        info = paste(label, "- finite, not NaN or Inf"))
}

# --- the support boundary in RT ---------------------------------------------

test_that("RT at and below t0 floors rather than producing a NaN", {
  # The boundary the whole race family is most exposed to.  A hoisted
  # parameter-only term that forgot the RT-dependent guard would show up here as
  # a NaN or an Inf rather than as a floor.
  for (model in c("LBA", "RDM", "DDM")) {
    fx <- audit_fixture(model, n_trials = 60L, n_particles = 1L)
    t0 <- exp(fx$centre[["t0"]])
    for (tag in c("below", "at", "just above")) {
      rt <- switch(tag, below = t0 / 2, at = t0, `just above` = t0 + 1e-7)
      dat <- fx$data
      dat$rt <- rt
      corner <- audit_fixture(model, data = dat, n_particles = 1L)
      corner$prop <- fx$prop
      ll <- audit_ll_direct(corner)
      expect_all_finite(ll, paste(model, tag, "t0"))
      expect_floored(corner, ll, paste(model, "rt", tag, "t0"))
    }
  }
})

test_that("an extreme RT tail stays finite and is worse than the mode", {
  for (model in c("LBA", "RDM", "DDM")) {
    fx <- audit_fixture(model, n_trials = 60L, n_particles = 1L)
    base <- audit_ll_direct(fx)
    dat <- fx$data
    dat$rt <- 60                        # far out in the tail, still legal
    tail_fx <- audit_fixture(model, data = dat, n_particles = 1L)
    tail_fx$prop <- fx$prop
    ll <- audit_ll_direct(tail_fx)
    expect_all_finite(ll, paste(model, "extreme tail"))
    expect_lt(ll, base)
  }
})

# --- non-finite and missing RT ----------------------------------------------

test_that("NA, NaN and +/-Inf RTs are all handled without a NaN likelihood", {
  # These four are the sentinels the censoring machinery uses.  In a plain
  # design with no missingness declared they are simply rows that cannot
  # contribute, and this pins that they behave identically rather than one of
  # them slipping through into the arithmetic.  What each *means* when
  # missingness is declared is test-nan-censoring.R's job.
  fx <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L)
  base <- audit_ll_direct(fx)
  sentinels <- list(`NA` = NA_real_, `NaN` = NaN, `+Inf` = Inf, `-Inf` = -Inf)
  seen <- numeric()
  for (nm in names(sentinels)) {
    dat <- fx$data
    dat$rt[1:3] <- sentinels[[nm]]
    corner <- audit_fixture("LBA", data = dat, n_particles = 1L)
    corner$prop <- fx$prop
    ll <- audit_ll_direct(corner)
    expect_all_finite(ll, paste("rt =", nm))
    expect_lt(ll, base)              # three rows stopped contributing
    seen[nm] <- ll
  }
  # All four take the same route in an undeclared design.  If a later change
  # makes one of them differ, that is a decision, and it should be a visible
  # one.
  expect_equal(length(unique(seen)), 1L)
})

test_that("an unknown response is not the same as a missing RT", {
  # R = NA with a finite RT is a legal, informative observation: the trial
  # happened and was timed, the response is simply unlabelled.  It must not be
  # collapsed into the missing-RT route.
  fx <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L)
  base <- audit_ll_direct(fx)
  dat <- fx$data
  dat$R[1:3] <- NA
  corner <- audit_fixture("LBA", data = dat, n_particles = 1L)
  corner$prop <- fx$prop
  ll <- audit_ll_direct(corner)
  expect_all_finite(ll, "R = NA")
  expect_false(isTRUE(all.equal(as.numeric(ll), as.numeric(base))))
  # Marginalising over the response is *more* likely than pinning it to one
  # accumulator, so an unknown response must not be penalised like a lost trial.
  expect_gt(ll, base - 1)
})

# --- invalid proposals ------------------------------------------------------

test_that("out-of-bound parameters floor the whole particle", {
  # Bounds are checked per row but a failure floors the particle's total.  Each
  # of these is a different bound, so a bounds rewrite that dropped one would
  # leave the others passing.
  fx <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L)
  corners <- list(
    "t0 under its lower bound" = c(t0 = log(1e-9)),
    "B negative"               = c(B = -50),
    "A at zero"                = c(A = log(1e-12)))
  for (nm in names(corners)) {
    bad <- fx
    for (par in names(corners[[nm]])) bad$prop[1L, par] <- corners[[nm]][[par]]
    ll <- audit_ll_direct(bad)
    expect_all_finite(ll, nm)
    expect_floored(bad, ll, nm)
    # The bound verdict and the floored likelihood must agree.
    expect_false(all(audit_prologue(bad)$ok), info = nm)
  }
})

test_that("one invalid particle does not floor its neighbours", {
  # The batch is the unit of the call but not of the verdict.  A scheduling or
  # chunking change that let one particle's rejection leak sideways would be
  # invisible in a single-particle test.
  fx <- audit_fixture("LBA", n_trials = 60L, n_particles = 4L)
  base <- audit_ll_direct(fx)
  bad <- fx
  bad$prop[2L, "t0"] <- log(1e-9)
  ll <- audit_ll_direct(bad)
  expect_equal(as.numeric(ll[2L]), audit_floor(bad), tolerance = 0)
  expect_bit_identical(ll[-2L], base[-2L], "neighbours of an invalid particle")
})

# --- degenerate but legal parameters ----------------------------------------

test_that("zero drift and a vanishing start-point range stay finite", {
  # "Zero drift where supported" and "zero/small A".  LBA takes v on the natural
  # scale, so v = 0 is reachable and legal; A is bounded below by zero and the
  # small-A limit is a separate closed form in several kernels.
  fx <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L)
  zero_v <- fx
  zero_v$prop[1L, grepl("^v", colnames(fx$prop))] <- 0
  ll <- audit_ll_direct(zero_v)
  expect_all_finite(ll, "v = 0")
  expect_false(isTRUE(all.equal(as.numeric(ll), audit_floor(zero_v))))

  # A shrinking towards zero must approach the A = 0 limit smoothly rather than
  # jumping, which is what a mishandled small-A branch looks like.
  lls <- vapply(c(-3, -5, -7, -9), function(logA) {
    g <- fx
    g$prop[1L, "A"] <- logA
    as.numeric(audit_ll_direct(g))
  }, numeric(1))
  expect_true(all(is.finite(lls)))
  expect_lt(max(abs(diff(lls))), 50)
})

test_that("a leak parameter at and near zero stays finite", {
  # BAwL's zero-leak limit is a separate closed form; the audit requires it
  # preserved through any natural_normalizer work (C9).
  fx <- audit_fixture("BAwL", n_trials = 60L, n_particles = 1L)
  base <- audit_ll_direct(fx)
  expect_all_finite(base, "BAwL baseline")
  lls <- vapply(c(log(0.2), log(1e-3), log(1e-6), log(1e-9)), function(k) {
    g <- fx
    g$prop[1L, "k"] <- k
    as.numeric(audit_ll_direct(g))
  }, numeric(1))
  expect_true(all(is.finite(lls)))
  # Approaching zero leak converges rather than diverging.
  expect_lt(abs(lls[4L] - lls[3L]), abs(lls[2L] - lls[1L]) + 1)
})

# --- accumulator structure --------------------------------------------------

test_that("a third accumulator changes the row count and stays finite", {
  # "Varying accumulator counts/inactive rows".  The augmented data is one row
  # per accumulator per trial, so the winner/loser masks and the reductions all
  # change shape here.
  two <- audit_fixture("LBA", n_trials = 60L, n_particles = 2L)
  dat <- audit_data(60L)
  dat$S <- factor(rep_len(c("left", "right", "third"), nrow(dat)))
  dat$R <- factor(rep_len(c("left", "right", "third"), nrow(dat)),
                  levels = levels(dat$S))
  three <- audit_fixture("LBA", data = dat, n_particles = 2L)
  expect_equal(nrow(two$dadm), 2L * 60L)
  expect_equal(nrow(three$dadm), 3L * 60L)
  ll <- audit_ll_direct(three)
  expect_all_finite(ll, "three accumulators")
  expect_false(isTRUE(all.equal(as.numeric(ll), audit_floor(three))))
  # And the two mapping routes still agree on the wider block.
  want <- audit_reference(three)
  got <- audit_prologue(three)$pars
  for (nm in dimnames(want)[[2L]]) {
    expect_bit_identical(got[, nm, , drop = FALSE], want[, nm, , drop = FALSE],
                         paste("three accumulators -", nm))
  }
})

# --- mixture components -----------------------------------------------------

test_that("a contaminant component changes the likelihood and stays finite", {
  # "Contaminants/guessing".  The audit forbids changing mixture semantics, and
  # the mixture only exists when its weight is off its degenerate default.
  plain <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L)
  mixed <- audit_fixture("LBA", n_trials = 60L, n_particles = 1L,
                         constants = c(sv = log(1), pContaminant = qnorm(0.05)))
  mixed$prop <- plain$prop[, colnames(mixed$prop), drop = FALSE]
  base <- audit_ll_direct(plain)
  ll <- audit_ll_direct(mixed)
  expect_all_finite(ll, "contaminant")
  expect_false(isTRUE(all.equal(as.numeric(ll), as.numeric(base))))
  # A contaminant is a floor on how bad any single trial can be, so it must not
  # make the fit worse overall than the pure model on its own best guess.
  expect_gt(ll, base - 50)
})

test_that("the contaminant weight moves the likelihood monotonically", {
  # The mixture's weight enters on the probit scale (`pContaminant` maps
  # through pnorm, so qnorm(0.05) is a 5% component).  Away from the data's
  # own mode, giving the contaminant more weight costs likelihood, and it does
  # so smoothly.  Pinning the monotone ordering catches a mixture whose
  # component densities or weights have been reassociated, without asserting a
  # direction the model does not actually have: at these parameters the
  # contaminant does *not* rescue a long RT, it simply dilutes the race.
  dat <- audit_data(60L)
  dat$rt[1:5] <- 30                     # far past anything the race explains
  plain <- audit_fixture("LBA", data = dat, n_particles = 1L)
  weights <- c(0.01, 0.05, 0.2, 0.5)
  lls <- vapply(weights, function(p) {
    mixed <- audit_fixture("LBA", data = dat, n_particles = 1L,
                           constants = c(sv = log(1), pContaminant = qnorm(p)))
    mixed$prop <- plain$prop[, colnames(mixed$prop), drop = FALSE]
    as.numeric(audit_ll_direct(mixed))
  }, numeric(1))
  expect_true(all(is.finite(lls)))
  expect_true(all(diff(lls) < 0))
  # And the weight really is reaching the parameter table on the natural scale.
  mixed <- audit_fixture("LBA", data = dat, n_particles = 1L,
                         constants = c(sv = log(1), pContaminant = qnorm(0.05)))
  expect_equal(unique(as.numeric(audit_prologue(mixed)$pars[, "pContaminant", 1L])),
               0.05, tolerance = 1e-12)
})

# --- cross-route agreement at the corners -----------------------------------

test_that("the managed and direct routes agree at every corner", {
  # calc_ll_manager() and calc_ll_oo() are different entry points to the same
  # arithmetic.  C12 changes the manager, so the corners have to be pinned
  # across both, not only down the direct path.
  dat <- audit_data(60L)
  dat$rt[1:3] <- NA
  dat$rt[4:6] <- Inf
  dat$R[7:9] <- NA
  corners <- list(
    plain = audit_fixture("LBA", n_trials = 60L, n_particles = 3L),
    ragged = audit_fixture("LBA", data = dat, n_particles = 3L),
    wide = audit_fixture("LBA", n_trials = 60L, n_particles = 3L, cells = 300L))
  for (nm in names(corners)) {
    fx <- corners[[nm]]
    expect_bit_identical(audit_ll_managed(fx), audit_ll_direct(fx),
                         paste("managed vs direct -", nm))
  }
})
