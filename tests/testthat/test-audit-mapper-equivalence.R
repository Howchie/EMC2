local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# Acceptance gate 2 of the architecture-efficiency audit: "Mapper-only changes
# must match mapped natural parameters and bounds against the existing
# independent get_pars_c_batch_wrapper_oo route."
#
# test-param-table-prologue.R already pins the ordinary design shapes.  What is
# added here is the matrix the audit's mapper commits (C7 dynamic cell scratch,
# C8 joint cell partition and scalar coefficients) will be judged against, and
# in particular the cell-count cutoff those commits exist to remove: 255, 256,
# 257 and much larger, plus the identical-input probe that crosses the cutoff
# without changing a single mapped value.
#
# These tests pass today.  That is the point -- they establish that the current
# routes agree, so that when C7 replaces the fixed 256-cell stack buffers the
# same assertions still have to hold.

RNGkind("L'Ecuyer-CMRG")

expect_mapper_agrees <- function(fx, designs = NULL, label = "") {
  got <- audit_prologue(fx, designs)$pars
  want <- audit_reference(fx, designs)
  # Compare column by column rather than as one flattened block, so a failure
  # names the parameter that moved instead of only a maximum difference.
  shared <- dimnames(want)[[2L]]
  testthat::expect_true(all(shared %in% dimnames(got)[[2L]]), info = label)
  testthat::expect_equal(dim(got)[c(1L, 3L)], dim(want)[c(1L, 3L)], info = label)
  for (nm in shared) {
    expect_bit_identical(got[, nm, , drop = FALSE], want[, nm, , drop = FALSE],
                         paste(label, "-", nm))
  }
}

# --- the cell-count cutoff --------------------------------------------------

test_that("mapping agrees across the old 256-cell cutoff", {
  # `EMC2_PT_MAX_CELLS = 256` used to be the admission test for the design-cell
  # path, and it was a cliff: above it, mapping fell back to
  # coefficient-by-trial loops and transforms and bounds fell back to row
  # resolution, for 113.50 ms against 54.83 ms on identical results.  C7
  # replaced it with a memory budget, so all three of these now take the same
  # route -- but they must still agree with the independent reference, and this
  # is where a regression that reintroduced a count would show up.
  for (cells in c(255L, 256L, 257L)) {
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 4L, cells = cells)
    expect_equal(nrow(audit_designs(fx)$v), cells)
    expect_mapper_agrees(fx, label = paste("cells =", cells))
  }
})
test_that("mapping agrees far above the cutoff", {
  # "Do not only raise the constant to 1,024 and move the problem": a fixture
  # well past any plausible replacement cap.
  fx <- audit_fixture("RDM", n_trials = 3000L, n_particles = 3L, cells = 1500L)
  expect_gt(nrow(audit_designs(fx)$v), 1024L)
  expect_mapper_agrees(fx, label = "cells = 1500")
})

test_that("an unreferenced design row changes no likelihood", {
  # The audit's stricter probe.  Append an *unused* row to a 256-row design:
  # `expand` never references it, so every trial, particle and mapped value is
  # unchanged and only the cell count moves.  Before C7 this crossed the
  # admission cliff and so compared the two routes; it now compares the cell
  # route with itself, which is still worth pinning -- a design row nobody reads
  # must not be able to move a number -- and the route comparison has moved to
  # the budget test below, where it can be made on demand.
  fx <- audit_fixture("RDM", n_trials = 800L, n_particles = 6L, cells = 256L)
  a <- audit_ll_args(fx)
  expect_equal(nrow(a$designs$v), 256L)
  b <- audit_pad_designs(a, "v", 1L)
  expect_equal(nrow(b$designs$v), 257L)

  la <- audit_call_ll(fx, a)
  lb <- audit_call_ll(fx, b)
  expect_informative_ll(fx, la, "cell route")
  expect_bit_identical(lb, la, "fallback route vs cell route")

  # And the mapped parameters themselves, not only the likelihood they reduce
  # to: a reduction can agree while the mapping underneath does not.
  pa <- do.call(EMC2:::pt_prologue_oo, c(a, list(return_pars = TRUE)))
  pb <- do.call(EMC2:::pt_prologue_oo, c(b, list(return_pars = TRUE)))
  expect_bit_identical(pa$pars, pb$pars, "mapped parameters across the cutoff")
  expect_identical(pa$ok, pb$ok)
})

test_that("the cutoff probe also holds with a wide parameter vector", {
  # The audit measures the cutoff twice: once with a small parameter vector and
  # many covariate cells, and once with a wide design matrix, where the effect
  # is far larger.  Both shapes must map identically.
  fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 3L, width = 256L)
  a <- audit_ll_args(fx)
  expect_equal(nrow(a$designs$v), 256L)
  expect_gt(ncol(fx$prop), 250L)
  b <- audit_pad_designs(a, "v", 1L)
  expect_bit_identical(audit_call_ll(fx, b), audit_call_ll(fx, a),
                       "wide design across the cutoff")
})

# --- the cell scratch budget ------------------------------------------------

test_that("the row route and the cell route are the same numbers", {
  # C7's gate.  The cell path is an optimisation: it must be indistinguishable
  # from computing every trial's coefficients one at a time.  Setting the
  # budget to zero refuses every design its cell scratch and drives the whole
  # pipeline -- mapping, transforms and bounds -- onto the general row route,
  # which is the comparison the old 256-cell cliff used to make by accident.
  on.exit(EMC2:::emc_pt_cell_budget(EMC2:::emc_pt_cell_budget()), add = TRUE)
  for (cells in c(8L, 255L, 256L, 257L, 1500L)) {
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 4L, cells = cells)
    with_cells <- audit_prologue(fx)$pars
    ll_cells <- audit_ll_direct(fx)
    expect_informative_ll(fx, ll_cells, paste("cells =", cells))

    old <- EMC2:::emc_pt_cell_budget(0)
    on.exit(EMC2:::emc_pt_cell_budget(old), add = TRUE)
    # A fresh fixture, because the design plan is decided once per table and
    # the budget is read when it is built.
    fy <- audit_fixture("RDM", n_trials = 600L, n_particles = 4L, cells = cells)
    with_rows <- audit_prologue(fy)$pars
    ll_rows <- audit_ll_direct(fy)
    EMC2:::emc_pt_cell_budget(old)

    expect_bit_identical(with_rows, with_cells, paste("mapped pars, cells =", cells))
    expect_bit_identical(ll_rows, ll_cells, paste("likelihood, cells =", cells))
  }
})

test_that("the budget is a memory limit that round-trips", {
  # It is a number of bytes, not a number of cells: what it admits depends on
  # what each cell costs, which is why nothing outside ParamTable.h has to know
  # a cell count.
  default <- EMC2:::emc_pt_cell_budget()
  on.exit(EMC2:::emc_pt_cell_budget(default), add = TRUE)
  expect_gt(default, 0)
  # Big enough that no design anyone fits is refused: three doubles per cell,
  # so the default admits hundreds of thousands of them.
  expect_gt(default / 24, 1e5)
  # The setter returns the previous value, so a caller can restore it.
  expect_identical(EMC2:::emc_pt_cell_budget(4096), default)
  expect_identical(EMC2:::emc_pt_cell_budget(default), 4096)
  expect_identical(EMC2:::emc_pt_cell_budget(), default)
  # And it refuses nonsense rather than silently taking it.
  expect_error(EMC2:::emc_pt_cell_budget(-1), "non-negative")
  expect_error(EMC2:::emc_pt_cell_budget(c(1, 2)), "one finite")
  expect_error(EMC2:::emc_pt_cell_budget(NaN), "finite")
  expect_error(EMC2:::emc_pt_cell_budget(1e300), "size_t")
  expect_identical(EMC2:::emc_pt_cell_budget(), default)
})

test_that("a design too wide for the budget still maps correctly", {
  # The budget refuses, the general route runs, and the answer is the same.
  # This is the case the old constant made unreachable without building a
  # design of a few hundred thousand cells.
  default <- EMC2:::emc_pt_cell_budget()
  on.exit(EMC2:::emc_pt_cell_budget(default), add = TRUE)
  fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 3L, cells = 300L)
  want <- audit_prologue(fx)$pars
  # Room for 100 cells at three doubles each: a 300-cell design does not fit.
  EMC2:::emc_pt_cell_budget(100 * 24)
  fy <- audit_fixture("RDM", n_trials = 600L, n_particles = 3L, cells = 300L)
  expect_bit_identical(audit_prologue(fy)$pars, want, "300 cells over a 100-cell budget")
  expect_mapper_agrees(fy, label = "300 cells over a 100-cell budget")
  EMC2:::emc_pt_cell_budget(default)
})

# --- constants --------------------------------------------------------------

test_that("duplicate constant names take the last value", {
  # `constants` is a named vector, and nothing dedupes it.  Whatever precedence
  # the mapper implements has to survive a mapper rewrite, so pin it rather
  # than leaving it to be rediscovered.
  fx <- audit_fixture("LBA", n_trials = 120L, n_particles = 3L,
                      constants = c(sv = log(1), sv = log(3)))
  pars <- audit_prologue(fx)$pars
  # The last value wins.  Compare against exp(log(3)), not 3: the round trip
  # through log is not exact for 3, and a literal here would pin a rounding
  # accident rather than the precedence rule this test is about.
  expect_true(all(pars[, "sv", ] == exp(log(3))))
  expect_false(any(pars[, "sv", ] == exp(log(1))))
  expect_mapper_agrees(fx, label = "duplicate constants")
})

test_that("a constant column is fixed while its neighbours are sampled", {
  # sv is invariant, so its bound check and mapping are hoisted out of the
  # particle loop; v is not.  Both must still be right per particle.
  fx <- audit_fixture("LBA", n_trials = 120L, n_particles = 5L,
                      formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                      constants = c(sv = log(2)))
  pars <- audit_prologue(fx)$pars
  expect_true(all(pars[, "sv", ] == 2))
  # v genuinely varies over particles, or "fixed vs sampled" is untested.
  expect_gt(length(unique(pars[1L, "v", ])), 1L)
  expect_mapper_agrees(fx, label = "fixed vs sampled")
})

# --- degenerate design shapes ----------------------------------------------

test_that("one-column and single-cell designs map identically", {
  fx <- audit_fixture("RDM", n_trials = 120L, n_particles = 4L,
                      formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1))
  for (nm in names(audit_designs(fx))) {
    expect_equal(ncol(audit_designs(fx)[[nm]]), 1L, info = nm)
  }
  expect_mapper_agrees(fx, label = "intercept only")
})

test_that("a design cell with no observed trials maps identically", {
  # design() drops coefficients for unobserved cells, so the interesting case is
  # a *level* that exists in the factor but appears in no row of this subject's
  # data.  The mapper still has to pick a representative for it.
  dat <- audit_data(200L)
  dat$E <- factor(dat$E, levels = c(levels(dat$E), "never"))
  expect_true("never" %in% levels(dat$E))
  expect_equal(sum(dat$E == "never"), 0L)
  fx <- audit_fixture("RDM", n_particles = 4L, data = dat,
                      formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1))
  expect_mapper_agrees(fx, label = "unobserved factor level")
})

test_that("mapping agrees on both the compressed and the expanded design", {
  # `.oo_expanded_designs(expand = TRUE)` takes the identity-expand branch, so
  # this is a second, structurally different route to the same numbers.
  for (cells in c(64L, 257L)) {
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 3L, cells = cells)
    expect_mapper_agrees(fx, audit_designs(fx, expand = FALSE),
                         paste("compressed, cells =", cells))
    expect_mapper_agrees(fx, audit_designs(fx, expand = TRUE),
                         paste("expanded, cells =", cells))
  }
})

# --- transforms and bounds --------------------------------------------------

test_that("split pre-sum transforms map identically across the cutoff", {
  # transform_utils.cpp carries its own 256-cell stack buffers (`ok_cell`, `v`),
  # so the transform lane has the same cliff as the mapping lane and needs the
  # same gate.  A pre-sum transform is the case that exercises it.
  for (cells in c(255L, 257L)) {
    fx <- audit_fixture("LBA", n_trials = 600L, n_particles = 3L, cells = cells,
                        pre_transform_terms = list(v = c("v", "v_x")))
    expect_mapper_agrees(fx, label = paste("split transform, cells =", cells))
  }
})

test_that("bound verdicts are per row on both sides of the cutoff", {
  # Bounds are the other consumer of the cell buffers.  Push t0 below its lower
  # bound and the whole particle must be rejected identically on both routes.
  for (cells in c(255L, 257L)) {
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 2L, cells = cells)
    fx$prop[1L, "t0"] <- log(1e-6)
    got <- audit_prologue(fx)
    expect_false(all(got$ok), info = paste("cells =", cells))
    # And the likelihood floors rather than returning something finite.
    ll <- audit_ll_direct(fx)
    expect_true(is.finite(ll[1L]))
    expect_lt(ll[1L], ll[2L])
  }
})

# --- trend fallback ---------------------------------------------------------

test_that("the trend lane still matches the reference route", {
  # A TrendRuntime disables the row-constant short cut, because the trend writes
  # into the parameter table outside the pipeline.  The audit requires the
  # general route to survive any cell work, so pin the two lanes together.
  dat <- audit_data(200L)
  # A covariate that actually varies.  With a flat covariate the trend adds
  # nothing and the comparison would hold for a trend lane that did nothing at
  # all, which is the failure mode this test is supposed to catch.
  set.seed(4)
  dat$drift <- as.numeric(scale(seq_len(nrow(dat))))
  trend <- make_trend(par_names = "B", cov_names = "drift", kernels = "lin_incr")
  fx <- audit_fixture("LBA", n_particles = 4L, data = dat,
                      formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                      covariates = "drift", trend = trend)
  expect_false(is.null(fx$model()$trend))
  pars <- audit_prologue(fx)$pars
  # B is B ~ 1, so without the trend it would be one value for every trial.
  expect_gt(length(unique(round(pars[, "B", 1L], 12))), 10L)
  expect_mapper_agrees(fx, label = "trend runtime")
})

test_that("a continuous covariate keeps the general row route honest", {
  # The refinement of a continuous covariate is bounded below by T, so there is
  # no reuse and C8 must defer to the row route.  Whatever it decides, the
  # mapped values cannot move, and v must genuinely vary by trial or the
  # fixture is not testing the row route at all.
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 3L, cells = 300L)
  pars <- audit_prologue(fx)$pars
  expect_gt(length(unique(round(pars[, "v", 1L], 12))), 100L)
  expect_mapper_agrees(fx, label = "continuous covariate")
})

# --- particle batching ------------------------------------------------------

test_that("a particle maps the same alone as inside a wide-cell batch", {
  # The planned, invariant and row-constant lanes all carry state forward from
  # the template particle.  Above the cell cutoff they carry different state,
  # so the audit's "reordered particles" gate needs checking on both sides.
  for (cells in c(200L, 400L)) {
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 5L, cells = cells)
    batch <- audit_prologue(fx)$pars
    for (i in seq_len(nrow(fx$prop))) {
      one <- fx
      one$prop <- fx$prop[i, , drop = FALSE]
      expect_bit_identical(audit_prologue(one)$pars, batch[, , i],
                           sprintf("cells = %d, particle %d", cells, i))
    }
  }
})
