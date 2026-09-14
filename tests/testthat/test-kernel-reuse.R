# C9's first deliverable: measure before specialising.
#
# The audit names three candidates for kernel-preparation reuse -- BAwL's
# `natural_normalizer`, RDM's scale/geometry, the DDM's `a/s`, `v/s`, `sv/s`
# scaling -- and says to instrument them first, then take one model at a time
# and abandon any whose setup outweighs the saved arithmetic. That decision
# needs a number, and the number is how often the arithmetic is recomputed
# against how often it changes.

test_that("each model declares the columns its subexpression reads", {
  # A declaration, not a heuristic: it says which columns the candidate
  # arithmetic is a function of, and nothing else.
  expect_identical(EMC2:::emc_kernel_reuse_columns("RDM"), c("v", "B", "A", "s"))
  # `t0` is deliberately absent: it enters only through `rt - t0`, which varies
  # per trial whatever the design says, so it is not reusable arithmetic.
  expect_false("t0" %in% EMC2:::emc_kernel_reuse_columns("RDM"))
  # natural_normalizer is a pnorm of v/sv and nothing else.
  expect_identical(EMC2:::emc_kernel_reuse_columns("LBA"), c("v", "sv"))
  expect_identical(EMC2:::emc_kernel_reuse_columns("BAwL"), c("v", "sv"))
  expect_identical(EMC2:::emc_kernel_reuse_columns("DDM"), c("a", "v", "sv", "s"))
  # Variants inherit their family's declaration.
  expect_identical(EMC2:::emc_kernel_reuse_columns("RDMSWTN"), c("v", "B", "A", "s"))
  # A model with no declared candidate is simply not measured.
  expect_identical(length(EMC2:::emc_kernel_reuse_columns("LNR")), 0L)
  expect_identical(length(EMC2:::emc_kernel_reuse_columns("nonsense")), 0L)
})

test_that("the flag round-trips and refuses nonsense", {
  default <- EMC2:::emc_kernel_stats()
  on.exit(EMC2:::emc_kernel_stats(default), add = TRUE)
  # Off by default: when it is on the caller pays a clock read per particle.
  expect_false(default)
  expect_identical(EMC2:::emc_kernel_stats(TRUE), FALSE)
  expect_identical(EMC2:::emc_kernel_stats(), TRUE)
  expect_identical(EMC2:::emc_kernel_stats(FALSE), TRUE)
  expect_error(EMC2:::emc_kernel_stats(NA), "TRUE or FALSE")
  expect_error(EMC2:::emc_kernel_stats(c(TRUE, TRUE)), "TRUE or FALSE")
})

test_that("nothing is counted while the flag is off", {
  on.exit({ EMC2:::emc_kernel_stats(FALSE); EMC2:::emc_kernel_stats_reset() },
          add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(FALSE)
  fx <- audit_fixture("RDM", n_trials = 300L, n_particles = 5L, cells = 4L)
  invisible(audit_ll_direct(fx))
  expect_identical(nrow(EMC2:::emc_kernel_stats_read()), 0L)
})

test_that("the counters measure recomputation against change", {
  on.exit({ EMC2:::emc_kernel_stats(FALSE); EMC2:::emc_kernel_stats_reset() },
          add = TRUE)
  for (model in c("RDM", "LBA", "DDM")) {
    EMC2:::emc_kernel_stats_reset()
    EMC2:::emc_kernel_stats(TRUE)
    fx <- audit_fixture(model, n_trials = 600L, n_particles = 10L, cells = 4L)
    invisible(audit_ll_direct(fx))
    EMC2:::emc_kernel_stats(FALSE)
    st <- EMC2:::emc_kernel_stats_read()
    expect_identical(nrow(st), 1L, info = model)
    expect_identical(st$calls[1L], 10, info = model)
    # Rows are per-trial evaluations and cells are what a cell-resolution
    # version would do, so reuse is rows/cells and can never be below 1.
    expect_gt(st$rows[1L], 0)
    expect_gt(st$cells[1L], 0)
    expect_gte(st$rows[1L] / st$cells[1L], 1, label = paste(model, "reuse"))
    expect_gte(st$seconds[1L], 0)
  }
})

test_that("reuse falls as the design gets wider, and that is the point", {
  # A narrow design recomputes the same handful of values for every trial; a
  # wide one has more distinct values to compute. If the ratio did not fall
  # with the cell count the measurement would not be measuring anything.
  on.exit({ EMC2:::emc_kernel_stats(FALSE); EMC2:::emc_kernel_stats_reset() },
          add = TRUE)
  ratio <- function(cells) {
    EMC2:::emc_kernel_stats_reset()
    EMC2:::emc_kernel_stats(TRUE)
    fx <- audit_fixture("RDM", n_trials = 600L, n_particles = 6L, cells = cells)
    invisible(audit_ll_direct(fx))
    EMC2:::emc_kernel_stats(FALSE)
    st <- EMC2:::emc_kernel_stats_read()
    st$rows[1L] / st$cells[1L]
  }
  narrow <- ratio(4L)
  wide <- ratio(64L)
  expect_gt(narrow, wide)
  # And the narrow case has enough reuse to be worth a commit's work at all,
  # which is the question C9 exists to answer.
  expect_gt(narrow, 10)
})

test_that("instrumentation changes no likelihood", {
  # It reads a clock and a cached partition; it must not touch a number.
  on.exit({ EMC2:::emc_kernel_stats(FALSE); EMC2:::emc_kernel_stats_reset() },
          add = TRUE)
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 400L, n_particles = 5L, cells = 8L)
    EMC2:::emc_kernel_stats(FALSE)
    off <- audit_ll_direct(fx)
    EMC2:::emc_kernel_stats(TRUE)
    on <- audit_ll_direct(fx)
    EMC2:::emc_kernel_stats(FALSE)
    expect_bit_identical(on, off, paste(model, "- instrumentation on vs off"))
  }
  EMC2:::emc_kernel_stats_reset()
})

test_that("the counters reset", {
  on.exit({ EMC2:::emc_kernel_stats(FALSE); EMC2:::emc_kernel_stats_reset() },
          add = TRUE)
  EMC2:::emc_kernel_stats_reset()
  EMC2:::emc_kernel_stats(TRUE)
  fx <- audit_fixture("RDM", n_trials = 200L, n_particles = 3L, cells = 4L)
  invisible(audit_ll_direct(fx))
  expect_gt(nrow(EMC2:::emc_kernel_stats_read()), 0L)
  EMC2:::emc_kernel_stats_reset()
  expect_identical(nrow(EMC2:::emc_kernel_stats_read()), 0L)
  EMC2:::emc_kernel_stats(FALSE)
})

test_that("the profile schema carries the kernel columns", {
  sch <- EMC2:::.emc_profile_schema
  for (nm in c("kernel_rows", "kernel_cells", "kernel_seconds")) {
    expect_true(nm %in% names(sch), info = nm)
    expect_identical(sch[[nm]]$group, "kernel", info = nm)
  }
  # A run that never switched instrumentation on carries no kernel numbers, not
  # zeros -- zero reuse and unmeasured are different claims.
  EMC2:::emc_kernel_stats_reset()
  expect_true(all(is.na(EMC2:::.emc_kernel_totals())))
  row <- EMC2:::.emc_profile_row(iteration = 1L, total = 0.1)
  expect_true(is.na(row$kernel_rows))
})

test_that("the profile report identifies the execution route", {
  expect_true("route" %in% names(EMC2:::.emc_profile_schema))
  row <- EMC2:::.emc_profile_row(iteration = 1L, route = "serial", total = 0.1)
  report <- capture.output(
    EMC2:::.emc_profile_report(row, elapsed = 0.1, drop_first = FALSE)
  )
  expect_true(any(grepl("route serial", report, fixed = TRUE)))
})
