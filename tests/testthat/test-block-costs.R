# C16's first half: the per-block costs, which no per-iteration record ever saw.
#
# `run_emc` builds chain and efficient proposals, runs convergence diagnostics,
# writes a checkpoint and joins the block into the fit -- all between blocks, all
# invisible to the sampler profile. They are C1 schema fields now.

test_that("the schema declares the block fields the audit reserved", {
  sch <- EMC2:::.emc_profile_schema
  for (nm in c("block_proposals", "block_eff_proposals", "block_check_progress",
               "block_save", "block_concat")) {
    expect_true(nm %in% names(sch), info = nm)
    expect_identical(sch[[nm]]$group, "block", info = nm)
    # Each is a component of the block, not of an iteration: a block is not an
    # iteration and nesting them under `total` would take percentages against
    # the wrong denominator.
    expect_identical(sch[[nm]]$parent, "block_total", info = nm)
  }
  expect_identical(sch[["block_total"]]$group, "block")
  expect_true(is.na(sch[["block_total"]]$parent))
})

test_that("nothing is recorded unless profiling is on", {
  old <- options(emc2.sampler_profile = FALSE)
  on.exit({ options(old); EMC2:::.emc_block_reset() }, add = TRUE)
  EMC2:::.emc_block_reset()
  EMC2:::.emc_block_record(block = 1L, stage = "burn", block_total = 1)
  expect_null(EMC2:::.emc_block_profile())
})

test_that("records accumulate and reset", {
  old <- options(emc2.sampler_profile = TRUE)
  on.exit({ options(old); EMC2:::.emc_block_reset() }, add = TRUE)
  EMC2:::.emc_block_reset()
  expect_null(EMC2:::.emc_block_profile())
  for (b in 1:3) {
    EMC2:::.emc_block_record(block = b, stage = "burn", block_iterations = 10L,
                             block_chains = 2L, block_total = 0.1 * b,
                             block_sample = 0.09 * b)
  }
  got <- EMC2:::.emc_block_profile()
  expect_identical(nrow(got), 3L)
  expect_identical(got$block, 1:3)
  expect_equal(got$block_total, c(0.1, 0.2, 0.3))
  EMC2:::.emc_block_reset()
  expect_null(EMC2:::.emc_block_profile())
})

test_that("an unknown block field is refused like any other", {
  old <- options(emc2.sampler_profile = TRUE)
  on.exit({ options(old); EMC2:::.emc_block_reset() }, add = TRUE)
  expect_error(EMC2:::.emc_block_record(block = 1L, block_concatenate = 1),
               "unknown profile field")
})

test_that("blocks report even with no per-iteration record", {
  # Which is every multi-block run: the per-iteration record is built inside a
  # forked chain and does not outlive `concat_emc`, while the block record is
  # kept in the master. Having one without the other is ordinary, not an error.
  old <- options(emc2.sampler_profile = TRUE)
  on.exit({ options(old); EMC2:::.emc_block_reset() }, add = TRUE)
  EMC2:::.emc_block_reset()
  EMC2:::.emc_block_record(block = 1L, stage = "burn", block_total = 2,
                           block_sample = 1.8, block_check_progress = 0.1)
  blocks <- EMC2:::.emc_block_profile()
  out <- capture.output(EMC2:::.emc_profile_report(NULL, elapsed = 4,
                                                   blocks = blocks))
  txt <- paste(out, collapse = " ")
  expect_match(txt, "per-block costs")
  expect_match(txt, "block_sample")
  expect_false(grepl("no profile recorded", txt))
  # And with neither, it still says so rather than printing an empty section.
  bare <- capture.output(EMC2:::.emc_profile_report(NULL))
  expect_match(paste(bare, collapse = " "), "no profile recorded")
})

test_that("the proposal split reports only what was measured", {
  # preburn builds neither kind of proposal, so both stay NA rather than
  # reading as "measured, and it took no time".
  old <- options(emc2.sampler_profile = TRUE)
  # A local handle: `EMC2:::x$y <- v` is a replacement call on `EMC2`, and R
  # has no `:::<-`.
  st <- EMC2:::.emc_profile_state
  on.exit({ options(old); st$proposal_split <- NULL }, add = TRUE)
  st$proposal_split <- NULL
  expect_true(is.na(EMC2:::.emc_proposal_split("chain", 1)))
  expect_true(is.na(EMC2:::.emc_proposal_split("eff", 1)))
  st$proposal_split <- c(chain = 0.25, eff = NA_real_)
  expect_identical(EMC2:::.emc_proposal_split("chain", 1), 0.25)
  expect_true(is.na(EMC2:::.emc_proposal_split("eff", 1)))
})
