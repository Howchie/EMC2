# The policy in R/failure_policy.R, pinned.
#
# The audit's point about this mechanism is that leaving it implicit is what let
# a diverging solve look like efficient sampling.  Writing the policy down only
# helps if the code follows it, and the two things that must never drift are:
# numerical rejection always continues, whatever the policy, and nothing at all
# aborts unless the user asked for it.

reset_policy_state <- function() {
  st <- EMC2:::.emc_profile_state
  st$announced <- NULL
  st$rejects <- NULL
  st$reject_msgs <- NULL
  invisible(NULL)
}

err <- function(msg) tryCatch(stop(msg), error = identity)

# --- resolving the policy ---------------------------------------------------

test_that("the policy defaults to report and refuses to be misconfigured", {
  old <- options(emc2.failure_policy = NULL)
  on.exit(options(old), add = TRUE)
  expect_identical(.emc_failure_policy(), "report")
  for (p in .EMC_FAILURE_POLICIES) {
    options(emc2.failure_policy = p)
    expect_identical(.emc_failure_policy(), p)
  }
  # A typo must not silently mean "strict": an unrecognised value falls back to
  # the default rather than to the most drastic reading of it.
  for (bad in list("Strict", "abort", "", NA_character_, 3L, c("report", "strict"))) {
    options(emc2.failure_policy = bad)
    expect_identical(.emc_failure_policy(), "report")
  }
})

# --- what aborts ------------------------------------------------------------

test_that("numerical rejection continues under every policy", {
  old <- options(emc2.failure_policy = NULL)
  on.exit(options(old), add = TRUE)
  # This is the PMwG invariant. A policy that could abort here would be a
  # policy that breaks the sampler.
  for (p in .EMC_FAILURE_POLICIES) {
    options(emc2.failure_policy = p)
    expect_identical(.emc_failure_action("numerical"), "continue")
  }
})

test_that("only strict aborts, and it aborts on everything else", {
  old <- options(emc2.failure_policy = NULL)
  on.exit(options(old), add = TRUE)
  others <- setdiff(.EMC_REJECT_CLASSES, "numerical")
  expect_identical(sort(others),
                   sort(c("infrastructure", "programming", "unknown")))
  for (p in c("silent", "report")) {
    options(emc2.failure_policy = p)
    for (cls in others) {
      expect_identical(.emc_failure_action(cls), "continue", info = paste(p, cls))
    }
  }
  options(emc2.failure_policy = "strict")
  for (cls in others) {
    expect_identical(.emc_failure_action(cls), "abort", info = cls)
  }
})

test_that("an abort carries the class it stopped for", {
  # So a handler further out, or a user reading the message, can tell a broken
  # worker from a bug in a model's R code without re-running the fit.
  cond <- tryCatch(.emc_failure_abort("programming", "particle",
                                      err("object 'foo' not found")),
                   error = identity)
  expect_s3_class(cond, "emc_failure")
  expect_s3_class(cond, "emc_failure_programming")
  expect_identical(cond$emc_class, "programming")
  expect_identical(cond$emc_source, "particle")
  expect_match(conditionMessage(cond), "object 'foo' not found")
  expect_match(conditionMessage(cond), "strict")
  expect_identical(conditionMessage(cond$emc_cause), "object 'foo' not found")
  # A condition with no message still produces a usable one.
  bare <- tryCatch(.emc_failure_abort("unknown", "gibbs", NULL), error = identity)
  expect_s3_class(bare, "emc_failure")
  expect_match(conditionMessage(bare), "no message")
})

# --- what the user sees -----------------------------------------------------

test_that("a non-numerical failure is announced once per source and class", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "report")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)

  expect_warning(.emc_failure_announce("infrastructure", "particle",
                                       err("cannot allocate vector")),
                 "infrastructure")
  # Once. A wide fit can fail thousands of times, and a warning per occurrence
  # would bury the one that matters.
  expect_silent(.emc_failure_announce("infrastructure", "particle",
                                      err("cannot allocate vector")))
  # A different class, or the same class somewhere else, is different news.
  expect_warning(.emc_failure_announce("programming", "particle",
                                       err("unused argument")), "programming")
  expect_warning(.emc_failure_announce("infrastructure", "gibbs",
                                       err("worker gave no reply")), "gibbs")
})

test_that("numerical rejection is never announced", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "report")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)
  expect_silent(.emc_failure_announce("numerical", "particle",
                                      err("chol failed")))
  options(emc2.failure_policy = "strict")
  expect_silent(.emc_failure_announce("numerical", "gibbs",
                                      err("not positive definite")))
})

test_that("silent says nothing at all", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "silent")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)
  expect_silent(.emc_failure_announce("programming", "particle",
                                      err("subscript out of bounds")))
  counts <- list(particle = c(numerical = 0L, infrastructure = 5L,
                              programming = 0L, unknown = 0L))
  expect_silent(.emc_failure_report_stage(counts, "burn"))
})

# --- the stage summary ------------------------------------------------------

test_that("the stage summary reports only what was not numerical rejection", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "report")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)

  # A stage of ordinary numerical rejection is not news, however many there
  # were: that is the sampler doing its job.
  only_numerical <- list(
    particle = c(numerical = 4000L, infrastructure = 0L, programming = 0L,
                 unknown = 0L),
    gibbs = c(numerical = 12L, infrastructure = 0L, programming = 0L,
              unknown = 0L))
  expect_null(.emc_failure_summary(only_numerical))
  expect_silent(.emc_failure_report_stage(only_numerical, "sample"))

  # One that was not is, whatever else happened alongside it.
  mixed <- only_numerical
  mixed$particle[["infrastructure"]] <- 7L
  mixed$gibbs[["programming"]] <- 2L
  summary <- .emc_failure_summary(mixed)
  expect_identical(nrow(summary), 2L)
  expect_false("numerical" %in% summary$class)
  # Ordered by how often, so the first line is the one worth chasing.
  expect_identical(summary$n, c(7L, 2L))
  expect_identical(summary$class, c("infrastructure", "programming"))
  expect_warning(.emc_failure_report_stage(mixed, "burn"),
                 "did not advance")
  expect_warning(.emc_failure_report_stage(mixed, "burn"), "burn stage")
})

test_that("the summary carries the first message of each class", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "report")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)
  # `.emc_reject_record` keeps the first message per (source, class); the
  # summary is where it becomes visible without a debugger attached.
  suppressWarnings({
    .emc_reject_record(err("cannot allocate vector of size 8 Gb"), "particle")
    .emc_reject_record(err("cannot allocate vector of size 1 Gb"), "particle")
  })
  summary <- .emc_failure_summary(.emc_reject_counts_all())
  expect_identical(nrow(summary), 1L)
  expect_identical(summary$n, 2L)
  expect_match(summary$message, "8 Gb")
})

# --- the catch sites --------------------------------------------------------

test_that("safe_new_particle repeats the previous state under report", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "report")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)
  pars <- list(alpha = matrix(c(1, 2), 2L, 1L))
  got <- NULL
  expect_warning({
    got <- with_mocked_bindings(
      safe_new_particle(s = 1L, data = NULL, pm_settings = list(x = 1),
                        prev_ll = -7, parameters = pars, stage = "burn",
                        type = "standard", tune = NULL),
      new_particle = function(...) stop("object 'foo' not found"),
      .package = "EMC2")
  }, "programming")
  # The old state, unchanged, which is what every caller expects.
  expect_identical(got$proposal, c(1, 2))
  expect_identical(got$ll, -7)
  expect_identical(got$pm_settings, list(x = 1))
  expect_identical(.emc_reject_counts("particle")[["programming"]], 1L)
})

test_that("safe_new_particle raises under strict, and only for the right class", {
  reset_policy_state()
  old <- options(emc2.failure_policy = "strict")
  on.exit({ options(old); reset_policy_state() }, add = TRUE)
  pars <- list(alpha = matrix(c(1, 2), 2L, 1L))
  call_it <- function(msg) {
    with_mocked_bindings(
      safe_new_particle(s = 1L, data = NULL, pm_settings = list(),
                        prev_ll = -7, parameters = pars, stage = "burn",
                        type = "standard", tune = NULL),
      new_particle = function(...) stop(msg), .package = "EMC2")
  }
  cond <- suppressWarnings(tryCatch(call_it("cannot allocate vector"),
                                    error = identity))
  expect_s3_class(cond, "emc_failure_infrastructure")

  # A proposal the model has no density for is still rejected, not raised --
  # under strict as under anything else.
  got <- call_it("chol(): matrix is not positive definite")
  expect_identical(got$proposal, c(1, 2))
  expect_identical(got$ll, -7)
})
