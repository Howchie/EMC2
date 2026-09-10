# The profiling record is a contract between the worker pool, run_stage() and
# two benchmark scripts.  These tests pin the parts of that contract that broke
# silently before it was declared in one place: a field name that no longer
# means what a reader thinks it means, and a metric that understates worker
# cost by the number of subjects per worker.

test_that("the schema is internally consistent", {
  sch <- .emc_profile_schema
  expect_true(length(sch) > 0)
  # The list is keyed by field name, and the report walks it by name.
  expect_identical(names(sch), unname(vapply(sch, `[[`, character(1), "name")))
  expect_false(any(duplicated(names(sch))))
  for (f in sch) {
    expect_true(f$type %in% c("numeric", "integer", "character"), info = f$name)
    expect_true(f$group %in% c("id", "time", "work", "wire", "memory", "health"),
                info = f$name)
    expect_true(nzchar(f$desc), info = f$name)
  }
  # Every declared parent must be a real field, or the report's percentages
  # would be taken against a column that does not exist.
  parents <- stats::na.omit(vapply(sch, function(f) f$parent, character(1)))
  expect_true(all(parents %in% names(sch)))
  # Only time fields nest.
  nested <- Filter(function(f) !is.na(f$parent), sch)
  expect_true(all(vapply(nested, function(f) f$group == "time", logical(1))))
})

test_that("rows carry every field with the declared type", {
  row <- .emc_profile_row()
  expect_identical(names(row), names(.emc_profile_schema))
  expect_identical(nrow(row), 1L)
  for (f in .emc_profile_schema) {
    expect_identical(class(row[[f$name]]),
                     switch(f$type, integer = "integer",
                            character = "character", "numeric"),
                     info = f$name)
    expect_true(is.na(row[[f$name]]), info = f$name)
  }
  filled <- .emc_profile_row(iteration = 4, total = 0.25, worker_cpu_max = 0.2)
  expect_identical(filled$iteration, 4L)
  expect_equal(filled$total, 0.25)
  # Unsupplied fields stay NA rather than 0: a missing measurement and a
  # measurement of zero are different claims.
  expect_true(is.na(filled$subject_cpu_max))
})

test_that("an undeclared field is an error, not an extra column", {
  # This is the whole point of the schema.  `worker_max` used to be written by
  # the pool, copied by run_stage() and printed by a benchmark under a label
  # that did not match it; re-introducing it by typo must not be silent.
  expect_error(.emc_profile_row(worker_max = 1), "unknown profile field")
  expect_error(.emc_profile_row(worker_sum = 1), "unknown profile field")
  expect_error(.emc_profile_row(1), "must be named")
})

test_that("binding rows keeps the schema and tolerates empties", {
  rows <- list(.emc_profile_row(iteration = 1), NULL, .emc_profile_row(iteration = 2))
  bound <- .emc_profile_bind(rows)
  expect_identical(nrow(bound), 2L)
  expect_identical(names(bound), names(.emc_profile_schema))
  expect_identical(bound$iteration, c(1L, 2L))
  empty <- .emc_profile_bind(list(NULL))
  expect_identical(nrow(empty), 0L)
  expect_identical(names(empty), names(.emc_profile_schema))
})

test_that("worker workload is the assigned partition, not the largest subject", {
  # The measurement contract this commit repairs.  Eight subjects of 10 ms
  # each, two workers: the slowest single subject is 10 ms, but each worker was
  # actually handed 40 ms of work.  Reading the former as "the slowest worker"
  # makes 30 ms of real work look like transport overhead.
  times <- rep(0.01, 8)
  part <- list(1:4, 5:8)
  worker_load <- vapply(part, function(subs) sum(times[subs]), numeric(1))
  expect_equal(max(times), 0.01)
  expect_equal(max(worker_load), 0.04)
  expect_equal(sum(worker_load), sum(times))
})

test_that("range helpers return NA rather than an infinity", {
  # max(numeric(0)) is -Inf with a warning, and an -Inf in a profile column
  # poisons every mean taken over it afterwards.
  expect_true(is.na(.emc_range_or_na(numeric(0), max)))
  expect_true(is.na(.emc_range_or_na(c(NA_real_, NA_real_), min)))
  expect_equal(.emc_range_or_na(c(NA, 2, 4), max), 4)
  expect_equal(.emc_range_or_na(c(NA, 2, 4), mean), 3)
})

test_that("failures are classified without changing rejection behaviour", {
  expect_identical(.emc_classify_failure(
    simpleError("the leading minor of order 2 is not positive definite")),
    "numerical")
  expect_identical(.emc_classify_failure(simpleError("cannot open the connection")),
                   "infrastructure")
  expect_identical(.emc_classify_failure(simpleError("could not find function \"f\"")),
                   "programming")
  expect_identical(.emc_classify_failure(simpleError("something unfamiliar")),
                   "unknown")
  expect_identical(.emc_classify_failure(NULL), "unknown")
  # An allocation failure inside a numerical routine is still infrastructure:
  # calling it numerical rejection is what would let a broken pool masquerade
  # as a badly mixing chain.
  expect_identical(.emc_classify_failure(
    simpleError("cannot allocate vector of size 4.0 Gb")), "infrastructure")
})

test_that("rejection counters accumulate per source and difference correctly", {
  old <- .emc_profile_state$rejects
  on.exit(.emc_profile_state$rejects <- old, add = TRUE)
  .emc_reject_reset()
  expect_true(all(.emc_reject_counts("particle") == 0L))

  before <- .emc_reject_counts("particle")
  .emc_reject_record(simpleError("chol failed"), "particle")
  .emc_reject_record(simpleError("chol failed"), "particle")
  .emc_reject_record(simpleError("could not find function \"g\""), "gibbs")
  expect_identical(unname(.emc_reject_counts("particle")["numerical"]), 2L)
  # Sources are independent: a Gibbs failure must not inflate the particle
  # count, since run_stage() reads only one of them and the workers report the
  # other.
  expect_identical(unname(.emc_reject_counts("particle")["programming"]), 0L)
  expect_identical(unname(.emc_reject_counts("gibbs")["programming"]), 1L)

  delta <- .emc_reject_delta(before, "particle")
  expect_identical(unname(delta["numerical"]), 2L)
  fields <- .emc_reject_row_fields(delta, "particle")
  expect_true("reject_particle_nume" %in% names(fields))
  expect_true(all(names(fields) %in% names(.emc_profile_schema)))
  row <- do.call(.emc_profile_row, fields)
  expect_identical(row$reject_particle_nume, 2L)
  .emc_reject_reset()
})

test_that("the reporter prints the corrected metric and survives empty input", {
  out <- capture.output(.emc_profile_report(NULL))
  expect_match(paste(out, collapse = " "), "no profile recorded")
  out <- capture.output(.emc_profile_report(.emc_profile_bind(list())))
  expect_match(paste(out, collapse = " "), "no profile recorded")

  rows <- lapply(1:3, function(i) .emc_profile_row(
    iteration = i, stage = "preburn", subjects = 8L, workers = 2L,
    workers_active = 2L, total = 0.04, gibbs = 0.001, particle = 0.035,
    response_wait = 0.034, fill = 0.001, subject_cpu_max = 0.01,
    subject_cpu_sum = 0.08, worker_cpu_max = 0.04, particles_sum = 800,
    particles_max = 100))
  out <- paste(capture.output(
    .emc_profile_report(.emc_profile_bind(rows), elapsed = 1.2,
                        drop_first = FALSE)), collapse = "\n")
  expect_match(out, "worker_cpu_max")
  expect_match(out, "slowest worker's assigned partition")
  expect_match(out, "subject_cpu_max")
  # The retired label must not come back through the formatter either.
  expect_false(grepl("worker_max", out, fixed = TRUE))
  # Nested time fields are reported against their own parent, so
  # response_wait's share is of `particle`, not of `total`.
  expect_match(out, "response_wait")
  # A health section only appears when something actually went wrong.
  expect_false(grepl("health", out, fixed = TRUE))

  rows[[2]] <- .emc_profile_row(iteration = 2, total = 0.04, fallback_subjects = 3,
                                degraded = 1)
  out <- paste(capture.output(
    .emc_profile_report(.emc_profile_bind(rows), drop_first = FALSE)),
    collapse = "\n")
  expect_match(out, "health")
  expect_match(out, "fallback_subjects")
})

test_that("expensive measurements stay behind their own switch", {
  op <- options(emc2.sampler_profile = TRUE, emc2.sampler_profile_expensive = FALSE)
  on.exit(options(op), add = TRUE)
  expect_true(.emc_profile_enabled())
  expect_false(.emc_profile_expensive())
  options(emc2.sampler_profile = FALSE, emc2.sampler_profile_expensive = TRUE)
  # Expensive implies enabled: it is a refinement of profiling, not a way to
  # pay for /proc reads in a run that records nothing.
  expect_false(.emc_profile_expensive())
})

test_that("tree memory is measured, and labelled with how", {
  mem <- .emc_profile_memory()
  expect_true(is.list(mem) && all(c("bytes", "method") %in% names(mem)))
  if (file.exists(file.path("/proc", Sys.getpid(), "smaps_rollup"))) {
    # PSS divides shared pages by the number of processes mapping them; summing
    # RSS over a forked tree counts the template's context once per worker.
    expect_identical(mem$method, "pss")
    expect_true(mem$bytes > 0)
    # Passing the current process explicitly is the same tree as the default.
    # Not bit-identical: R allocates between the two reads, so compare loosely.
    both <- .emc_profile_memory(Sys.getpid())
    expect_equal(both$bytes, mem$bytes, tolerance = 0.05)
    # Duplicate pids must not be counted twice.
    dup <- .emc_profile_memory(rep(Sys.getpid(), 4L))
    expect_equal(dup$bytes, mem$bytes, tolerance = 0.05)
  } else {
    expect_true(is.na(mem$bytes) || mem$bytes > 0)
  }
})
