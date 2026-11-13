local_rng_guard()

test_that("rebuild claims preserve order, NULL values and subject streams", {
  skip_on_os("windows")
  set.seed(71)
  seeds <- .emc_subject_streams(31L)
  X <- stats::setNames(seq_along(seeds), paste0("subject", seq_along(seeds)))
  work <- function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    if (i == 7L) return(NULL)
    list(draw = rnorm(4), seed = .Random.seed)
  }
  expected <- .emc_with_preserved_rng(lapply(X, work))
  before <- .Random.seed
  for (i in seq_len(12L)) {
    got <- .emc_with_preserved_rng(.emc_rebuild_lapply(X, work, 2L, rev(seq_along(X))))
    expect_identical(got, expected)
    expect_identical(.Random.seed, before)
    expect_false(.emc_rebuild_queue$busy)
  }
})

test_that("rebuild queues propagate worker errors and can be reused", {
  skip_on_os("windows")
  expect_error(suppressWarnings(.emc_rebuild_lapply(1:12, function(i) {
    if (i == 3L) stop("rebuild failed")
    i
  }, 2L)), "rebuild failed")
  expect_false(.emc_rebuild_queue$busy)
  expect_identical(.emc_rebuild_lapply(1:12, identity, 2L), as.list(1:12))
})

test_that("single-core and short rebuilds keep the serial path", {
  expect_identical(.emc_rebuild_lapply(integer(), identity, 1L), list())
  expect_identical(.emc_rebuild_lapply(1:3, identity, 1L), as.list(1:3))
})

test_that("measured rebuild costs select serial only when it has a margin", {
  old_costs <- .emc_rebuild_queue$costs
  old_decision <- .emc_rebuild_queue$last_decision
  on.exit({
    .emc_rebuild_queue$costs <- old_costs
    .emc_rebuild_queue$last_decision <- old_decision
  }, add = TRUE)
  .emc_rebuild_queue$costs <- list()
  cheap_id <- .emc_rebuild_cost_id("dispatch_test", 2L)
  expensive_id <- .emc_rebuild_cost_id("dispatch_expensive", 2L)
  .emc_rebuild_queue$costs[[cheap_id]] <- list(
    job_seconds = rep(.001, 8L), startup_seconds = .05)
  .emc_rebuild_queue$costs[[expensive_id]] <- list(
    job_seconds = rep(.2, 8L), startup_seconds = .005)

  cheap <- .emc_rebuild_prediction("dispatch_test", 8L, 2L)
  expensive <- .emc_rebuild_prediction("dispatch_expensive", 8L, 2L)
  expect_lt(cheap$serial, .9 * cheap$parallel)
  expect_gt(expensive$serial, .9 * expensive$parallel)

  local_mocked_bindings(
    auto_mclapply = function(...) stop("unexpected parallel dispatch"),
    .package = "EMC2")
  seen <- integer()
  X <- stats::setNames(as.list(seq_len(8L)), paste0("job", seq_len(8L)))
  set.seed(108)
  before <- .Random.seed
  result <- .emc_rebuild_lapply(X, function(i) {
    seen <<- c(seen, i)
    runif(1)
    i
  }, 2L, order = 8:1, cost_key = "dispatch_test")
  expect_identical(result, X)
  expect_identical(seen, 8:1)
  expect_identical(.Random.seed, before)
  expect_identical(.emc_rebuild_queue$last_decision$route, "serial")
})

test_that("serial rebuild errors propagate and restore the caller RNG", {
  old_costs <- .emc_rebuild_queue$costs
  on.exit(.emc_rebuild_queue$costs <- old_costs, add = TRUE)
  .emc_rebuild_queue$costs <- list()
  .emc_rebuild_queue$costs[[.emc_rebuild_cost_id("dispatch_error", 2L)]] <-
    list(job_seconds = rep(.001, 8L), startup_seconds = .05)
  set.seed(109)
  before <- .Random.seed
  expect_error(.emc_rebuild_lapply(1:8, function(i) {
    runif(1)
    if (i == 3L) stop("serial rebuild failed")
    i
  }, 2L, cost_key = "dispatch_error"), "serial rebuild failed")
  expect_identical(.Random.seed, before)
  expect_false(isTRUE(.emc_rebuild_queue$busy))
})

test_that("parallel rebuild timing calibrates job and startup costs", {
  skip_on_os("windows")
  old_costs <- .emc_rebuild_queue$costs
  old_decision <- .emc_rebuild_queue$last_decision
  on.exit({
    .emc_rebuild_queue$costs <- old_costs
    .emc_rebuild_queue$last_decision <- old_decision
  }, add = TRUE)
  .emc_rebuild_queue$costs <- list()
  result <- .emc_rebuild_lapply(1:12, function(i) {
    Sys.sleep(.002)
    i
  }, 2L, order = 12:1, cost_key = "dispatch_calibration")
  expect_identical(result, as.list(1:12))
  model <- .emc_rebuild_queue$costs[[
    .emc_rebuild_cost_id("dispatch_calibration", 2L)]]
  expect_length(model$job_seconds, 12L)
  expect_true(all(is.finite(model$job_seconds) & model$job_seconds >= 0))
  expect_true(is.finite(model$startup_seconds) && model$startup_seconds >= 0)
  expect_identical(.emc_rebuild_queue$last_decision$route, "parallel")
})


test_that("direct array slices match the split index implementation", {
  previous <- function(a, n) setNames(lapply(
    split(a, arrayInd(seq_along(a), dim(a))[, n]),
    array, dim = dim(a)[-n], dimnames(a)[-n]), dimnames(a)[[n]])
  for (dims in list(c(3L, 4L), c(2L, 3L, 5L), c(1L, 4L, 2L), c(2L, 1L, 3L, 4L))) {
    a <- array(seq_len(prod(dims)), dims)
    for (named in c(FALSE, TRUE)) {
      if (named) dimnames(a) <- lapply(dims, function(n) paste0("x", seq_len(n)))
      for (n in seq_along(dims)) expect_identical(split_along_dim(a, n), previous(a, n))
    }
  }
  a <- array(c(NA_real_, NaN, Inf, -Inf, -0, 1), c(1L, 2L, 3L))
  for (n in 1:3) expect_identical(split_along_dim(a, n), previous(a, n))
})

test_that("a job timed in the calling process overrides inflated fork timings", {
  skip_on_os("windows")
  old_costs <- .emc_rebuild_queue$costs
  old_decision <- .emc_rebuild_queue$last_decision
  on.exit({
    .emc_rebuild_queue$costs <- old_costs
    .emc_rebuild_queue$last_decision <- old_decision
  }, add = TRUE)
  .emc_rebuild_queue$costs <- list()
  # Fork timings that favour forking.
  seed_costs <- function(key) .emc_rebuild_queue$costs[[.emc_rebuild_cost_id(key, 2L)]] <-
    list(job_seconds = rep(.01, 12L), startup_seconds = .005)
  seed_costs("probe_cheap")
  seen <- integer()
  set.seed(110)
  before <- .Random.seed
  result <- .emc_rebuild_lapply(1:12, function(i) {
    seen <<- c(seen, i)
    runif(1)
    i
  }, 2L, order = 12:1, cost_key = "probe_cheap")
  expect_identical(result, as.list(1:12))
  # Every job ran once, here, in the requested order.
  expect_identical(seen, 12:1)
  expect_identical(.Random.seed, before)
  expect_identical(.emc_rebuild_queue$last_decision$route, "serial")

  # Jobs that are slow here too are still forked; the timed job is not rerun.
  seed_costs("probe_slow")
  seen <- integer()
  result <- .emc_rebuild_lapply(1:12, function(i) {
    seen <<- c(seen, i)
    Sys.sleep(.02)
    i
  }, 2L, order = 12:1, cost_key = "probe_slow")
  expect_identical(result, as.list(1:12))
  expect_identical(seen, 12L)
  expect_identical(.emc_rebuild_queue$last_decision$route, "parallel")
  # The stored timing settles later rounds without another job run here.
  seen <- integer()
  result <- .emc_rebuild_lapply(1:12, function(i) {
    seen <<- c(seen, i)
    Sys.sleep(.02)
    i
  }, 2L, order = 12:1, cost_key = "probe_slow")
  expect_identical(result, as.list(1:12))
  expect_length(seen, 0L)
})
