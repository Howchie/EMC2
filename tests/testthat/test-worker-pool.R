# Persistent per-block worker pool for the particle step (R/chain_pool.R).

test_that("LPT partitioning balances by cost, not by count", {
  part <- EMC2:::.emc_lpt_partition(c(10, 1, 1, 1, 1, 1, 1, 1, 1, 1), 2)
  expect_length(part, 2)
  expect_setequal(unlist(part), 1:10)
  loads <- vapply(part, function(p) sum(c(10, rep(1, 9))[p]), numeric(1))
  # A contiguous equal-count split would be 14 vs 5; LPT gets within one unit.
  expect_lte(diff(range(loads)), 1)

  # Every worker gets a slot even when there is less work than workers.
  part <- EMC2:::.emc_lpt_partition(c(3, 1), 4)
  expect_length(part, 4)
  expect_setequal(unlist(part), 1:2)
  expect_true(all(vapply(part, is.numeric, logical(1))))
})

test_that("one worker takes everything", {
  expect_equal(EMC2:::.emc_lpt_partition(c(5, 2, 9), 1), list(1:3))
})

test_that("subject cost proxies the rows the likelihood loops over", {
  d <- list(data.frame(a = 1:10), data.frame(a = 1:3))
  expect_equal(EMC2:::.emc_subject_cost(d), c(10, 3))
  # Joint models carry a list of dadms per subject.
  dj <- list(list(data.frame(a = 1:4), data.frame(a = 1:6)))
  expect_equal(EMC2:::.emc_subject_cost(dj), 10)
})

test_that("each subject gets an independent stream and the master moves past them", {
  RNGkind("L'Ecuyer-CMRG"); set.seed(3)
  st <- EMC2:::.emc_subject_streams(4)
  expect_length(st, 4)
  expect_equal(length(unique(vapply(st, function(s) paste(s, collapse = ","),
                                    character(1)))), 4)
  master <- get(".Random.seed", envir = globalenv())
  expect_false(any(vapply(st, identical, logical(1), master)))
})

test_that("the pool is refused where it cannot work", {
  expect_null(EMC2:::.emc_wpool_start(1, list()))
  expect_null(EMC2:::.emc_wpool_start(0, list()))
})

test_that("workers survive a block, compute, and shut down", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  # A stand-in context: the pool machinery does not care what the work is, only
  # that .emc_wpool_compute() can do it, so exercise the transport directly.
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not fork a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)
  expect_equal(pool$n, 2)
  expect_true(dir.exists(pool$dir))

  # Round trip twice, to show the workers persist rather than being one-shot.
  for (k in 1:2) {
    for (w in 1:2) {
      serialize(list(subs = w, echo = k), pool$wcs[[w]]); flush(pool$wcs[[w]])
    }
    got <- lapply(1:2, function(w) unserialize(pool$rcs[[w]]))
    # No real ctx, so compute() errors and the worker reports it rather than dying.
    expect_true(all(vapply(got, function(g) !is.null(g$failed), logical(1))))
  }
  EMC2:::.emc_wpool_stop(pool)
  expect_false(dir.exists(pool$dir))
  on.exit(NULL)
})

test_that("a pooled fit runs, and repeats itself exactly", {
  skip_on_os("windows")
  skip_on_cran()
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))
  fit_once <- function() {
    RNGkind("L'Ecuyer-CMRG"); set.seed(42)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 2, compress = TRUE))
    suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 2,
      cores_per_chain = 2, step_size = 10, max_tries = 1, verbose = FALSE,
      stop_criteria = list(iter = 10, max_gd = Inf, min_unique = 0, min_es = 0)))
  }
  a <- fit_once()
  expect_equal(a[[1]]$samples$idx, 11)
  expect_true(all(is.finite(a[[1]]$samples$subj_ll[, 11])))
  # Per-subject RNG streams make the pooled path independent of scheduling,
  # which the default mc.cores-derived path is not.
  expect_identical(a[[1]]$samples$alpha, fit_once()[[1]]$samples$alpha)
})

test_that("a fit does not move when the core count does", {
  skip_on_os("windows")
  skip_on_cran()
  # The invariant the whole design rests on, and what makes it safe for a chain
  # to grow its pool onto a finished sibling's cores mid-block.  Covers the
  # start points too: those used to be seeded from `mc.cores`.
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))
  fit_on <- function(cores_per_chain) {
    RNGkind("L'Ecuyer-CMRG"); set.seed(11)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 2, compress = TRUE))
    suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 2,
      cores_per_chain = cores_per_chain, step_size = 8, max_tries = 1,
      verbose = FALSE,
      stop_criteria = list(iter = 8, max_gd = Inf, min_unique = 0, min_es = 0)))
  }
  one <- fit_on(1); three <- fit_on(3)
  expect_identical(one[[1]]$samples$alpha, three[[1]]$samples$alpha)
  expect_identical(one[[2]]$samples$alpha, three[[2]]$samples$alpha)
})

test_that("messages larger than the pipe buffer still round-trip", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  # A group covariance is n_pars x n_pars, so the per-iteration message passes
  # the 64 kB pipe buffer at about 35 parameters.  A non-blocking writer fails
  # outright there and costs the block its parallelism.
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not fork a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)

  big <- list(subs = 1L, group_chol = matrix(0, 400, 400))
  expect_gt(length(serialize(big, NULL)), 65536L)
  expect_silent({ serialize(big, pool$wcs[[1]]); flush(pool$wcs[[1]]) })
  # No real ctx, so the worker reports an error -- but it received the whole
  # message, which is the point.
  expect_false(is.null(unserialize(pool$rcs[[1]])$failed))
})

test_that("fitting leaves the caller's generator as it found it", {
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  # What the caller's generator looks like after exactly one draw off it.
  RNGkind("Mersenne-Twister"); set.seed(5)
  invisible(sample.int(.Machine$integer.max, 1L))
  expected <- get(".Random.seed", envir = globalenv())

  RNGkind("Mersenne-Twister"); set.seed(5)
  st <- EMC2:::.emc_subject_streams(3)

  expect_length(st, 3)
  expect_identical(RNGkind()[1], "Mersenne-Twister")   # not left on L'Ecuyer
  expect_identical(get(".Random.seed", envir = globalenv()), expected)
  # The streams themselves are still L'Ecuyer, and still distinct.
  expect_length(unique(vapply(st, paste, character(1), collapse = ",")), 3L)
})

test_that("a pool marked dead still returns every subject's result", {
  skip_on_os("windows")
  # alive = FALSE means nothing is sent; the master must compute the lot itself
  # so that a broken pool costs speed and not an iteration.
  calls <- new.env(parent = emptyenv()); calls$n <- 0L
  fake <- function(msg, ctx) {
    calls$n <- calls$n + 1L
    list(props = matrix(msg$subs, ctx$n_pars + 1L, length(msg$subs)),
         pm = as.list(msg$subs), seeds = as.list(msg$subs),
         times = rep(1, length(msg$subs)))
  }
  pool <- list(n = 2, alive = FALSE, wcs = list(), rcs = list())
  ctx <- list(n_pars = 2)
  part <- list(c(1L, 3L), 2L)
  res <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_iter(pool, ctx, part, NULL, NULL,
                           list(NULL, NULL, NULL), c(0, 0, 0),
                           list(NULL, NULL, NULL)),
    .emc_wpool_compute = fake, .package = "EMC2")
  expect_equal(calls$n, 2L)          # one fallback per non-empty partition
  expect_equal(dim(res$props), c(3L, 3L))
  expect_equal(res$props[1, ], c(1, 2, 3))
  expect_false(res$alive)
})

# --- growing the pool onto cores freed by finished chains --------------------

test_that("the arena is only built when it can help", {
  expect_null(EMC2:::.emc_core_ctl(1, 8))                     # nothing to donate
  expect_null(EMC2:::.emc_core_ctl(3, 0))
})

test_that("a chain claims the cores released by finished siblings", {
  skip_on_os("windows")
  ctl <- EMC2:::.emc_core_ctl(3, 8)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  expect_equal(ctl$total, 24)

  # Nobody has finished: stay on the static share.
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 8)

  file.create(file.path(ctl$dir, "done_1"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 12)   # 24 / 2 remaining

  file.create(file.path(ctl$dir, "done_2"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 24)   # last chain takes the lot

  # Never drops below the chain's own share, whatever the bookkeeping says.
  file.create(file.path(ctl$dir, "done_3"))
  expect_gte(EMC2:::.emc_cores_now(ctl, 8), 8)
})

test_that("a chain releases its cores even when it fails", {
  skip_on_os("windows")
  ctl <- EMC2:::.emc_core_ctl(3, 4)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  boom <- function() {
    on.exit(EMC2:::.emc_core_release(ctl), add = TRUE)
    stop("chain failed")
  }
  expect_error(boom(), "chain failed")
  expect_length(list.files(ctl$dir), 1)
})

test_that("absent control leaves the static budget untouched", {
  expect_equal(EMC2:::.emc_cores_now(NULL, 8), 8)
  expect_silent(EMC2:::.emc_core_release(NULL))
})

test_that("growing the pool keeps the workers that were already running", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not fork a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)

  grown <- EMC2:::.emc_wpool_grow(pool, 4, list(tag = "ctx"))
  expect_equal(grown$n, 4)
  expect_identical(grown$wcs[[1]], pool$wcs[[1]])
  pool <- grown
  # Shrinking is not a thing: workers are only ever added within a block.
  expect_equal(EMC2:::.emc_wpool_grow(pool, 1, list(tag = "ctx"))$n, 4)

  for (w in 1:4) { serialize(list(subs = w), pool$wcs[[w]]); flush(pool$wcs[[w]]) }
  got <- lapply(1:4, function(w) unserialize(pool$rcs[[w]]))
  expect_true(all(vapply(got, function(g) !is.null(g$failed), logical(1))))
})
