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

test_that("worker model closures retain only the realised specification", {
  payload <- raw(1024^2)
  model <- local({
    unused <- payload
    spec <- list(c_name = "demo", p_types = c(a = 1))
    function() spec
  })
  slim <- EMC2:::.emc_wpool_slim_model(model)
  expect_identical(slim(), model())
  expect_lt(length(serialize(slim, NULL)), length(serialize(model, NULL)) / 10)
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
      EMC2:::.emc_wpool_send(pool$wcs[[w]], list(subs = w, echo = k))
    }
    got <- lapply(1:2, function(w) EMC2:::.emc_wpool_recv(pool$rcs[[w]]))
    # No real ctx, so compute() errors and the worker reports it rather than dying.
    expect_true(all(vapply(got, function(g) !is.null(g$failed), logical(1))))
  }
  EMC2:::.emc_wpool_stop(pool)
  expect_false(dir.exists(pool$dir))
  on.exit(NULL)
})

test_that("spawned workers descend from a history-free template", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  skip_if(!file.exists("/proc/self/stat"), "process ancestry needs procfs")
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not start a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)
  skip_if(!identical(pool$backend, "spawn"),
          "installed clean-process backend is not active")

  parent_pid <- function(pid) {
    stat <- readLines(sprintf("/proc/%d/stat", pid), n = 1L)
    # The executable name is parenthesised and may contain spaces; fields after
    # the final ')' begin with state and parent PID.
    tail <- sub("^.*\\) ", "", stat)
    as.integer(strsplit(tail, " ", fixed = TRUE)[[1L]][2L])
  }
  worker_pids <- vapply(pool$jobs, function(x) x$pid, integer(1))
  expect_true(all(vapply(worker_pids, parent_pid, integer(1)) ==
                    pool$template$job$pid))
  expect_false(any(vapply(worker_pids, parent_pid, integer(1)) == Sys.getpid()))
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

  # Size alone is not what used to break this.  unserialize() on a connection
  # asks for a whole XDR chunk at a time, so a double vector past 8096 elements
  # (64,768 bytes) short-read and killed the worker mid-message however the rest
  # of the message was shaped -- an alpha block for a wide group_design gets
  # there long before anything looks large.  Cover both shapes and both sides.
  cases <- list(
    small = list(subs = 1L, pad = 1:10),
    doubles = list(subs = 1L, group_chol = matrix(0, 400, 400)),
    # sample.int(), not seq_len(): a compact ALTREP sequence serialises to a few
    # bytes and would never reach the integer path the bug lived on.
    integers = list(subs = 1L, pad = sample.int(1000L, 200000L, replace = TRUE)),
    just_over = list(subs = 1L, pad = numeric(8097)),
    mixed = list(subs = 1L, a = matrix(rnorm(40000), 200),
                 b = sample.int(1000L, 50000L, replace = TRUE),
                 c = serialize(matrix(0, 300, 300), NULL))
  )
  expect_gt(length(serialize(cases$doubles, NULL)), 65536L)
  for (nm in names(cases)) {
    expect_true(EMC2:::.emc_wpool_send(pool$wcs[[1]], cases[[nm]]), label = nm)
    # No real ctx, so the worker reports an error -- but it received the whole
    # message and is still listening, which is the point.
    expect_false(is.null(EMC2:::.emc_wpool_recv(pool$rcs[[1]])$failed),
                 label = nm)
  }
})

test_that("the frame helpers survive short reads and short writes", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  # Straight at the transport, with no pool around it: a FIFO hands back at most
  # one 8 kB buffer per read, so anything bigger arrives in pieces.
  dir <- tempfile("emc_frame_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "p")
  skip_if(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE) != 0)

  obj <- list(x = matrix(rnorm(160000), 400),
              i = sample.int(1000L, 120000L, replace = TRUE),
              s = "tail marker")
  child <- parallel::mcparallel({
    wc <- fifo(path, "wb", blocking = TRUE)
    # Two frames back to back, to show the reader stops at the frame boundary
    # rather than swallowing whatever else is in the pipe.
    EMC2:::.emc_wpool_send(wc, obj)
    EMC2:::.emc_wpool_send(wc, list(second = TRUE))
    close(wc)
    invisible(NULL)
  }, detached = FALSE)
  rc <- fifo(path, "rb", blocking = TRUE)
  first <- EMC2:::.emc_wpool_recv(rc)
  second <- EMC2:::.emc_wpool_recv(rc)
  eof <- EMC2:::.emc_wpool_recv(rc)
  close(rc)
  suppressWarnings(parallel::mccollect(child, wait = TRUE))

  expect_identical(first, obj)
  expect_identical(second, list(second = TRUE))
  expect_null(eof)   # the writer closed; recv reports it rather than hanging
})

test_that("a message that arrives in pieces is still read whole", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  # The real failure was never about size on its own.  serialize() on a
  # connection emits a run of small writes, unserialize() a matching run of
  # reads, and any read that outran the writer was a hard error -- so a message
  # only had to be split across two scheduling slices to kill the worker, which
  # is why 3 kB messages died as readily as 3 MB ones once the box was busy.
  dir <- tempfile("emc_frame_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "p")
  skip_if(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE) != 0)

  obj <- list(a = matrix(rnorm(300), 30),
              b = sample.int(1000L, 200L, replace = TRUE), pm = as.list(1:20))
  payload <- serialize(obj, NULL)
  frame <- c(writeBin(as.double(length(payload)), raw(), size = 8,
                      endian = "little"), payload)
  child <- parallel::mcparallel({
    wc <- fifo(path, "wb", blocking = TRUE)
    cuts <- unique(c(seq(1, length(frame), length.out = 8), length(frame) + 1))
    for (i in seq_len(length(cuts) - 1L)) {
      writeBin(frame[floor(cuts[i]):(floor(cuts[i + 1L]) - 1L)], wc)
      flush(wc)
      Sys.sleep(0.02)
    }
    close(wc)
    invisible(NULL)
  }, detached = FALSE)
  rc <- fifo(path, "rb", blocking = TRUE)
  got <- EMC2:::.emc_wpool_recv(rc)
  close(rc)
  suppressWarnings(parallel::mccollect(child, wait = TRUE))
  expect_identical(got, obj)
})

test_that("a garbled frame length is refused rather than allocated on", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  dir <- tempfile("emc_frame_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "p")
  skip_if(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE) != 0)

  child <- parallel::mcparallel({
    wc <- fifo(path, "wb", blocking = TRUE)
    writeBin(as.double(2^40), wc, size = 8, endian = "little")
    flush(wc); close(wc); invisible(NULL)
  }, detached = FALSE)
  rc <- fifo(path, "rb", blocking = TRUE)
  got <- EMC2:::.emc_wpool_recv(rc)
  close(rc)
  suppressWarnings(parallel::mccollect(child, wait = TRUE))
  expect_null(got)
})

test_that("the shared broadcast has one complete file round trip", {
  skip_on_os("windows")
  dir <- tempfile("emc_shared_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  shared <- list(group_var = diag(3), idx_list = list(c(TRUE, FALSE, TRUE)))
  path <- EMC2:::.emc_wpool_shared_write(dir, serialize(shared, NULL))
  expect_true(file.exists(path))
  expect_identical(EMC2:::.emc_wpool_shared_read(path), shared)
  unlink(path)
  expect_false(file.exists(path))
})

test_that("an iteration shares one path and cleans it after replies", {
  skip_on_os("windows")
  dir <- tempfile("emc_shared_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  set.seed(19)
  sent_paths <- character()
  recv <- function(...) list(props = matrix(0, 3L, 1L), pm = list(NULL),
                             seeds = list(get(".Random.seed", envir = globalenv())),
                             times = 1)
  send <- function(con, msg) {
    sent_paths <<- c(sent_paths, msg$shared_file)
    expect_true(file.exists(msg$shared_file))
    TRUE
  }
  pool <- list(n = 2L, dir = dir, jobs = list(), wcs = list(1, 2),
               rcs = list(1, 2), alive = TRUE, done = NULL)
  pars <- list(alpha = matrix(0, 2L, 2L), subj_mu = matrix(0, 2L, 2L),
               tvar = diag(2))
  cache <- EMC2:::build_group_chol_cache(diag(2), list(c(TRUE, TRUE)))
  ctx <- list(n_pars = 2L, type = "standard")
  res <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_iter(
      pool, ctx, list(1L, 2L), pars, cache, list(NULL, NULL), c(0, 0),
      list(get(".Random.seed", envir = globalenv()),
           get(".Random.seed", envir = globalenv()))),
    .emc_wpool_send = send, .emc_wpool_recv = recv, .package = "EMC2"
  )
  expect_true(res$alive)
  expect_length(unique(sent_paths), 1L)
  expect_false(file.exists(sent_paths[[1L]]))
})

test_that("a worker reports a missing shared file and keeps listening", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not fork a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)
  missing <- file.path(pool$dir, "not-published.bin")
  bad <- list(kind = "particle", shared_file = missing,
              notify = FALSE, w = 1L)
  for (i in 1:2) {
    expect_true(EMC2:::.emc_wpool_send(pool$wcs[[1]], bad))
    got <- EMC2:::.emc_wpool_recv(pool$rcs[[1]])
    expect_match(got$failed, "shared broadcast")
  }
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
  fake <- function(msg, ctx, shared = NULL) {
    calls$n <- calls$n + 1L
    list(props = matrix(msg$subs, ctx$n_pars + 1L, length(msg$subs)),
         pm = as.list(msg$subs), seeds = as.list(msg$subs),
         times = rep(1, length(msg$subs)))
  }
  pool <- list(n = 2, alive = FALSE, wcs = list(), rcs = list())
  ctx <- list(n_pars = 2, type = "standard")
  part <- list(c(1L, 3L), 2L)
  pars <- list(alpha = matrix(0, 2, 3), subj_mu = matrix(0, 2, 3),
               tvar = diag(2))
  res <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_iter(pool, ctx, part, pars, NULL,
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
  expect_null(EMC2:::.emc_core_ctl(3, 8, cores_for_chains = 1))
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

test_that("reallocation respects the outer chain concurrency budget", {
  skip_on_os("windows")
  # Four chains, only two active at once: queued replacements do not make the
  # denominator smaller.  The global budget is 2 * 8, not 4 * 8.
  ctl <- EMC2:::.emc_core_ctl(4, 8, cores_for_chains = 2)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  expect_equal(ctl$total, 16)
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 8)

  file.create(file.path(ctl$dir, "done_1"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 8)
  file.create(file.path(ctl$dir, "done_2"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 8)
  file.create(file.path(ctl$dir, "done_3"))
  expect_equal(EMC2:::.emc_cores_now(ctl, 8), 16)
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

# --- recycling the workers to bound copy-on-write divergence -----------------

test_that("recycling fires on its period and leaves everything else alone", {
  skip_on_os("windows")
  ctx <- list(tag = "ctx")
  live <- list(n = 4L, alive = TRUE)
  # Iteration 1 has nothing to recycle, and only every `every`-th after it does.
  fired <- vapply(1:21, function(i) {
    called <- FALSE
    testthat::with_mocked_bindings(
      { EMC2:::.emc_wpool_recycle(live, i, 10L, ctx); called },
      .emc_wpool_stop = function(pool) { called <<- TRUE; invisible(NULL) },
      .emc_wpool_start = function(n, ctx) list(n = n, alive = TRUE),
      .package = "EMC2")
  }, logical(1))
  expect_equal(which(fired), c(11L, 21L))

  # Off switch, and the guards that must never re-fork.
  expect_identical(EMC2:::.emc_wpool_recycle(live, 11L, 0L, ctx), live)
  expect_identical(EMC2:::.emc_wpool_recycle(live, 11L, NA_integer_, ctx), live)
  expect_identical(EMC2:::.emc_wpool_recycle(list(n = 1L, alive = TRUE), 11L, 10L, ctx),
                   list(n = 1L, alive = TRUE))
  # A pool that has already fallen back *is* rebuilt at a recycle point: it is
  # safe for the same reason growing is (every subject owns its stream), and
  # leaving it dead costs the rest of the block its parallelism.
  dead <- list(n = 4L, alive = FALSE)
  revived <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_recycle(dead, 11L, 10L, ctx),
    .emc_wpool_stop = function(pool) invisible(NULL),
    .emc_wpool_start = function(n, ctx) list(n = n, alive = TRUE),
    .package = "EMC2")
  expect_true(revived$alive)
  expect_equal(revived$n, 4L)

  # But a rebuild that cannot fork backs off rather than retrying every period.
  still_dead <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_recycle(dead, 11L, 10L, ctx),
    .emc_wpool_stop = function(pool) invisible(NULL),
    .emc_wpool_start = function(n, ctx) NULL,
    .package = "EMC2")
  expect_false(still_dead$alive)
  expect_gt(still_dead$rebuild_skip, 0L)

  expect_null(EMC2:::.emc_wpool_recycle(NULL, 11L, 10L, ctx))
})

test_that("clean workers use a less aggressive recycle default", {
  expect_identical(EMC2:::.emc_wpool_recycle_default(list(backend = "spawn")), 50L)
  expect_identical(EMC2:::.emc_wpool_recycle_default(list(backend = "fork")), 10L)
  expect_identical(EMC2:::.emc_wpool_recycle_default(NULL), 10L)
})

test_that("recycling preserves the worker count, including donated cores", {
  skip_on_os("windows")
  asked <- NULL
  got <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_recycle(list(n = 12L, alive = TRUE), 11L, 10L, list()),
    .emc_wpool_stop = function(pool) invisible(NULL),
    .emc_wpool_start = function(n, ctx) { asked <<- n; list(n = n, alive = TRUE) },
    .package = "EMC2")
  expect_equal(asked, 12L)     # not the original share -- what it had grown to
  expect_equal(got$n, 12L)
})

test_that("a failed re-fork degrades instead of losing the block", {
  skip_on_os("windows")
  # The old workers are already stopped by the time the re-fork is attempted,
  # so there is nothing to fall back onto except the master.
  got <- suppressWarnings(testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_recycle(list(n = 4L, alive = TRUE), 11L, 10L, list()),
    .emc_wpool_stop = function(pool) invisible(NULL),
    .emc_wpool_start = function(n, ctx) NULL,
    .package = "EMC2"))
  expect_false(got$alive)
  expect_equal(got$n, 4L)
  # And .emc_wpool_iter() then computes every subject in the master.
  expect_null(got$dir)
})

test_that("real workers survive being recycled and still compute", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  pool <- EMC2:::.emc_wpool_start(2, list(tag = "ctx"))
  skip_if(is.null(pool), "could not fork a pool here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)
  old_dir <- pool$dir
  old_backend <- pool$backend
  old_pids <- vapply(pool$jobs, function(x) x$pid, integer(1))
  old_template_pid <- if (identical(old_backend, "spawn")) {
    pool$template$job$pid
  } else {
    NA_integer_
  }

  pool <- EMC2:::.emc_wpool_recycle(pool, 11L, 10L, list(tag = "ctx"))
  expect_true(isTRUE(pool$alive))
  expect_equal(pool$n, 2)
  if (identical(old_backend, "spawn")) {
    # Clean workers are replaced by the same pristine template: context is not
    # reloaded, and chain history is still absent from their ancestry.
    expect_identical(pool$dir, old_dir)
    expect_equal(pool$template$job$pid, old_template_pid)
    expect_false(any(vapply(pool$jobs, function(x) x$pid, integer(1)) %in%
                     old_pids))
  } else {
    expect_false(dir.exists(old_dir))
  }
  expect_true(dir.exists(pool$dir))

  for (w in 1:2) EMC2:::.emc_wpool_send(pool$wcs[[w]], list(subs = w))
  got <- lapply(1:2, function(w) EMC2:::.emc_wpool_recv(pool$rcs[[w]]))
  expect_true(all(vapply(got, function(g) !is.null(g$failed), logical(1))))
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

  for (w in 1:4) EMC2:::.emc_wpool_send(pool$wcs[[w]], list(subs = w))
  got <- lapply(1:4, function(w) EMC2:::.emc_wpool_recv(pool$rcs[[w]]))
  expect_true(all(vapply(got, function(g) !is.null(g$failed), logical(1))))
})

# --- fork failures degrade instead of killing the chain ---------------------

test_that("a fork failure during spawn returns NULL rather than erroring", {
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  # fork() fails with EAGAIN at the process limit; every other failure in
  # chain_pool.R returns NULL so the caller can fall back, and this must too.
  expect_silent(
    res <- testthat::with_mocked_bindings(
      EMC2:::.emc_wpool_start(3L, list(n_pars = 2)),
      mcparallel = function(...) stop("unable to fork, possible reason: ",
                                      "Resource temporarily unavailable"),
      .package = "parallel"
    )
  )
  expect_null(res)
})

test_that("a fork failure while recycling hands back a dead pool, not an error", {
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  ctx <- list(n_pars = 2)
  pool <- EMC2:::.emc_wpool_start(2L, ctx)
  skip_if(is.null(pool), "no pool available here")
  # Recycle stops the old workers before forking new ones, so an error escaping
  # here would lose the block's work with nothing left to run on.
  fresh <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_recycle(pool, 11L, 10L, ctx),
    mcparallel = function(...) stop("unable to fork"),
    .package = "parallel"
  )
  expect_false(fresh$alive)
  expect_equal(fresh$n, 2L)
})

test_that("a failed grow backs off instead of retrying every iteration", {
  skip_on_os("windows")
  pool <- list(n = 2L, dir = tempfile(), jobs = list(), wcs = list(),
               rcs = list(), alive = TRUE)
  dir.create(pool$dir)
  on.exit(unlink(pool$dir, recursive = TRUE), add = TRUE)
  calls <- 0L
  grow <- function(p) testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_grow(p, 4L, list(n_pars = 2)),
    mcparallel = function(...) { calls <<- calls + 1L; stop("unable to fork") },
    .package = "parallel"
  )
  p <- grow(pool)
  expect_equal(p$n, 2L)          # unchanged, and no error
  expect_equal(calls, 1L)
  expect_gt(p$grow_skip, 0L)
  # The next attempts are skipped outright rather than re-forking.
  before <- calls
  p <- grow(p); p <- grow(p)
  expect_equal(calls, before)
})

test_that("a dead pool is rebuilt at a recycle point when forking works again", {
  skip_on_os("windows")
  ctx <- list(n_pars = 2)
  dead <- list(n = 2L, dir = NULL, jobs = list(), wcs = list(), rcs = list(),
               alive = FALSE)
  fresh <- EMC2:::.emc_wpool_recycle(dead, 11L, 10L, ctx)
  skip_if(is.null(fresh$wcs) || !length(fresh$wcs), "no pool available here")
  on.exit(EMC2:::.emc_wpool_stop(fresh), add = TRUE)
  # Previously `alive = FALSE` was terminal: the rest of the block ran serially
  # in the master even though the machine could fork perfectly well.
  expect_true(fresh$alive)
  expect_equal(fresh$n, 2L)
})

# --- shared group draw is serialised once for all workers -------------------

test_that("the population covariance is reused from the shared cache", {
  population_var <- diag(2)
  shared <- list(group_var = population_var, idx_list = list(c(TRUE, TRUE)))
  raw <- serialize(shared, NULL)
  expect_identical(unserialize(raw), shared)
  # The raw-byte form remains a valid direct-message compatibility path; normal
  # iterations publish these bytes once and put only a file path on each wire.
  msg <- list(subs = 1:2, shared = raw, alpha = matrix(0, 2, 2),
              population_mu = matrix(0, 2, 2), pm = list(NULL, NULL),
              prev_ll = c(0, 0), seeds = list(NULL, NULL))
  expect_type(msg$shared, "raw")
  expect_identical(unserialize(msg$shared)$group_var, population_var)
})

test_that("only the covariance goes on the wire, never its factorisation", {
  # `.chol_factor()` keeps the root and its inverse, so the cache a worker needs
  # is three P x P matrices while the covariance it is derived from is one.
  # Sending the small one is worth several times its own size in pipe traffic
  # on a wide model, so the wire shape is pinned here rather than left to drift.
  set.seed(4)
  p <- 12L
  A <- matrix(rnorm(p * p), p, p)
  population_var <- crossprod(A) + diag(p)
  idx_list <- list(rep(TRUE, p))
  cache <- EMC2:::build_group_chol_cache(population_var, idx_list)
  wire <- serialize(list(group_var = population_var, idx_list = idx_list), NULL)
  expect_lt(length(wire), length(serialize(cache, NULL)) / 2)
  # ... and it is enough: the worker's reconstruction is the master's object.
  expect_identical(
    EMC2:::build_group_chol_cache(unserialize(wire)$group_var,
                                  unserialize(wire)$idx_list),
    cache)
})

test_that("compute() accepts pre-decoded, encoded, and file shared parts", {
  # The master's fallback path already holds the decoded object; the workers
  # only ever have the bytes.  Both must reach the same arguments.
  seen <- new.env(parent = emptyenv())
  ctx <- list(n_pars = 1L, data = list(a = 1), model = NULL, stage = "sample",
              type = "standard", tune = list(), marginalise = NULL,
              r_cores = 1L, chol_caches = list(a = NULL))
  # The two forms deliberately differ: the master already holds the built cache
  # and hands it over, while the worker is sent the covariance alone and has to
  # rebuild it.  Both must arrive at the same arguments.
  population_var <- matrix(9)
  idx_list <- list(TRUE)
  cache <- EMC2:::build_group_chol_cache(population_var, idx_list)
  shared <- list(group_var = population_var, idx_list = idx_list,
                 group_chol = cache)
  dir <- tempfile("emc_shared_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  shared_file <- EMC2:::.emc_wpool_shared_write(
    dir, serialize(shared[c("group_var", "idx_list")], NULL))
  msg <- list(subs = "a",
              shared = serialize(shared[c("group_var", "idx_list")], NULL),
              alpha = matrix(3, 1, 1), population_mu = matrix(7, 1, 1),
              pm = list(NULL),
              prev_ll = 0, seeds = list(get(".Random.seed", envir = globalenv())))
  fake <- function(s, data, pm_settings, ..., parameters, group_chol,
                   current_alpha, population_mu, population_var) {
    seen$parameters <- parameters; seen$group_chol <- group_chol
    seen$current_alpha <- current_alpha
    seen$population_mu <- population_mu
    seen$population_var <- population_var
    list(proposal = 0, ll = 0, pm_settings = NULL)
  }
  msg_file <- msg
  msg_file$shared <- NULL
  msg_file$shared_file <- shared_file
  msgs <- list(msg, msg_file)
  for (i in seq_along(msgs)) {
    current <- msgs[[i]]
    seen$parameters <- NULL; seen$group_chol <- NULL
    testthat::with_mocked_bindings(
      EMC2:::.emc_wpool_compute(current, ctx,
                                shared = if (i == 1L) NULL else shared),
      safe_new_particle = fake, .package = "EMC2"
    )
    expect_null(seen$parameters)
    expect_identical(seen$group_chol, cache)
    expect_equal(seen$current_alpha, 3)
    expect_equal(seen$population_mu, 7)
    expect_equal(seen$population_var, population_var)
  }
})

# --- single-subject likelihood mode -----------------------------------------

test_that("routing probes both arms and then keeps the cheaper one", {
  EMC2:::.emc_ll_route_reset()
  n <- EMC2:::.EMC_LL_PROBE_N
  # Arms are interleaved: serial, pool, serial, pool, ... Each call is kept
  # under .EMC_LL_OBVIOUS so the comparison, not the shortcut, decides.
  for (i in seq_len(n)) {
    expect_false(EMC2:::.emc_ll_route_use_pool())
    EMC2:::.emc_ll_route_record(FALSE, 0.4, 10L)
    expect_true(EMC2:::.emc_ll_route_use_pool())
    EMC2:::.emc_ll_route_record(TRUE, 0.04, 10L)
  }
  # The pool arm was 10x cheaper per particle, so it is kept.
  expect_true(EMC2:::.emc_ll_route_use_pool())
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "pool")

  # And the reverse: a pool that loses is dropped, ties included.
  EMC2:::.emc_ll_route_reset()
  for (i in seq_len(n)) {
    EMC2:::.emc_ll_route_record(FALSE, 0.1, 10L)
    EMC2:::.emc_ll_route_record(TRUE, 0.1, 10L)
  }
  expect_false(EMC2:::.emc_ll_route_use_pool())
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "serial")
  EMC2:::.emc_ll_route_reset()
})

test_that("an expensive serial call takes the pool without finishing the probe", {
  # A call costing this much cannot lose to a round trip measured in
  # milliseconds, and probing it would spend five more such calls on one core.
  EMC2:::.emc_ll_route_reset()
  EMC2:::.emc_ll_route_record(FALSE, EMC2:::.EMC_LL_OBVIOUS, 10L)
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "pool")
  expect_true(EMC2:::.emc_ll_route_use_pool())
  # A cheap call still has to earn it.
  EMC2:::.emc_ll_route_reset()
  EMC2:::.emc_ll_route_record(FALSE, EMC2:::.EMC_LL_OBVIOUS / 2, 10L)
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "probe")
  EMC2:::.emc_ll_route_reset()
})

test_that("a settled route is re-probed rather than kept forever", {
  EMC2:::.emc_ll_route_reset()
  n <- EMC2:::.EMC_LL_PROBE_N
  for (i in seq_len(n)) {
    EMC2:::.emc_ll_route_record(FALSE, 0.4, 10L)
    EMC2:::.emc_ll_route_record(TRUE, 0.04, 10L)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "pool")
  for (i in seq_len(EMC2:::.EMC_LL_REPROBE)) {
    EMC2:::.emc_ll_route_record(TRUE, 0.04, 10L)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "probe")
  EMC2:::.emc_ll_route_reset()
})

test_that("splitting a subject's particles over the pool changes no number", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("mkfifo")))
  dat <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))
  emc <- suppressMessages(make_emc(dat, des, type = "single"))
  s1 <- emc[[1]]
  ctx <- list(data = s1$data, model = EMC2:::.emc_wpool_slim_model(s1$model),
              n_pars = s1$n_pars, r_cores = 1L, type = s1$type)
  pool <- EMC2:::.emc_wpool_start(3, ctx)
  skip_if(is.null(pool), "could not start a pool here")
  on.exit({ EMC2:::.emc_ll_pool_clear(); EMC2:::.emc_wpool_stop(pool) }, add = TRUE)

  p <- sampled_pars(des, doMap = FALSE)
  set.seed(11)
  props <- matrix(rep(p, each = 20), nrow = 20, dimnames = list(NULL, names(p)))
  props <- props + matrix(rnorm(length(props), 0, 0.05), nrow = 20)
  ref <- EMC2:::calc_ll_manager(props, s1$data[[1]], s1$model, r_cores = 1)

  EMC2:::.emc_ll_pool_set(pool, 1L)
  got <- EMC2:::.emc_wpool_ll(props, 1L)
  expect_false(is.null(got))
  expect_equal(got, ref)

  # An uneven split must still land every row in its own place.
  odd <- props[1:7, , drop = FALSE]
  expect_equal(EMC2:::.emc_wpool_ll(odd, 1L),
               EMC2:::calc_ll_manager(odd, s1$data[[1]], s1$model, r_cores = 1))

  # Too few rows to be worth a round trip: the caller is told to do it itself.
  expect_null(EMC2:::.emc_wpool_ll(props[1:2, , drop = FALSE], 1L))
  # A subject the workers were not started for never reaches them.
  expect_equal(EMC2:::calc_ll_pooled(props, s1$data[[1]], s1$model, s = 99L), ref)
})

test_that("the likelihood pool only serves the subject it was registered for", {
  EMC2:::.emc_ll_pool_clear()
  # No pool at all: the wrapper is the plain manager, including its own split.
  expect_null(EMC2:::.emc_pool_state$ll_subject)
  expect_null(EMC2:::.emc_wpool_ll(matrix(0, 4, 2), 1L))
})

test_that("the work queue returns exactly what the static split returns", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("mkfifo")))
  dat <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))
  emc <- suppressMessages(make_emc(dat, des, type = "single"))
  s1 <- emc[[1]]
  ctx <- list(data = s1$data, model = EMC2:::.emc_wpool_slim_model(s1$model),
              n_pars = s1$n_pars, r_cores = 1L, type = s1$type)
  pool <- EMC2:::.emc_wpool_start(4, ctx)
  skip_if(is.null(pool), "could not start a pool here")
  on.exit({ EMC2:::.emc_ll_pool_clear(); EMC2:::.emc_wpool_stop(pool) }, add = TRUE)
  skip_if(is.null(pool$done), "no completion channel here")
  EMC2:::.emc_ll_pool_set(pool, 1L)

  p <- sampled_pars(des, doMap = FALSE)
  set.seed(11)
  mk <- function(n) {
    m <- matrix(rep(p, each = n), nrow = n, dimnames = list(NULL, names(p)))
    m + matrix(rnorm(length(m), 0, 0.4), nrow = n)
  }
  # Alternate the two paths.  A token left unread by one round, or a reply left
  # in a pipe, would be collected by the *next* call as though it answered that
  # one -- so the failure this guards against shows up a call later, not here.
  for (i in 1:8) {
    props <- mk(c(24L, 40L, 64L)[(i %% 3L) + 1L])
    ref <- EMC2:::calc_ll_manager(props, s1$data[[1]], s1$model, r_cores = 1)
    expect_equal(EMC2:::.emc_wpool_ll(props, 1L, dynamic = (i %% 2L == 0L)), ref)
  }
  expect_true(isTRUE(EMC2:::.emc_pool_state$ll_pool$alive))

  # Cutting finer must not move a number either.
  props <- mk(50L)
  ref <- EMC2:::calc_ll_manager(props, s1$data[[1]], s1$model, r_cores = 1)
  for (g in c(1, 4, 32)) {
    withr::local_options(emc2.ll_queue_grain = g)
    expect_equal(EMC2:::.emc_wpool_ll(props, 1L, dynamic = TRUE), ref)
  }

  # Too few rows per worker for a queue to have anything to schedule: the
  # static split answers instead, and says so.
  small <- mk(3L * pool$n)
  expect_equal(EMC2:::.emc_wpool_ll(small, 1L, dynamic = TRUE),
               EMC2:::calc_ll_manager(small, s1$data[[1]], s1$model, r_cores = 1))
  expect_false(EMC2:::.emc_pool_state$ll_dynamic)
})

test_that("the split strategy is probed only once the pool arm has won", {
  EMC2:::.emc_ll_route_reset()
  # While serial-versus-pool is still open, every pooled call is static: the
  # outer comparison must not be made against an arm that is itself alternating.
  expect_false(EMC2:::.emc_ll_route_use_dynamic())
  EMC2:::.emc_ll_route_record(FALSE, EMC2:::.EMC_LL_OBVIOUS, 10L)
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "pool")
  expect_false(EMC2:::.emc_ll_route_use_dynamic())   # static first, then alternate

  # A queue that wins is kept; the arms interleave, as above, because one call's
  # cost varies with the particles it draws.
  for (i in 1:EMC2:::.EMC_LL_PROBE_N) {
    EMC2:::.emc_ll_route_record(TRUE, 1.0, 10L, dynamic = FALSE)
    EMC2:::.emc_ll_route_record(TRUE, 0.5, 10L, dynamic = TRUE)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$split, "dynamic")
  expect_true(EMC2:::.emc_ll_route_use_dynamic())

  # A queue that loses is dropped -- the LNR/LBA case, where the round trips
  # cost more than the whole likelihood.
  EMC2:::.emc_ll_route_reset()
  EMC2:::.emc_ll_route_record(FALSE, EMC2:::.EMC_LL_OBVIOUS, 10L)
  for (i in 1:EMC2:::.EMC_LL_PROBE_N) {
    EMC2:::.emc_ll_route_record(TRUE, 0.5, 10L, dynamic = FALSE)
    EMC2:::.emc_ll_route_record(TRUE, 1.0, 10L, dynamic = TRUE)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$split, "static")
  expect_false(EMC2:::.emc_ll_route_use_dynamic())
  EMC2:::.emc_ll_route_reset()
})

test_that("a settled split decision survives an outer re-probe", {
  EMC2:::.emc_ll_route_reset()
  EMC2:::.emc_ll_route_record(FALSE, EMC2:::.EMC_LL_OBVIOUS, 10L)
  for (i in 1:EMC2:::.EMC_LL_PROBE_N) {
    EMC2:::.emc_ll_route_record(TRUE, 1.0, 10L, dynamic = FALSE)
    EMC2:::.emc_ll_route_record(TRUE, 0.5, 10L, dynamic = TRUE)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$split, "dynamic")
  # The split is measured only while the pool arm runs, so wiping it whenever
  # the serial-versus-pool decision cycles would strand the probe.
  for (i in 1:EMC2:::.EMC_LL_REPROBE) {
    EMC2:::.emc_ll_route_record(TRUE, 0.5, 10L, dynamic = TRUE)
  }
  expect_identical(EMC2:::.emc_pool_state$ll_route$mode, "probe")
  expect_identical(EMC2:::.emc_pool_state$ll_route$split, "dynamic")
  EMC2:::.emc_ll_route_reset()
})

test_that("growing the pool keeps the completion channel", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  ctx <- list(tag = "ctx")
  pool <- EMC2:::.emc_wpool_start(2, ctx)
  skip_if(is.null(pool), "no pool available here")
  on.exit(EMC2:::.emc_wpool_stop(pool), add = TRUE)
  skip_if(is.null(pool$done), "no completion channel here")
  grown <- EMC2:::.emc_wpool_grow(pool, 4, ctx)
  # Rebuilding the pool list without this dropped the FIFO on the floor: it
  # leaked, and the likelihood queue quietly fell back to a static split for the
  # rest of the block -- just when donated cores made scheduling matter most.
  expect_identical(grown$done, pool$done)
  pool <- grown
})
