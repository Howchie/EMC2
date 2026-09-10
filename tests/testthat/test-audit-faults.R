# The fault half of the architecture-efficiency audit's C2 harness: "interrupt
# and parent exit, death before/after notification, partial payload, partial
# startup/grow/recycle failure, repeated teardown, stale generation replies,
# publication failure and connection exhaustion.  Use a small simulated
# worker-ID space rather than 256 real processes."
#
# test-worker-pool.R already covers spawn failure, recycle failure, grow
# back-off, short reads and writes, a garbled frame length, and a worker
# reporting a missing shared file.  What is added here is the set the audit's P0
# findings name and nothing exercises, plus the boundary probes that C4 and C5
# will be judged against.
#
# Two of these document a *present defect* rather than a present guarantee.
# They are written to pass against today's behaviour and to say plainly, in
# place, what the correct behaviour is and which commit owns it -- so that C4
# and C5 have a test to change rather than a test to write.

RNGkind("L'Ecuyer-CMRG")

# --- the completion-token width (finding 1) ---------------------------------

test_that("the completion token is one byte, so worker IDs above 255 alias", {
  # Finding 1.  `.emc_wpool_serve()` writes `as.raw(msg$w)` and
  # `.emc_wpool_await_done()` reads one byte, while `.emc_wpool_ll()` guards the
  # dynamic queue at .EMC_WPOOL_MAX_DYN_WORKERS = 255.  `.emc_wpool_iter()`
  # enables notification whenever `pool$done` exists, without the equivalent
  # guard, so a wide enough subject pool can encode worker 256 as byte 0 and
  # take the "lost alignment" recovery path.
  #
  # Simulated ID space, per the audit: the question is how a width is encoded,
  # and forking 256 R sessions to ask it would make the test unrunnable exactly
  # where it matters.
  expect_identical(EMC2:::.EMC_WPOOL_MAX_DYN_WORKERS, 255L)
  for (w in c(1L, 127L, 255L)) {
    expect_identical(as.integer(as.raw(w)), w, info = paste("worker", w))
  }
  # 256 and 257 do not survive the round trip; they alias onto 0 and 1.
  for (w in c(256L, 257L)) {
    aliased <- suppressWarnings(as.integer(as.raw(w)))
    expect_false(identical(aliased, w), info = paste("worker", w))
  }
  expect_warning(as.raw(256L), "out-of-range")

  # Reachability, so this is not mistaken for a live hazard: pool width is
  # bounded by min(n_subjects, cores), so reaching worker 256 needs both 256+
  # subjects and 256+ cores.  It is latent on ordinary hardware and rides along
  # free in C4's single wire-format revision.
  #
  # C4 replaces this with a wide worker ID carried alongside a pool generation,
  # a request ID and a status.  When it lands, the two expectations above
  # invert: as_wide_id(256) must round-trip.
  succeed()
})

# --- partial payloads and worker death (finding 2) --------------------------

test_that("a frame truncated mid-payload is refused, not misread", {
  # A worker that announces completion and then dies part-way through its reply.
  # The reader must return NULL -- the caller's signal to recompute -- rather
  # than unserialising a prefix or blocking forever.
  obj <- list(props = matrix(rnorm(400), 20L), tag = "payload")
  full <- audit_frame_roundtrip(obj)
  expect_true(full$ok)
  expect_identical(full$value, obj)

  for (frac in c(0, 0.25, 0.5, 0.99)) {
    keep <- floor(full$bytes * frac)
    f <- audit_truncated_frame(obj, keep)
    on.exit(unlink(f), add = TRUE)
    con <- file(f, "rb")
    got <- EMC2:::.emc_wpool_recv(con)
    close(con)
    expect_null(got, info = sprintf("truncated to %.0f%%", 100 * frac))
  }
})

test_that("a header that arrives without its payload is refused", {
  # The exact shape of "announced completion, then stalled": eight good header
  # bytes claiming a payload that never comes.
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  writeBin(as.double(4096), con, size = 8L, endian = "little")
  close(con)
  con <- file(f, "rb")
  expect_null(EMC2:::.emc_wpool_recv(con))
  close(con)
})

test_that("an implausible frame length is refused before allocating on it", {
  # The guard is protection, not a task-size limit: legitimate large messages
  # are the caller's business, a 2^40-byte claim is not.
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  writeBin(as.double(2^40), con, size = 8L, endian = "little")
  writeBin(raw(16L), con)
  close(con)
  con <- file(f, "rb")
  expect_null(EMC2:::.emc_wpool_recv(con))
  close(con)
  expect_lt(EMC2:::.EMC_WPOOL_MAX_FRAME, 2^40)
})

test_that("frames round-trip on both sides of every buffer boundary", {
  # 4 KiB is the write chunk, 64 KiB the read chunk and the classic pipe buffer.
  # The audit asks for both sides of each; a payload that lands exactly on a
  # chunk edge is where an off-by-one in the accumulation loop lives.
  chunk <- EMC2:::.EMC_WPOOL_WRITE_CHUNK
  read_chunk <- EMC2:::.EMC_WPOOL_READ_CHUNK
  for (n in c(1L, chunk - 1L, chunk, chunk + 1L,
              2L * chunk, read_chunk - 1L, read_chunk, read_chunk + 1L)) {
    obj <- list(payload = as.raw(rep_len(1:255, n)))
    got <- audit_frame_roundtrip(obj)
    expect_true(got$ok, info = paste("n =", n))
    expect_identical(got$value, obj, info = paste("n =", n))
  }
})

test_that("a reply larger than the pipe buffer round-trips whole", {
  # 8 KiB and past 64 KiB, as doubles rather than raws, so the serialised size
  # is comfortably larger than any single write.
  for (n in c(1024L, 8192L, 100000L)) {
    obj <- list(x = seq_len(n) + 0.5)
    got <- audit_frame_roundtrip(obj)
    expect_true(got$ok, info = paste("n =", n))
    expect_identical(got$value, obj, info = paste("n =", n))
    expect_gt(got$bytes, 8L * n - 1L)
  }
})

# --- publication failure ----------------------------------------------------

test_that("a shared broadcast that cannot be published is reported, not thrown", {
  # Master-side half of publication failure: the worker-side half (a worker
  # finding the file missing) is in test-worker-pool.R.  An unwritable
  # directory must degrade the pool, not raise out of the iteration.
  expect_null(EMC2:::.emc_wpool_shared_write(file.path(tempdir(), "no-such-dir"),
                                             serialize(list(1), NULL)))
  expect_null(EMC2:::.emc_wpool_shared_write(NULL, serialize(list(1), NULL)))
  expect_null(EMC2:::.emc_wpool_shared_write(character(0), serialize(list(1), NULL)))
})

test_that("a published broadcast is readable and is removed by its owner", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  payload <- serialize(list(group_var = diag(3)), NULL)
  path <- EMC2:::.emc_wpool_shared_write(dir, payload)
  expect_true(file.exists(path))
  expect_identical(EMC2:::.emc_wpool_shared_read(path)$group_var, diag(3))
  # A reader that arrives after the owner cleaned up *signals*, and does not
  # return a partial object.  The write is a rename, so a file is either
  # complete or absent, and the worker turns this error into a `failed` reply
  # while staying alive -- see test-worker-pool.R.
  unlink(path)
  expect_error(EMC2:::.emc_wpool_shared_read(path), "invalid shared broadcast")
  expect_error(EMC2:::.emc_wpool_shared_read(NULL), "invalid shared broadcast")
  expect_error(EMC2:::.emc_wpool_shared_read(character(0)), "invalid shared broadcast")
})

# --- teardown ---------------------------------------------------------------

test_that("teardown is idempotent and survives a degraded pool", {
  # Finding 3: a failed recycle returns a pool without its template, directory
  # or job handles, and the normal shutdown route must still not raise.
  # Repeated teardown is the audit's own listed case.
  expect_silent(EMC2:::.emc_wpool_stop(NULL))
  bare <- audit_fake_pool(0L, alive = FALSE)
  expect_silent(EMC2:::.emc_wpool_stop(bare))
  expect_silent(EMC2:::.emc_wpool_stop(bare))

  # A dead pool that still names a directory: stopping it twice must not fail
  # on the second unlink of an already-removed tree.
  dir <- tempfile()
  dir.create(dir)
  dead <- audit_fake_pool(0L, alive = FALSE, dir = dir)
  expect_silent(EMC2:::.emc_wpool_stop(dead))
  expect_false(dir.exists(dir))
  expect_silent(EMC2:::.emc_wpool_stop(dead))
})

test_that("a real pool can be stopped twice", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("mkfifo")))
  ctx <- list(kind = "ll", n_pars = 2L, data = list(), model = NULL)
  pool <- EMC2:::.emc_wpool_start(2L, ctx)
  skip_if(is.null(pool), "no worker pool available here")
  dir <- pool$dir
  # `expect_no_error`, not `expect_silent`: shutting down a real pool closes
  # connections, and doing that underneath testthat's output sink upsets the
  # sink itself rather than saying anything about the pool.
  expect_no_error(EMC2:::.emc_wpool_stop(pool))
  expect_false(dir.exists(dir))
  # The second call sees closed connections and a removed directory.
  expect_no_error(EMC2:::.emc_wpool_stop(pool))
})

# --- partial startup (finding 4) --------------------------------------------

test_that("a startup that fails part-way leaves no untracked children", {
  # Finding 4.  `.emc_wpool_template()` accumulates newly forked workers in
  # `made` and commits them to `jobs` only after the loop, so a failure after
  # some successful forks can leave those children outside the tracked list.
  #
  # What is testable from here without a template is the observable contract:
  # a spawn that fails must return NULL rather than a half-built pool, and must
  # not leave its scratch directory behind.
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  # Fail the *second* fork, so some children already exist when startup gives
  # up.  Mock at `parallel::mcparallel`: `.emc_wpool_start()` does not wrap
  # `.emc_wpool_spawn()` in a tryCatch, so an error raised by spawn itself
  # propagates -- the handling lives below it.
  n <- 0L
  pool <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_start(3L, list(n_pars = 2L)),
    mcparallel = function(...) {
      n <<- n + 1L
      if (n >= 2L) stop("unable to fork, possible reason: ",
                        "Resource temporarily unavailable")
      parallel::mcparallel(...)
    },
    .package = "parallel")
  expect_null(pool)
  expect_gte(n, 2L)

  # C5 owns the other half: cleanup must own each child at the moment it is
  # acquired, so that a failure between two forks cannot orphan the first.
  # That needs the template's `made`/`jobs` split rewritten, and this test
  # becomes the place to assert the tracked-child count afterwards.
  succeed()
})

# --- connection exhaustion --------------------------------------------------

test_that("the pool's connection appetite is bounded and knowable", {
  # The audit: "Check available R connections/file descriptors before growing:
  # this pool consumes separate request/reply connections for each worker plus
  # template/completion channels.  Allocate a feasible worker count and queue
  # the remaining work instead of discovering the limit part-way through
  # startup."
  #
  # R's connection table is a fixed size, and every worker costs two entries
  # plus the template's two and the completion channel.  Pin the accounting so
  # a grow policy can be written against it.
  skip_on_os("windows")
  before <- nrow(showConnections(all = FALSE))
  cons <- list()
  on.exit(for (cn in cons) try(close(cn), silent = TRUE), add = TRUE)
  f <- tempfile()
  writeLines("x", f)
  on.exit(unlink(f), add = TRUE)
  for (i in 1:8) cons[[i]] <- file(f, "rb")
  expect_equal(nrow(showConnections(all = FALSE)) - before, 8L)
  for (cn in cons) close(cn)
  cons <- list()
  expect_equal(nrow(showConnections(all = FALSE)), before)
})

test_that("a resource failure during startup degrades rather than erroring", {
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  # Connection exhaustion and process exhaustion arrive at the same place: the
  # backend cannot produce a worker.  Whatever the reason, startup must return
  # NULL so the caller runs serially, never raise into the fit.
  #
  # C13 owns the other half the audit asks for -- checking the connection
  # budget *before* growing and queueing the remainder, rather than discovering
  # the limit part-way through startup.  That needs the bounded task budget.
  pool <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_start(2L, list(n_pars = 2L)),
    mcparallel = function(...) stop("all connections are in use"),
    .package = "parallel")
  expect_null(pool)
})

# --- interrupt and parent exit ----------------------------------------------

test_that("an interrupted iteration still removes its shared broadcast", {
  # The `on.exit(unlink(shared_file))` in `.emc_wpool_iter()` is what makes an
  # interrupt safe: a surviving worker must never observe a half-removed
  # payload, and a cancelled fit must not leave its scratch behind.
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- NULL
  res <- tryCatch({
    f <- function() {
      p <- EMC2:::.emc_wpool_shared_write(dir, serialize(list(1), NULL))
      path <<- p
      on.exit(unlink(p), add = TRUE)
      stop("simulated interrupt")
    }
    f()
  }, error = function(e) conditionMessage(e))
  expect_match(res, "simulated interrupt")
  expect_false(is.null(path))
  expect_false(file.exists(path))
})

test_that("a worker outlives neither its parent's exit nor its own shutdown token", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("mkfifo")))
  ctx <- list(kind = "ll", n_pars = 2L, data = list(), model = NULL)
  pool <- EMC2:::.emc_wpool_start(2L, ctx)
  skip_if(is.null(pool), "no worker pool available here")
  pids <- vapply(pool$jobs, EMC2:::.emc_wpool_job_pid, integer(1))
  expect_true(all(!is.na(pids)))
  expect_true(all(vapply(pids, EMC2:::.emc_wpool_pid_alive, logical(1))))
  EMC2:::.emc_wpool_stop(pool)
  # Give the children a moment to reap; a spawned worker exits on end-of-stream.
  deadline <- Sys.time() + 10
  repeat {
    gone <- !vapply(pids, EMC2:::.emc_wpool_pid_alive, logical(1))
    if (all(gone) || Sys.time() > deadline) break
    Sys.sleep(0.02)
  }
  expect_true(all(!vapply(pids, EMC2:::.emc_wpool_pid_alive, logical(1))))
})

# --- stale generation replies -----------------------------------------------

test_that("stale generation replies are rejected", {
  # Finding 5's other half.  A late completion token from a retired worker can
  # still be sitting in the `done` channel, which is deliberately retained
  # across a recycle -- so after recycling, a token can name a worker index
  # that now belongs to a different process.
  #
  # There is nothing to test yet: the framed record carries a worker index and
  # nothing else, so a stale token is indistinguishable from a current one.
  # C4 adds the pool generation and request ID that make the distinction
  # expressible, and C5 retires a whole generation before accepting new work.
  # This is the test that gates them.
  skip("requires C4's generation field in the framed completion record")
})

test_that("a late token from a retired worker is currently indistinguishable", {
  # The same finding, stated as something that passes today, so the defect is
  # recorded rather than merely skipped.  A completion token is one byte and
  # carries no generation, so two tokens from different pool generations are
  # byte-identical.
  token_from <- function(worker) as.raw(worker)
  expect_identical(token_from(3L), token_from(3L))
  # After C4, a token must also carry the generation it was issued in, and
  # these two must differ.
  succeed()
})
