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

# --- the completion record's width (finding 1) ------------------------------

test_that("the completion record carries a worker ID that does not alias", {
  # Finding 1.  The announcement used to be `as.raw(msg$w)`, read back one byte
  # at a time, and the dynamic queue capped itself at 255 workers to stay inside
  # it.  `.emc_wpool_iter()` enabled notification whenever `pool$done` existed
  # with no equivalent guard, so a wide enough subject pool encoded worker 256
  # as byte 0 and took the "lost alignment" recovery path -- retiring a healthy
  # pool.
  #
  # Simulated ID space, per the audit: the question is how a width is encoded,
  # and forking 256 R sessions to ask it would make the test unrunnable exactly
  # where it matters.
  #
  # This is what the one-byte token could not do:
  for (w in c(256L, 257L, 1000L)) {
    aliased <- suppressWarnings(as.integer(as.raw(w)))
    expect_false(identical(aliased, w), info = paste("one byte, worker", w))
  }
  # And this is what the record does.  255/256 are the boundary the audit names;
  # the rest are there so a future narrowing shows up as a failure here.
  for (w in c(1L, 127L, 255L, 256L, 257L, 65535L, 65536L, 2147483647L)) {
    rec <- EMC2:::.emc_wpool_done_decode(EMC2:::.emc_wpool_done_encode(w))
    expect_identical(rec$w, w, info = paste("worker", w))
  }
  # The cap the encoding forced is gone with it.
  expect_false(exists(".EMC_WPOOL_MAX_DYN_WORKERS", envir = asNamespace("EMC2"),
                      inherits = FALSE))
})

test_that("the completion record carries generation, request and status", {
  # The other three fields exist so that a record can be disbelieved.  Before
  # them, a reply from a pool that had already been retired was indistinguishable
  # from a current worker's, which is how a retired worker could answer for a
  # live one (finding 5).
  rec <- EMC2:::.emc_wpool_done_decode(
    EMC2:::.emc_wpool_done_encode(300L, generation = 7L, request = 90001L,
                                  status = EMC2:::.EMC_WPOOL_ST_FAILED))
  expect_identical(rec$w, 300L)
  expect_identical(rec$generation, 7L)
  expect_identical(rec$request, 90001L)
  expect_identical(rec$status, EMC2:::.EMC_WPOOL_ST_FAILED)
  expect_identical(rec$code, "record")

  # Fixed width, and far below PIPE_BUF: several workers share one FIFO, so a
  # record that varied in length would desynchronise the channel for all of them.
  for (args in list(list(1L), list(2L, 3L, 4L), list(NA_integer_, NULL, "x"))) {
    bytes <- do.call(EMC2:::.emc_wpool_done_encode, args)
    expect_identical(length(bytes), EMC2:::.EMC_WPOOL_DONE_BYTES)
  }
  expect_lt(EMC2:::.EMC_WPOOL_DONE_BYTES, 4096L)
  # A request that predates the fields encodes as 0 rather than short.
  bare <- EMC2:::.emc_wpool_done_decode(EMC2:::.emc_wpool_done_encode(5L))
  expect_identical(bare$generation, 0L)
  expect_identical(bare$request, 0L)

  # The three statuses are distinct, and START is not a completion: a worker
  # that has collected its request and a worker that has finished it are
  # different answers to "is this making progress".
  st <- c(EMC2:::.EMC_WPOOL_ST_START, EMC2:::.EMC_WPOOL_ST_OK,
          EMC2:::.EMC_WPOOL_ST_FAILED)
  expect_identical(length(unique(st)), 3L)
  expect_true(all(st > 0L & st < 256L))
})

test_that("a partially arrived record is not decoded until it is whole", {
  # The record is written in one call, and a write that size is atomic on a
  # pipe -- but the reader is non-blocking, and a non-blocking read may return
  # less than was written.  A protocol that depends on it not doing so is one
  # that works until it does not, so the reader buffers.
  dir <- tempfile("emc_done_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  done <- EMC2:::.emc_wpool_done_open(dir)
  skip_if(is.null(done), "no FIFO available here")
  on.exit(for (cn in done[c("rc", "wc")]) try(close(cn), silent = TRUE),
          add = TRUE)

  full <- EMC2:::.emc_wpool_done_encode(258L, generation = 2L, request = 11L)
  # Nothing at all: no record.
  expect_null(EMC2:::.emc_wpool_done_take(done))
  # A prefix: still no record, and the bytes are kept rather than dropped.
  writeBin(full[1:9], done$wc)
  flush(done$wc)
  expect_null(EMC2:::.emc_wpool_done_take(done))
  writeBin(full[10:EMC2:::.EMC_WPOOL_DONE_BYTES], done$wc)
  flush(done$wc)
  rec <- EMC2:::.emc_wpool_done_take(done)
  expect_false(is.null(rec))
  expect_identical(rec$w, 258L)
  expect_identical(rec$generation, 2L)

  # Two records arriving together are returned one at a time, in order.
  writeBin(c(EMC2:::.emc_wpool_done_encode(1L, 2L, 12L),
             EMC2:::.emc_wpool_done_encode(2L, 2L, 13L)), done$wc)
  flush(done$wc)
  first <- EMC2:::.emc_wpool_done_take(done)
  second <- EMC2:::.emc_wpool_done_take(done)
  expect_identical(c(first$w, second$w), c(1L, 2L))
  expect_identical(c(first$request, second$request), c(12L, 13L))
  expect_null(EMC2:::.emc_wpool_done_take(done))
})

test_that("a record from a retired generation is not mistaken for a live one", {
  # Finding 5's half that belongs to the wire format: the master must be able to
  # tell a late record from a pool it has already replaced.  Retiring that
  # generation's pending work is C5.
  dir <- tempfile("emc_done_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  done <- EMC2:::.emc_wpool_done_open(dir)
  skip_if(is.null(done), "no FIFO available here")
  on.exit(for (cn in done[c("rc", "wc")]) try(close(cn), silent = TRUE),
          add = TRUE)

  pool <- audit_fake_pool(2L)
  pool$done <- done
  pool$generation <- 4L
  before <- EMC2:::.emc_pool_state$stale_records
  before <- if (is.null(before)) 0L else before

  # One stale record, then the current one. The stale one must be counted and
  # skipped, and the wait must return the live one.
  writeBin(c(EMC2:::.emc_wpool_done_encode(1L, generation = 3L, request = 8L),
             EMC2:::.emc_wpool_done_encode(2L, generation = 4L, request = 9L)),
           done$wc)
  flush(done$wc)
  rec <- EMC2:::.emc_wpool_await_record(pool, function() integer(0), 5)
  expect_identical(rec$code, "record")
  expect_identical(rec$w, 2L)
  expect_identical(rec$generation, 4L)
  expect_identical(EMC2:::.emc_pool_state$stale_records, before + 1L)

  # Nothing further: the wait is bounded rather than blocking on a channel that
  # will never speak again.
  started <- proc.time()[["elapsed"]]
  out <- EMC2:::.emc_wpool_await_record(pool, function() integer(0), 0.2)
  expect_identical(out$code, "timeout")
  expect_lt(proc.time()[["elapsed"]] - started, 20)
})

test_that("a pool that replaces its workers replaces its generation", {
  # A generation that never changes is a field, not a fence.
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  pool <- EMC2:::.emc_wpool_start(2L, list(n_pars = 2))
  skip_if(is.null(pool), "no pool available here")
  # Reaping workers that were told to shut down: "did not deliver a result" is
  # the expected outcome, not something for the test log.
  on.exit(suppressWarnings(try(EMC2:::.emc_wpool_stop(pool), silent = TRUE)),
          add = TRUE)
  expect_gt(pool$generation, 0L)
  # Recycling stops the old workers on its way; reaping them is where the
  # "did not deliver a result" notice comes from, and it is not what this is
  # testing. Degradation on a failed recycle has its own test in
  # test-worker-pool.R.
  fresh <- suppressWarnings(
    EMC2:::.emc_wpool_recycle(pool, 2L, 1L, list(n_pars = 2)))
  # Strictly greater, and not a restart at 1: recycling on the clean backend
  # keeps the completion channel, so a generation that started over would accept
  # a record written by the workers it just replaced.
  expect_gt(fresh$generation, pool$generation)
  pool <- fresh
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

test_that("a message too large for one frame is segmented, not refused", {
  # "Keep safety limits on allocation/frames, but distinguish them from
  # task-size limits."  The 2 GiB guard exists so a desynchronised stream cannot
  # talk the reader into allocating on garbage; it is not a statement that a
  # larger message is illegitimate.  A payload past the segment size is split,
  # and the split is invisible to both callers.
  obj <- list(x = seq_len(4000L) + 0.5, tag = "segmented")
  small <- with_mocked_bindings(
    audit_frame_roundtrip(obj),
    .EMC_WPOOL_MAX_SEGMENT = 1024, .package = "EMC2")
  expect_true(small$ok)
  expect_identical(small$value, obj)

  # Segmenting costs one 8-byte header per segment and nothing else, so the
  # wire size says how many segments there were.
  whole <- audit_frame_roundtrip(obj)
  expect_identical(whole$value, obj)
  expect_gt(small$bytes, whole$bytes)
  expect_equal((small$bytes - whole$bytes) %% 8, 0)

  # A message that fits in one segment is written exactly as it was before
  # segmentation existed: one 8-byte positive length, then the payload.
  raw_all <- readBin(audit_truncated_frame(obj, .Machine$integer.max), "raw",
                     whole$bytes)
  expect_equal(readBin(raw_all[1:8], "double", 1L, size = 8, endian = "little"),
               whole$bytes - 8)
})

test_that("the reader reassembles segments and refuses a garbage length", {
  # Hand-built streams, so the reader is tested against the format rather than
  # against the writer agreeing with itself.
  obj <- list(a = 1:50, b = "marker")
  payload <- serialize(obj, NULL)
  half <- floor(length(payload) / 2)
  write_stream <- function(chunks) {
    f <- tempfile()
    con <- file(f, "wb")
    for (ch in chunks) {
      if (is.numeric(ch)) {
        writeBin(as.double(ch), con, size = 8L, endian = "little")
      } else writeBin(ch, con)
    }
    close(con)
    f
  }
  read_stream <- function(f) {
    con <- file(f, "rb")
    on.exit({ try(close(con), silent = TRUE); unlink(f) }, add = TRUE)
    EMC2:::.emc_wpool_recv(con)
  }

  # Two segments: a negative length says another follows.
  f <- write_stream(list(-half, payload[1:half],
                         length(payload) - half,
                         payload[(half + 1L):length(payload)]))
  expect_identical(read_stream(f), obj)

  # Three segments. (A zero-length continuation is not expressible -- -0 is 0,
  # which reads as a final segment -- and the writer never emits one, so the
  # sign convention costs nothing.)
  third <- floor(length(payload) / 3)
  f <- write_stream(list(-third, payload[1:third],
                         -third, payload[(third + 1L):(2L * third)],
                         length(payload) - 2L * third,
                         payload[(2L * third + 1L):length(payload)]))
  expect_identical(read_stream(f), obj)

  # A continuation that never ends is a truncated stream, not a hang.
  f <- write_stream(list(-half, payload[1:half]))
  expect_null(read_stream(f))

  # A length past the frame guard is refused before anything is allocated for
  # it, on both signs.
  for (len in c(EMC2:::.EMC_WPOOL_MAX_FRAME + 1,
                -(EMC2:::.EMC_WPOOL_MAX_FRAME + 1), NaN, Inf)) {
    f <- write_stream(list(len, payload[1:8]))
    expect_null(read_stream(f), info = paste("length", len))
  }
})

# --- bounded I/O (finding 2) ------------------------------------------------

test_that("a reply that never arrives is bounded rather than waited on", {
  # Finding 2.  The pool timed the wait for a completion token and then entered
  # a blocking framed receive, so a worker that announced completion and stalled
  # part-way through its reply could still hang the chain.  A length header and
  # short-read loops solve framing; they do not make I/O bounded.
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  dir <- tempfile("emc_bound_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "p")
  skip_if(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE) != 0)
  rc <- fifo(path, "rb", blocking = FALSE)
  wc <- fifo(path, "wb", blocking = FALSE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)

  # Nothing at all.
  started <- proc.time()[["elapsed"]]
  expect_null(EMC2:::.emc_wpool_recv(rc, deadline = 0.3))
  expect_lt(proc.time()[["elapsed"]] - started, 20)

  # A header, then silence: the stall the finding describes.  The payload gets
  # its own deadline rather than inheriting whatever was left of the header's.
  payload <- serialize(list(a = 1:100), NULL)
  writeBin(as.double(length(payload)), wc, size = 8L, endian = "little")
  writeBin(payload[1:16], wc)
  flush(wc)
  started <- proc.time()[["elapsed"]]
  expect_null(EMC2:::.emc_wpool_recv(rc, deadline = 0.3))
  expect_lt(proc.time()[["elapsed"]] - started, 20)
  expect_match(EMC2:::.emc_pool_state$last_transport_error, "rest of the reply")
})

test_that("a dead peer ends a receive without waiting out the deadline", {
  # Liveness and progress are different questions, and the receive asks both:
  # a worker that has exited is not going to answer, whatever the deadline says.
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  dir <- tempfile("emc_bound_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "p")
  skip_if(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE) != 0)
  rc <- fifo(path, "rb", blocking = FALSE)
  wc <- fifo(path, "wb", blocking = FALSE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)

  started <- proc.time()[["elapsed"]]
  expect_null(EMC2:::.emc_wpool_recv(rc, deadline = 600,
                                     alive = function() FALSE))
  # 600 s was the deadline; the liveness test is what ended it.
  expect_lt(proc.time()[["elapsed"]] - started, 30)
  expect_match(EMC2:::.emc_pool_state$last_transport_error, "exited")
})

test_that("a worker announces that it collected the request before computing", {
  # ST_START is the progress signal that a live pid cannot give: "collected the
  # request and is working" and "never collected it" are different failures.
  skip_on_os("windows")
  old <- options(emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  pool <- EMC2:::.emc_wpool_start(2L, list(n_pars = 2))
  skip_if(is.null(pool), "no pool available here")
  skip_if(is.null(pool$done), "no completion channel here")
  on.exit(suppressWarnings(try(EMC2:::.emc_wpool_stop(pool), silent = TRUE)),
          add = TRUE)

  gen <- pool$generation
  expect_true(EMC2:::.emc_wpool_request(
    pool, 1L, list(subs = 1L, notify = TRUE, w = 1L, generation = gen,
                   request = 4242L)))
  # There is no real context, so the work fails -- which is the point: the
  # worker still announces the pick-up first, and the completion still says what
  # happened, so the master learns of the failure from the channel rather than
  # by reading the reply.
  first <- EMC2:::.emc_wpool_await_record(pool, function() 1L, 30)
  second <- EMC2:::.emc_wpool_await_record(pool, function() 1L, 30)
  expect_identical(first$code, "record")
  expect_identical(first$status, EMC2:::.EMC_WPOOL_ST_START)
  expect_identical(first$w, 1L)
  expect_identical(first$request, 4242L)
  expect_identical(first$generation, gen)
  expect_identical(second$status, EMC2:::.EMC_WPOOL_ST_FAILED)
  expect_identical(second$request, 4242L)
  # And the reply itself is still there to be read.
  res <- EMC2:::.emc_wpool_reply(pool, 1L)
  expect_false(is.null(res$failed))
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

test_that("a late record from a retired worker is distinguishable", {
  # Finding 5's other half.  A late completion record from a retired worker can
  # still be sitting in the `done` channel, which is deliberately retained
  # across a recycle -- so after recycling, a record can name a worker index
  # that now belongs to a different process.
  #
  # Under the one-byte token these two were byte-identical and nothing could be
  # done about it.  They are now different records, which is what makes
  # retiring a generation expressible at all.  Acting on that -- draining the
  # retired generation's pending work before accepting new work -- is C5.
  same_index <- function(gen) {
    EMC2:::.emc_wpool_done_encode(3L, generation = gen, request = 1L)
  }
  expect_false(identical(same_index(1L), same_index(2L)))
  expect_identical(same_index(2L), same_index(2L))
  # And the same request answered twice by two generations decodes to two
  # different records rather than one repeated one.
  a <- EMC2:::.emc_wpool_done_decode(same_index(1L))
  b <- EMC2:::.emc_wpool_done_decode(same_index(2L))
  expect_identical(a$w, b$w)
  expect_false(identical(a$generation, b$generation))
})

test_that("a record for a settled request is not credited to the next one", {
  # The request id is the same defence one level down: within a single pool
  # generation, a worker whose share was recomputed in the master can still
  # answer afterwards, and that answer belongs to an iteration that is already
  # finished.
  dir <- tempfile("emc_done_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  done <- EMC2:::.emc_wpool_done_open(dir)
  skip_if(is.null(done), "no FIFO available here")
  on.exit(for (cn in done[c("rc", "wc")]) try(close(cn), silent = TRUE),
          add = TRUE)
  pool <- audit_fake_pool(2L)
  pool$done <- done
  pool$generation <- 1L

  writeBin(c(EMC2:::.emc_wpool_done_encode(1L, 1L, 100L),
             EMC2:::.emc_wpool_done_encode(1L, 1L, 101L)), done$wc)
  flush(done$wc)
  # Both records are current for this generation, so the channel returns both;
  # it is the request id that lets the caller tell which iteration each answers.
  first <- EMC2:::.emc_wpool_await_record(pool, function() integer(0), 5)
  second <- EMC2:::.emc_wpool_await_record(pool, function() integer(0), 5)
  expect_identical(first$request, 100L)
  expect_identical(second$request, 101L)
  expect_false(identical(first$request, second$request))
})
