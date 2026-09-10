# ===========================================================================
# Persistent worker pool for the particle step
# ===========================================================================
#
# Workers persist for an MCMC block and rebalance subject partitions each
# iteration. Named pipes carry shared block state and per-iteration work.
#
# The clean template forks workers before sampler history exists. Each subject
# carries its own L'Ecuyer stream, so results do not depend on worker count or
# partitioning; unavailable pools fall back to mcmapply.

# Store the last worker failure for diagnostics.
.emc_pool_state <- new.env(parent = emptyenv())

# Losing the pool preserves results and runs the remainder serially. Warn
# immediately because deferred warnings from forked workers may be lost.
.emc_wpool_degraded <- function(msg) {
  .emc_pool_state$last_error <- msg
  if (!isTRUE(.emc_pool_state$warned)) {
    .emc_pool_state$warned <- TRUE
    warning("EMC2 worker pool lost (", msg, "); the rest of this block runs ",
            "serially in the chain process. Results are unaffected, but it ",
            "will be slow.", call. = FALSE, immediate. = TRUE)
  }
  invisible(NULL)
}

# --- framed message transport ------------------------------------------------
#
# serialize()/unserialize() cannot be pointed straight at a FIFO.  A read from a
# pipe returns whatever happens to have arrived, and R treats that short read as
# a hard error -- "error reading from connection" -- rather than asking for the
# rest.  Sending a message as serialize(msg, con) was therefore a race the whole
# time: it emits a run of small writes, unserialize() issues a matching run of
# reads, and any read that outran the writer killed the worker in mid-message.
# The master, still writing into a pipe with no reader left, then took the
# SIGPIPE.  Both halves usually won the race, because a message below the 64 kB
# pipe buffer lands there before the worker is scheduled -- but a busy box loses
# it at 3 kB, and above the buffer the master has to interleave, so it stops
# being a race and starts being a certainty.  Two chunk sizes make it certain on
# their own: unserialize() asks for a whole XDR chunk at a time, 64,768 bytes
# for a double vector and 32,384 for an integer one, against reads that never
# return more than one 8 kB connection buffer.
#
# Messages are framed instead: an 8-byte little-endian payload length followed
# by the serialised bytes, both moved by loops that expect short reads and short
# writes.  Serialisation itself happens in memory, where it is exact, and the
# connection only ever sees raw bytes.  This is also about twice as fast as the
# old path on a typical iteration message, which made many more, smaller writes.
#
# Correctness does not rest on the framing being right.  A desynchronised stream
# fails the length sanity check or fails to unserialise, both of which return
# NULL, and every caller treats a NULL reply as "recompute this share here".

# Request size for a read; a pipe returns what it has, so this is only a ceiling.
.EMC_WPOOL_READ_CHUNK <- 65536L
# PIPE_BUF.  A blocking write this size completes in one go, so a partial write
# needs a signal to land in a window that is only a few microseconds wide.
.EMC_WPOOL_WRITE_CHUNK <- 4096L
# No pool message comes near this.  It exists so a desynchronised stream cannot
# talk the reader into allocating on a garbage length.  It is a safety limit,
# not a task-size limit, and the difference matters: a message that legitimately
# exceeds it is segmented rather than refused.
.EMC_WPOOL_MAX_FRAME <- 2^31
# Segment size for a payload larger than one frame.  A segment header carries a
# negative length to mean "this many bytes, and another segment follows"; the
# last segment's length is positive.  A message that fits in one segment is
# therefore byte-identical to what this sent before segmentation existed.
.EMC_WPOOL_MAX_SEGMENT <- 2^30

# Record why a transport call failed so the degradation warning can name it.
.emc_wpool_transport_fail <- function(msg) {
  .emc_pool_state$last_transport_error <- msg
  FALSE
}

.emc_wpool_transport_error <- function(default) {
  msg <- .emc_pool_state$last_transport_error
  if (is.null(msg) || !nzchar(msg)) default else msg
}

# writeBin() reports a short write as a warning and returns nothing, so promote
# that warning: a half-written frame is a dead channel, not a slow one.
#
# Every write is PIPE_BUF or smaller, which is what keeps a message from being
# interleaved with another writer's on a shared FIFO.
#
# On `deadline`, and what it does not do.  It is checked between chunks, so it
# bounds a transfer that is progressing too slowly -- a peer draining the pipe a
# few bytes at a time.  It cannot bound a write that is already blocked inside
# the kernel, because R offers no way to ask whether a connection is writable
# and no way to interrupt a write in progress; a peer that is alive and has
# stopped reading entirely will still hold the final chunk. That case is what
# the receive-side deadline and the worker's own liveness are for. Saying so
# here because "sends have deadlines" would otherwise read as a stronger
# guarantee than the platform can give.
.emc_wpool_put_bytes <- function(con, bytes, deadline = NULL) {
  n <- length(bytes)
  from <- 1
  started <- proc.time()[["elapsed"]]
  while (from <= n) {
    if (!is.null(deadline) && from > 1 &&
        proc.time()[["elapsed"]] - started > deadline) {
      stop(sprintf("the request was still being written after %.0f s", deadline),
           call. = FALSE)
    }
    to <- min(n, from + .EMC_WPOOL_WRITE_CHUNK - 1)
    withCallingHandlers(
      writeBin(bytes[from:to], con),
      warning = function(w) stop(conditionMessage(w), call. = FALSE))
    from <- to + 1
  }
  invisible(NULL)
}

# NULL means the bytes did not arrive.
#
# Two connection modes, and the difference is the whole reason this takes a
# deadline.  On a *blocking* connection an empty read is the pipe closing, and
# that is how both a clean shutdown and a dead worker present themselves.  On a
# *non-blocking* one an empty read means only "nothing yet", which is
# indistinguishable from end of stream -- R offers no way to tell them apart --
# so the question has to be put to the peer instead: is the process still
# alive, and has the deadline passed?  That is what makes a receive bounded.
# Waiting for a completion record first and then entering an unbounded blocking
# read, which is what this did before, bounds only the announcement.
#
# `deadline` is seconds of elapsed time, measured from entry, and is checked
# between reads. It bounds the wait for bytes that have not arrived; it cannot
# preempt a read already blocked inside the kernel, which is why the connection
# has to be non-blocking for the bound to mean anything.
.emc_wpool_get_bytes <- function(con, n, deadline = NULL, alive = NULL) {
  if (n <= 0) return(raw(0))
  out <- raw(n)
  got <- 0
  started <- proc.time()[["elapsed"]]
  wait <- 0
  idle <- 0L
  while (got < n) {
    chunk <- tryCatch(
      readBin(con, "raw", min(n - got, .EMC_WPOOL_READ_CHUNK)),
      error = function(e) {
        # R reports "nothing available yet" on a non-blocking connection as an
        # error rather than as a zero-length read, which is the same idiom the
        # completion channel already has to use.  With a deadline in force this
        # is therefore a retry; without one the connection is blocking and an
        # error is an error.
        # ... unless the connection itself has gone, in which case retrying
        # would spin out the whole deadline waiting on a handle that can never
        # produce anything.
        if (!is.null(deadline) &&
            isTRUE(tryCatch(isOpen(con), error = function(e2) FALSE))) {
          return(raw(0))
        }
        .emc_wpool_transport_fail(conditionMessage(e))
        NULL
      })
    if (is.null(chunk)) return(NULL)
    k <- length(chunk)
    if (k == 0L) {
      if (is.null(deadline)) {
        .emc_wpool_transport_fail(if (got == 0) {
          "the other end of the pipe closed"
        } else "the pipe closed part-way through a message")
        return(NULL)
      }
      idle <- idle + 1L
      # Liveness is cheaper to ask than to guess, but it is a /proc read: only
      # ask it once the spin has clearly failed.
      if ((idle %% 64L) == 0L && !is.null(alive) && !alive()) {
        .emc_wpool_transport_fail(if (got == 0) {
          "the worker exited without replying"
        } else "the worker exited part-way through a message")
        return(NULL)
      }
      if (proc.time()[["elapsed"]] - started > deadline) {
        .emc_wpool_transport_fail(sprintf(
          "no %s within %.0f s", if (got == 0) "reply" else "rest of the reply",
          deadline))
        return(NULL)
      }
      if (idle > 16L) { Sys.sleep(wait); wait <- min(max(wait * 2, 5e-5), 2e-3) }
      next
    }
    out[(got + 1):(got + k)] <- chunk
    got <- got + k
  }
  out
}

# Every write to a worker's request pipe goes through here, so the bound and
# the failure record are in one place rather than at each call site.
.emc_wpool_request <- function(pool, w, msg, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_send(pool$wcs[[w]], msg, deadline = deadline)
}

# Every read of a worker's reply pipe goes through here.  Those connections are
# non-blocking, so a bare `.emc_wpool_recv()` on one would read "nothing has
# arrived yet" as "the worker is gone" and retire a healthy pool on the first
# poll.  Routing them through one function is also what stops a future change
# of connection mode from being a silent behaviour change at six call sites.
.emc_wpool_reply <- function(pool, w, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_recv(pool$rcs[[w]], deadline = deadline,
                  alive = .emc_wpool_alive_fn(pool, w))
}

# Whether worker `w`'s process is still there. NULL when the pool has no job
# handle for it, which the read loop treats as "cannot tell" rather than dead.
.emc_wpool_alive_fn <- function(pool, w) {
  job <- if (!is.null(pool$jobs) && w <= length(pool$jobs)) pool$jobs[[w]] else NULL
  if (is.null(job)) return(NULL)
  function() {
    pid <- .emc_wpool_job_pid(job)
    is.null(pid) || is.na(pid) || .emc_wpool_pid_alive(pid)
  }
}

.emc_wpool_send <- function(con, obj, deadline = NULL) {
  payload <- serialize(obj, NULL)
  total <- length(payload)
  tryCatch({
    from <- 1
    repeat {
      to <- min(total, from + .EMC_WPOOL_MAX_SEGMENT - 1)
      last <- to >= total
      # Negative length means "another segment follows". A message that fits in
      # one segment is written exactly as it always was.
      len <- if (last) to - from + 1 else -(to - from + 1)
      .emc_wpool_put_bytes(con, writeBin(as.double(len), raw(), size = 8,
                                         endian = "little"), deadline)
      .emc_wpool_put_bytes(con, payload[from:to], deadline)
      if (last) break
      from <- to + 1
    }
    flush(con)
    TRUE
  }, error = function(e) .emc_wpool_transport_fail(conditionMessage(e)))
}

.emc_wpool_recv <- function(con, deadline = NULL, alive = NULL) {
  parts <- list()
  repeat {
    header <- .emc_wpool_get_bytes(con, 8, deadline, alive)
    if (is.null(header)) return(NULL)
    n <- readBin(header, "double", 1L, size = 8, endian = "little")
    more <- is.finite(n) && n < 0
    if (more) n <- -n
    # Validated before anything is allocated for it: a desynchronised stream
    # must not be able to talk the reader into reserving a garbage length.
    if (!is.finite(n) || n < 0 || n > .EMC_WPOOL_MAX_FRAME) {
      .emc_wpool_transport_fail("implausible message length on the pipe")
      return(NULL)
    }
    # The header arrived, so the payload is owed to us: it gets its own
    # deadline rather than the remainder of the header's.
    payload <- .emc_wpool_get_bytes(con, n, deadline, alive)
    if (is.null(payload)) return(NULL)
    if (!length(parts) && !more) {
      parts <- list(payload)
      break
    }
    parts[[length(parts) + 1L]] <- payload
    if (!more) break
  }
  payload <- if (length(parts) == 1L) parts[[1L]] else unlist(parts, use.names = FALSE)
  tryCatch(unserialize(payload), error = function(e) {
    .emc_wpool_transport_fail(conditionMessage(e))
    NULL
  })
}

# The shared group covariance is identical for every worker in an iteration.
# Keep its serialised bytes in one file and put only the path in each request.
# Write to a private temporary name and rename it into place: workers can only
# observe a complete file, even on a filesystem where opening a path races with
# the writer.  The caller owns the final path and removes it after all replies
# (or serial fallbacks) have finished using the decoded object.
.emc_wpool_shared_write <- function(dir, bytes) {
  if (!is.character(dir) || length(dir) != 1L || !dir.exists(dir)) {
    .emc_wpool_transport_fail("shared broadcast directory is unavailable")
    return(NULL)
  }
  if (!is.raw(bytes)) {
    .emc_wpool_transport_fail("shared broadcast payload is not raw bytes")
    return(NULL)
  }
  if (length(bytes) > .EMC_WPOOL_MAX_FRAME ||
      length(bytes) > .Machine$integer.max) {
    .emc_wpool_transport_fail("shared broadcast payload is too large")
    return(NULL)
  }
  path <- tempfile(pattern = "shared_", tmpdir = dir, fileext = ".bin")
  tmp <- paste0(path, ".tmp")
  con <- NULL
  complete <- FALSE
  on.exit({
    if (!is.null(con)) try(close(con), silent = TRUE)
    if (!complete) unlink(c(path, tmp))
  }, add = TRUE)
  tryCatch({
    con <- file(tmp, open = "wb")
    withCallingHandlers(
      writeBin(bytes, con),
      warning = function(w) stop(conditionMessage(w), call. = FALSE))
    flush(con)
    close(con)
    con <- NULL
    if (!file.rename(tmp, path)) {
      stop("could not publish shared broadcast file", call. = FALSE)
    }
    complete <- TRUE
    path
  }, error = function(e) {
    .emc_wpool_transport_fail(conditionMessage(e))
    NULL
  })
}

# Read a published shared payload and validate the complete file before
# unserialising it.  A missing/truncated file is a worker-local failure; the
# surrounding worker handler reports it and the master recomputes that share
# from its already-decoded cache.
.emc_wpool_shared_read <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("invalid shared broadcast path", call. = FALSE)
  }
  size <- file.info(path)$size
  if (length(size) != 1L || is.na(size) || !is.finite(size) || size < 0 ||
      size > .EMC_WPOOL_MAX_FRAME || size > .Machine$integer.max) {
    stop("invalid shared broadcast file", call. = FALSE)
  }
  con <- file(path, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  payload <- .emc_wpool_get_bytes(con, as.integer(size))
  if (is.null(payload) || length(payload) != as.integer(size)) {
    stop("shared broadcast file was truncated", call. = FALSE)
  }
  tryCatch(unserialize(payload), error = function(e) {
    stop(conditionMessage(e), call. = FALSE)
  })
}

# How long to wait for a worker's reply before treating it as wedged.
#
# The pool's recovery paths all key on a worker EXITING: a dead worker closes
# its pipe, the read returns zero bytes and the master recomputes that share
# serially.  A worker that is alive and stuck emits no such signal, and the
# master's blocking read then waits forever -- a fit that looks like it is
# running, with a frozen progress bar and no indication of which subject is
# responsible.  That is not hypothetical: a single PDE solve diverging at the
# edge of the RLF/FPE parameter space has wedged a sampler before, and those
# solvers sit on this path.
#
# The deadline is deliberately loose: 20x the median observed receive time,
# floored at two minutes, so a merely slow iteration can never trip it.
# options(emc2.worker_timeout = <seconds>) overrides it; Inf disables it.
.EMC_WPOOL_TIMEOUT_FLOOR <- 120
.EMC_WPOOL_TIMEOUT_COLD <- 600
# Starting a clean R process, loading the namespace and forking workers from it.
# Generous, because a cold library on a busy filesystem is genuinely slow, but
# bounded: this used to be a blocking read with nothing behind it.
.EMC_WPOOL_TIMEOUT_TEMPLATE <- 300
# One-shot wrapper jobs run whole model comparisons, not one iteration, so the
# sampler's history says nothing about how long they should take.  A safety
# limit is not a task-size limit: this exists so a dead wrapper cannot hang the
# session, not to decide that an hour of legitimate work is too much.
.EMC_WPOOL_TIMEOUT_LOO <- 24 * 3600

# `adapt = FALSE` for work the sampler's own receive history says nothing
# about: scaling a model comparison by the median particle-step time would be
# using one task's size to bound another's.
.emc_wpool_deadline <- function(cold = .EMC_WPOOL_TIMEOUT_COLD, adapt = TRUE) {
  opt <- getOption("emc2.worker_timeout", NA_real_)
  if (!is.na(opt)) return(as.numeric(opt))
  h <- .emc_pool_state$iter_receive_times
  if (!isTRUE(adapt) || is.null(h) || !length(h)) return(cold)
  max(.EMC_WPOOL_TIMEOUT_FLOOR, 20 * stats::median(h))
}

.emc_wpool_record_receive <- function(seconds) {
  h <- c(.emc_pool_state$iter_receive_times, seconds)
  .emc_pool_state$iter_receive_times <- utils::tail(h, 32L)
  invisible(NULL)
}

# Wait for the next completion record on the shared non-blocking `done`
# channel.  Returns one of
#
#   list(code = "record", w, generation, request, status)
#   list(code = "dead",  w)      a worker we are still waiting on has died
#   list(code = "timeout")       `deadline` seconds with every pending worker
#                                alive and none of them announcing anything
#
# `pending` is a function rather than a vector because the wait loops here: the
# set of workers still outstanding is decided by the caller and can only be
# read live.
#
# Records from another pool generation are dropped rather than acted on.  They
# are the late replies of a pool that has already been retired, and before the
# generation field existed they were indistinguishable from a current worker's
# -- which is how a retired worker could answer for a live one.  Counting them
# is deliberate: retiring the generation they belong to is C5's work, and a
# count is what shows whether it is needed.
.emc_wpool_await_record <- function(pool, pending, deadline = NULL) {
  started <- proc.time()[["elapsed"]]
  idle <- 0L
  wait <- 0
  repeat {
    rec <- .emc_wpool_done_take(pool$done)
    if (!is.null(rec)) {
      if (identical(rec$generation, .emc_wpool_one_int(pool$generation))) {
        return(rec)
      }
      .emc_pool_state$stale_records <-
        (if (is.null(.emc_pool_state$stale_records)) 0L
         else .emc_pool_state$stale_records) + 1L
      next
    }
    idle <- idle + 1L
    if ((idle %% 64L) == 0L) {
      for (w in pending()) {
        pid <- .emc_wpool_job_pid(pool$jobs[[w]])
        if (!is.null(pid) && !.emc_wpool_pid_alive(pid)) {
          return(list(code = "dead", w = w))
        }
      }
      if (!is.null(deadline) &&
          proc.time()[["elapsed"]] - started > deadline) {
        return(list(code = "timeout"))
      }
    }
    if (idle > 16L) { Sys.sleep(wait); wait <- min(max(wait * 2, 5e-5), 2e-3) }
  }
}

# Spawn requires an installed namespace and serialisable context. Development
# namespaces and custom trend kernels use fork, where external pointers survive.
.emc_wpool_backend <- function(ctx) {
  requested <- match.arg(getOption("emc2.worker_backend", "spawn"),
                         c("spawn", "fork"))
  if (requested == "fork") return("fork")
  ns_path <- tryCatch(getNamespaceInfo(asNamespace("EMC2"), "path"),
                      error = function(e) "")
  installed <- nzchar(ns_path) &&
    file.exists(file.path(ns_path, "Meta", "package.rds"))
  rscript <- file.path(R.home("bin"), "Rscript")
  if (!installed || !file.exists(rscript) || .emc_wpool_has_custom_ptr(ctx)) {
    return("fork")
  }
  "spawn"
}

.emc_wpool_has_custom_ptr <- function(ctx) {
  has_ptr <- function(x) {
    if (typeof(x) == "externalptr") return(TRUE)
    values <- if (is.list(x) || is.pairlist(x)) as.list(x) else list()
    attrs <- attributes(x)
    any(vapply(c(values, if (is.null(attrs)) list() else attrs), has_ptr,
               logical(1)))
  }
  model <- ctx$model
  spec <- if (is.function(model)) {
    tryCatch(model(), error = function(e) NULL)
  } else if (is.list(model) && length(model) &&
             all(vapply(model, is.function, logical(1)))) {
    lapply(model, function(x) tryCatch(x(), error = function(e) NULL))
  } else {
    model
  }
  has_ptr(spec)
}

# Reclose each model around only its realised specification so the clean
# template receives no unused constructor state. The likelihood-facing contract
# remains a callable returning the model list.
.emc_wpool_slim_model <- function(model) {
  if (is.function(model)) {
    spec <- model()
    out <- function() NULL
    body(out) <- quote(spec)
    environment(out) <- list2env(list(spec = spec), parent = emptyenv())
    return(out)
  }
  if (is.list(model)) {
    out <- lapply(model, .emc_wpool_slim_model)
    attributes(out) <- attributes(model)
    return(out)
  }
  model
}

.emc_wpool_pid_alive <- function(pid) {
  if (length(pid) != 1L || is.na(pid)) return(FALSE)
  alive <- isTRUE(tryCatch(tools::pskill(pid, 0L),
                           error = function(e) FALSE))
  if (!alive) return(FALSE)
  # kill(pid, 0) reports TRUE for a zombie retained by the clean template.
  # Treat that process as dead so a missing completion token triggers the
  # serial correctness fallback instead of an infinite FIFO wait.
  # A process that exits between the signal probe and this read is normal, and
  # readLines() warns as well as errors when the file has gone.
  stat <- suppressWarnings(
    tryCatch(readLines(sprintf("/proc/%d/stat", as.integer(pid)), n = 1L),
             error = function(e) character()))
  if (length(stat)) {
    tail <- sub("^.*\\) ", "", stat)
    state <- strsplit(tail, " ", fixed = TRUE)[[1L]][1L]
    if (identical(state, "Z")) return(FALSE)
  }
  TRUE
}

.emc_wpool_job_pid <- function(job) {
  if (!is.null(job$pid) && !is.na(job$pid)) return(as.integer(job$pid))
  if (!is.null(job$pidfile) && file.exists(job$pidfile)) {
    pid <- suppressWarnings(as.integer(tryCatch(readLines(job$pidfile, n = 1L),
                                                error = function(e) NA_character_)))
    if (length(pid) == 1L && !is.na(pid)) return(pid)
  }
  NA_integer_
}

.emc_wpool_terminate_jobs <- function(jobs, wait = TRUE, terminate = TRUE) {
  if (!length(jobs)) return(invisible(NULL))
  external <- vapply(jobs, function(x) isTRUE(x$external), logical(1))
  forked <- jobs[!external]
  spawned <- jobs[external]
  if (length(forked)) {
    if (terminate) {
      for (job in forked) try(tools::pskill(job$pid), silent = TRUE)
    }
    try(parallel::mccollect(forked, wait = wait), silent = TRUE)
  }
  if (length(spawned)) {
    pids <- vapply(spawned, .emc_wpool_job_pid, integer(1))
    if (wait) {
      # proc.time(), not Sys.time(): a deadline measured on a clock the
      # administrator can set backwards is not a deadline.
      deadline <- proc.time()[["elapsed"]] + 2
      while (any(vapply(pids, .emc_wpool_pid_alive, logical(1))) &&
             proc.time()[["elapsed"]] < deadline) Sys.sleep(0.005)
    }
    live <- pids[vapply(pids, .emc_wpool_pid_alive, logical(1))]
    if (terminate || wait) {
      for (pid in live) try(tools::pskill(pid), silent = TRUE)
    }
  }
  invisible(NULL)
}

# --- pool lifecycle ---------------------------------------------------------

# Start workers for `idx` and connect to their named pipes. The request FIFO is
# opened non-blocking for the handshake, then reopened blocking for payloads:
# large serialised messages may exceed the pipe buffer. Replies are read before
# the next send, so blocking writes cannot deadlock.
.emc_wpool_spawn <- function(dir, idx, ctx, backend = "fork", template = NULL) {
  req <- file.path(dir, sprintf("req%d", idx))
  ans <- file.path(dir, sprintf("ans%d", idx))
  made <- tryCatch(
    system2("mkfifo", shQuote(c(req, ans)), stdout = FALSE, stderr = FALSE),
    warning = function(e) 1L, error = function(e) 1L
  )
  if (!identical(as.integer(made), 0L)) return(NULL)

  n <- length(idx)
  jobs <- vector("list", n)
  wcs <- vector("list", n)
  rcs <- vector("list", n)
  success <- FALSE
  on.exit({
    if (!success) {
      for (cn in c(wcs, rcs)) if (!is.null(cn)) try(close(cn), silent = TRUE)
      live <- jobs[!vapply(jobs, is.null, logical(1))]
      .emc_wpool_terminate_jobs(live, wait = TRUE)
      # Remove partial FIFOs so a later grow can recreate the worker indices.
      unlink(c(req, ans))
    }
  }, add = TRUE)

  # Process creation can fail under resource limits; return NULL so callers
  # preserve the serial fallback rather than aborting the chain.
  started <- tryCatch({
    if (identical(backend, "spawn")) {
      if (is.null(template) || !isTRUE(template$alive)) {
        stop("clean worker template is not available")
      }
      if (!.emc_wpool_send(template$wc, list(command = "spawn", idx = idx,
                                             req = req, ans = ans))) {
        stop(.emc_wpool_transport_error("template stopped accepting requests"))
      }
      reply <- .emc_wpool_recv(
        template$rc, deadline = .EMC_WPOOL_TIMEOUT_TEMPLATE,
        alive = local({ j <- template$job; function() {
          pid <- .emc_wpool_job_pid(j)
          is.null(pid) || is.na(pid) || .emc_wpool_pid_alive(pid)
        } }))
      if (!isTRUE(reply$ok) || length(reply$pids) != n) {
        stop(if (is.null(reply$error)) "template could not fork workers" else reply$error)
      }
      jobs <- lapply(as.integer(reply$pids), function(pid) {
        list(pid = pid, external = TRUE)
      })
    } else {
      for (w in seq_len(n)) {
        jobs[[w]] <- local({
          i <- w
          parallel::mcparallel(.emc_wpool_serve(req[i], ans[i], ctx),
                               detached = FALSE)
        })
      }
    }
    TRUE
  }, error = function(e) FALSE)
  if (!started) return(NULL)

  tryCatch({
    for (w in seq_len(n)) {
      # Poll for the worker's FIFO with bounded exponential backoff.
      f <- NULL
      deadline <- proc.time()[["elapsed"]] + 20  # startup timeout ceiling
      wait <- 0.0002
      repeat {
        if (identical(backend, "spawn")) {
          if (!.emc_wpool_pid_alive(jobs[[w]]$pid)) {
            stop("Worker process died before pool initialization")
          }
        } else {
          res <- parallel::mccollect(jobs[[w]], wait = FALSE, timeout = 0)
          if (!is.null(res)) stop("Worker process died before pool initialization")
        }
        f <- tryCatch(suppressWarnings(fifo(req[w], "wb", blocking = FALSE)),
                      error = function(e) NULL)
        if (!is.null(f)) break
        if (proc.time()[["elapsed"]] > deadline) break
        Sys.sleep(wait)
        wait <- min(wait * 2, 0.01)
      }
      if (is.null(f)) stop("Timeout waiting for worker process to initialize")
      # Replace the non-blocking handshake handle with the blocking transport.
      blocking <- fifo(req[w], "wb", blocking = TRUE)
      close(f)
      wcs[[w]] <- blocking
    }
    for (w in seq_len(n)) {
      # Non-blocking, because a deadline on a blocking read is decoration: R
      # cannot interrupt a read that is already inside the kernel.  Every read
      # of these connections therefore passes a deadline and a liveness test;
      # see .emc_wpool_get_bytes().
      rcs[[w]] <- fifo(ans[w], "rb", blocking = FALSE)
    }
    success <- TRUE
  }, error = function(e) NULL)

  if (!success) return(NULL)
  list(jobs = jobs, wcs = wcs, rcs = rcs)
}

# Start a pristine process that loads package/context once, then forks the
# persistent workers. It provides a clean parent for recycling without
# reserialising or reloading the context.
.emc_wpool_template_start <- function(dir, ctx_file, lib) {
  req <- file.path(dir, "template_req")
  ans <- file.path(dir, "template_ans")
  made <- tryCatch(
    system2("mkfifo", shQuote(c(req, ans)), stdout = FALSE, stderr = FALSE),
    warning = function(e) 1L, error = function(e) 1L
  )
  if (!identical(as.integer(made), 0L)) return(NULL)

  boot <- file.path(dir, "template.rds")
  pidfile <- file.path(dir, "template.pid")
  log <- file.path(dir, "template.log")
  saveRDS(list(req = req, ans = ans, ctx = ctx_file, lib = lib,
               pidfile = pidfile), boot, compress = FALSE)
  # Give the template and its children a recognisable process title for
  # diagnostics; sanitise the directory tag before embedding it in the command.
  tag <- sprintf("EMC2-worker-pool[%s|master=%d]",
                 gsub("[^A-Za-z0-9._-]", "_", basename(dir)), Sys.getpid())
  expr <- paste0(
    "EMC2.worker<-'", tag, "';",
    "b<-readRDS(commandArgs(TRUE)[1L]);",
    "writeLines(as.character(Sys.getpid()),b$pidfile);",
    ".libPaths(unique(c(b$lib,.libPaths())));",
    "ns<-loadNamespace('EMC2');",
    "get('.emc_wpool_template',ns)(b)"
  )
  status <- tryCatch(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(expr), shQuote(boot)),
    stdout = log, stderr = log, wait = FALSE
  ), error = function(e) 1L)
  if (!identical(as.integer(status), 0L)) {
    unlink(c(req, ans, boot, pidfile, log))
    return(NULL)
  }
  job <- list(pid = NA_integer_, pidfile = pidfile, log = log, external = TRUE)
  wc <- rc <- NULL
  success <- FALSE
  on.exit({
    if (!success) {
      for (cn in list(wc, rc)) if (!is.null(cn)) try(close(cn), silent = TRUE)
      .emc_wpool_terminate_jobs(list(job), wait = TRUE)
      unlink(c(req, ans, boot, pidfile, log))
    }
  }, add = TRUE)

  f <- NULL
  deadline <- proc.time()[["elapsed"]] + 20
  wait <- 0.0002
  repeat {
    pid <- .emc_wpool_job_pid(job)
    if (!is.na(pid)) {
      job$pid <- pid
      if (!.emc_wpool_pid_alive(pid)) break
    }
    f <- tryCatch(suppressWarnings(fifo(req, "wb", blocking = FALSE)),
                  error = function(e) NULL)
    if (!is.null(f) || proc.time()[["elapsed"]] > deadline) break
    Sys.sleep(wait)
    wait <- min(wait * 2, 0.01)
  }
  if (is.null(f)) return(NULL)
  wc <- fifo(req, "wb", blocking = TRUE)
  close(f)
  # Non-blocking, like a worker's reply pipe and for the same reason: every
  # receive on it is bounded, so a template that starts and then wedges cannot
  # hang the pool it was started for.  Growing and recycling both go through
  # here, so an unbounded read is a hang in the middle of a fit, not only at
  # start-up.
  rc <- fifo(ans, "rb", blocking = FALSE)
  alive <- local({ j <- job; function() {
    pid <- .emc_wpool_job_pid(j)
    is.null(pid) || is.na(pid) || .emc_wpool_pid_alive(pid)
  } })
  ready <- .emc_wpool_recv(rc, deadline = .EMC_WPOOL_TIMEOUT_TEMPLATE,
                           alive = alive)
  if (!isTRUE(ready$ready)) return(NULL)
  job$pid <- as.integer(ready$pid)
  success <- TRUE
  list(wc = wc, rc = rc, job = job, alive = TRUE)
}

# `ctx` is constant for the block. Spawn serialises it once to the template;
# fork workers inherit it directly. `allow_one` is used by the one-shot LOO
# wrapper, which needs a single clean process around loo's own inner workers;
# the normal sampler/comparison pool still refuses a one-worker pool.
.emc_wpool_start <- function(n_workers, ctx, allow_one = FALSE) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L ||
      (!isTRUE(allow_one) && n_workers <= 1L)) return(NULL)
  if (Sys.info()[["sysname"]] == "Windows") return(NULL)
  if (!nzchar(Sys.which("mkfifo"))) return(NULL)
  # After the platform checks, so a machine with no pool at all does not report
  # a connection problem it does not have -- and before anything is created, so
  # a width the connection table cannot hold is never half-built.
  n_workers <- .emc_wpool_feasible_workers(n_workers)
  if (n_workers < 1L || (!isTRUE(allow_one) && n_workers <= 1L)) return(NULL)

  dir <- tempfile("emc_wpool_")
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) return(NULL)

  backend <- .emc_wpool_backend(ctx)
  ctx_file <- NULL
  lib <- NULL
  if (identical(backend, "spawn")) {
    ctx_file <- file.path(dir, "context.rds")
    saved <- tryCatch({ saveRDS(ctx, ctx_file, compress = FALSE); TRUE },
                      error = function(e) FALSE)
    if (!saved) backend <- "fork"
    if (identical(backend, "spawn")) {
      ns_path <- getNamespaceInfo(asNamespace("EMC2"), "path")
      lib <- dirname(ns_path)
    }
  }

  template <- NULL
  if (identical(backend, "spawn")) {
    template <- .emc_wpool_template_start(dir, ctx_file, lib)
    if (is.null(template)) backend <- "fork"
  }

  # Per pool, not per session: a later block that degrades must say so too.
  .emc_pool_state$warned <- FALSE
  spawned <- .emc_wpool_spawn(dir, seq_len(n_workers), ctx, backend,
                              template = template)
  if (is.null(spawned)) {
    if (!is.null(template)) .emc_wpool_template_stop(template)
    unlink(dir, recursive = TRUE)
    return(NULL)
  }

  # Generation 1. Every replacement of these workers increments it, and every
  # request and completion record carries it, so a record written by a pool
  # that has already been retired is recognisable instead of plausible.
  list(n = n_workers, dir = dir, jobs = spawned$jobs, wcs = spawned$wcs,
       rcs = spawned$rcs, alive = TRUE, backend = backend,
       generation = .emc_wpool_next_generation(),
       ctx_file = ctx_file, lib = lib, template = template,
       done = .emc_wpool_done_open(dir))
}

# Completion FIFO shared by every worker. A worker announces itself on it
# before writing its (possibly large) reply, which is what lets the master wait
# for whoever finished rather than blocking on one worker's pipe.
#
# The announcement used to be one byte holding the worker index. One byte says
# only "it is me", and it cannot even say that past 255: `as.raw(256)` is 0 with
# a warning, so a 256-worker pool would announce worker 0 and the master would
# read it as lost alignment and retire a healthy pool. The dynamic queue capped
# itself at 255 workers to avoid the question; the subject pool never did.
#
# The record below answers four questions instead of one, and each of the other
# three is load-bearing:
#
#   worker      which pipe holds the reply, in a width that does not run out
#   generation  which pool this belongs to, so a reply from a pool that has
#               already been retired is recognisable rather than plausible
#   request     which iteration it answers, for the same reason
#   status      what happened, so the master knows a failure is coming before
#               it commits to reading the reply -- and, as ST_START, that the
#               request was picked up at all
#
# ST_START is progress; a live pid is not. They are different questions: a
# worker that never collected its request and a worker wedged inside the
# likelihood are both alive, and only the first is worth retrying elsewhere.
#
# Sixteen bytes, little-endian, fixed width:
#   1-4    worker id       int32, 1-based
#   5-8    pool generation int32
#   9-12   request id      int32
#   13     status          raw
#   14-16  reserved, zero
#
# That is far below PIPE_BUF, so a blocking write of a whole record cannot be
# interleaved with another worker's. The reader accumulates anyway: a
# non-blocking read is permitted to return less than was written, and a
# protocol that depends on it not doing so is one that works until it does not.
.EMC_WPOOL_DONE_BYTES <- 16L
.EMC_WPOOL_ST_START <- 1L
.EMC_WPOOL_ST_OK <- 2L
.EMC_WPOOL_ST_FAILED <- 3L

# A request that predates these fields, or a corrupted one, encodes as 0
# rather than shortening the record: a record of the wrong length would
# desynchronise the channel for every worker sharing it.
.emc_wpool_one_int <- function(x) {
  x <- suppressWarnings(as.integer(x))
  if (!length(x) || is.na(x[[1L]])) 0L else x[[1L]]
}

.emc_wpool_done_encode <- function(w, generation = 0L, request = 0L,
                                   status = .EMC_WPOOL_ST_OK) {
  out <- raw(.EMC_WPOOL_DONE_BYTES)
  out[1:12] <- writeBin(vapply(list(w, generation, request),
                               .emc_wpool_one_int, integer(1)),
                        raw(), size = 4L, endian = "little")
  out[13] <- as.raw(.emc_wpool_one_int(status) %% 256L)
  out
}

# Monotonic within this process, across every pool it ever builds.  Taking a
# fresh generation per pool rather than restarting at 1 is what makes the fence
# work at the one place it matters: recycling on the clean backend keeps the
# directory and the completion channel, so a record written by the workers that
# were just replaced is still in the pipe.  A pool that started again at 1
# would accept it.
.emc_wpool_next_generation <- function() {
  n <- .emc_pool_state$generation_seq
  n <- if (is.null(n)) 1L else n + 1L
  if (n >= .Machine$integer.max) n <- 1L
  .emc_pool_state$generation_seq <- n
  n
}

# Monotonic within this process, which is the right scope: the master and the
# workers it owns are the only readers, and the generation distinguishes pools
# after a wrap.
.emc_wpool_next_request <- function() {
  n <- .emc_pool_state$request_seq
  n <- if (is.null(n)) 1L else n + 1L
  if (n >= .Machine$integer.max) n <- 1L
  .emc_pool_state$request_seq <- n
  n
}

.emc_wpool_done_decode <- function(bytes) {
  ints <- readBin(bytes[1:12], "integer", 3L, size = 4L, endian = "little")
  list(code = "record", w = ints[[1L]], generation = ints[[2L]],
       request = ints[[3L]], status = as.integer(bytes[[13L]]))
}

# One complete record, or NULL if nothing whole has arrived yet. Whatever is
# available is drained into the pool's buffer on every call, so several records
# that arrived together cost one read between them.
.emc_wpool_done_take <- function(done) {
  buf <- done$buf
  chunk <- tryCatch(readBin(done$rc, "raw", .EMC_WPOOL_READ_CHUNK),
                    error = function(e) raw(0))
  if (length(chunk)) buf$bytes <- c(buf$bytes, chunk)
  if (length(buf$bytes) < .EMC_WPOOL_DONE_BYTES) return(NULL)
  rec <- .emc_wpool_done_decode(buf$bytes[seq_len(.EMC_WPOOL_DONE_BYTES)])
  buf$bytes <- buf$bytes[-seq_len(.EMC_WPOOL_DONE_BYTES)]
  rec
}

# The reader is non-blocking to detect worker death; keeping a master writer
# open prevents spurious EOF. NULL disables notification and selects the
# static split.
.emc_wpool_done_open <- function(dir) {
  path <- file.path(dir, "done")
  made <- tryCatch(
    system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE),
    warning = function(e) 1L, error = function(e) 1L
  )
  if (!identical(as.integer(made), 0L)) return(NULL)
  rc <- tryCatch(fifo(path, "rb", blocking = FALSE), error = function(e) NULL)
  if (is.null(rc)) return(NULL)
  # Opening for write needs the reader to exist already, which it now does.
  wc <- tryCatch(fifo(path, "wb", blocking = FALSE), error = function(e) NULL)
  if (is.null(wc)) { try(close(rc), silent = TRUE); return(NULL) }
  # Partial records live here between polls. An environment because the pool is
  # copied by value into every helper that reads from it.
  buf <- new.env(parent = emptyenv())
  buf$bytes <- raw(0)
  list(rc = rc, wc = wc, buf = buf)
}

# --- ownership and retirement -----------------------------------------------

# Everything the retired generation left on the completion channel.
#
# The channel is deliberately retained across a recycle, so it is the one place
# where a record written by workers that no longer exist can outlive them.  The
# generation field lets the master recognise such a record; retiring the
# generation is throwing them away before the replacement workers are given any
# work, so that "records waiting" and "this iteration's workers have answered"
# cannot be confused.  A partial record goes too: whoever was part-way through
# writing it is gone.
.emc_wpool_retire <- function(pool) {
  if (is.null(pool) || is.null(pool$done)) return(0L)
  n <- 0L
  repeat {
    if (is.null(.emc_wpool_done_take(pool$done))) break
    n <- n + 1L
  }
  pool$done$buf$bytes <- raw(0)
  if (n) {
    .emc_pool_state$retired_records <-
      (if (is.null(.emc_pool_state$retired_records)) 0L
       else .emc_pool_state$retired_records) + n
  }
  n
}

# R's connection table is a fixed number of slots, and this pool takes two per
# worker plus the completion channel and the template handshake.  A wide core
# budget can therefore exhaust it -- 64 workers is 132 connections against a
# default limit of 128 -- and the failure would land part-way through startup,
# after some workers had already been forked.
#
# Allocate a feasible width instead and let the partition give the remaining
# subjects to the workers that do exist.  That is the same work, queued.
.EMC_WPOOL_CONN_PER_WORKER <- 2L
.EMC_WPOOL_CONN_OVERHEAD <- 6L   # completion channel (2), template (2), slack (2)

.emc_wpool_conn_limit <- function() {
  limit <- suppressWarnings(as.integer(
    Sys.getenv("R_MAX_NUM_CONNECTIONS", "128")))
  if (is.na(limit) || limit <= 0L) 128L else limit
}

.emc_wpool_feasible_workers <- function(n_workers) {
  n_workers <- as.integer(n_workers)
  # One is checked too: growing by one repeatedly reaches the same limit as
  # asking for the whole width at once.
  if (is.na(n_workers) || n_workers < 1L) return(n_workers)
  open <- tryCatch(nrow(showConnections(all = TRUE)),
                   error = function(e) NA_integer_)
  if (is.na(open)) return(n_workers)
  room <- .emc_wpool_conn_limit() - open - .EMC_WPOOL_CONN_OVERHEAD
  feasible <- room %/% .EMC_WPOOL_CONN_PER_WORKER
  if (feasible >= n_workers) return(n_workers)
  feasible <- max(0L, as.integer(feasible))
  .emc_pool_state$conn_capped <- TRUE
  if (!isTRUE(.emc_pool_state$conn_warned)) {
    .emc_pool_state$conn_warned <- TRUE
    warning(sprintf(
      paste0("EMC2 worker pool limited to %d workers instead of %d: R allows ",
             "%d open connections and this pool needs %d per worker. Results ",
             "are unaffected. Raise R_MAX_NUM_CONNECTIONS to use more."),
      feasible, n_workers, .emc_wpool_conn_limit(),
      .EMC_WPOOL_CONN_PER_WORKER), call. = FALSE, immediate. = TRUE)
  }
  feasible
}

# Grow the pool into cores released by completed chains. Per-subject streams make
# changing worker count and assignment deterministic.
.emc_wpool_grow <- function(pool, n_total, ctx) {
  n_total <- as.integer(n_total)
  if (is.null(pool) || !isTRUE(pool$alive) || is.na(n_total)) return(pool)
  n_total <- min(n_total, pool$n + .emc_wpool_feasible_workers(n_total - pool$n))
  n_new <- n_total - pool$n
  if (n_new <= 0L) return(pool)

  # Do not retry a failed grow every iteration: handshake failures can be slow.
  # Exponential backoff leaves the current pool width in place.
  if (isTRUE(pool$grow_skip > 0L)) {
    pool$grow_skip <- pool$grow_skip - 1L
    return(pool)
  }

  spawned <- .emc_wpool_spawn(pool$dir, pool$n + seq_len(n_new), ctx,
                              backend = if (is.null(pool$backend)) "fork" else pool$backend,
                              template = pool$template)
  if (is.null(spawned)) {
    fails <- if (is.null(pool$grow_fails)) 1L else pool$grow_fails + 1L
    pool$grow_fails <- fails
    pool$grow_skip <- bitwShiftL(1L, min(fails, 6L))  # exponential skip interval
    return(pool)
  }

  # Keep the existing completion FIFO so added workers join the same queue --
  # and, because it is kept, empty it: this is a new generation, and nothing
  # written under the old one answers work that has not been sent yet.
  .emc_wpool_retire(pool)
  list(n = n_total, dir = pool$dir, jobs = c(pool$jobs, spawned$jobs),
       wcs = c(pool$wcs, spawned$wcs), rcs = c(pool$rcs, spawned$rcs),
       alive = TRUE, grow_fails = 0L, grow_skip = 0L,
       generation = .emc_wpool_next_generation(),
       backend = pool$backend, ctx_file = pool$ctx_file, lib = pool$lib,
       template = pool$template, done = pool$done)
}

# Replace workers with a fresh fork every `every` iterations. Recycling limits
# copy-on-write growth from dirty block-context pages; the pristine template
# avoids inheriting sampler history. Subject-local RNG streams make replacement
# deterministic. A non-positive period disables recycling.
.emc_wpool_recycle <- function(pool, i, every, ctx) {
  if (is.null(pool)) return(pool)
  if (is.na(every) || every <= 0L) return(pool)
  # Skip the degenerate first iteration and single-worker pool.
  if (i <= 1L || pool$n <= 1L) return(pool)
  if (((i - 1L) %% every) != 0L) return(pool)

  target <- pool$n
  # A dead pool can recover at a recycle point, but repeated fork failures use
  # backoff rather than retrying every period.
  dead <- !isTRUE(pool$alive)
  if (dead && isTRUE(pool$rebuild_skip > 0L)) {
    pool$rebuild_skip <- pool$rebuild_skip - 1L
    return(pool)
  }

  if (identical(pool$backend, "spawn") && !is.null(pool$template) &&
      isTRUE(pool$template$alive)) {
    .emc_wpool_stop_workers(pool)
    unlink(file.path(pool$dir,
                     c(sprintf("req%d", seq_len(target)),
                       sprintf("ans%d", seq_len(target)))))
    spawned <- .emc_wpool_spawn(pool$dir, seq_len(target), ctx,
                                backend = "spawn", template = pool$template)
    fresh <- if (is.null(spawned)) NULL else {
      pool$jobs <- spawned$jobs
      pool$wcs <- spawned$wcs
      pool$rcs <- spawned$rcs
      pool$alive <- TRUE
      pool$generation <- .emc_wpool_next_generation()
      # The replaced workers may have left records on the channel this pool
      # keeps; retire them before the replacements are given anything to do.
      .emc_wpool_retire(pool)
      pool
    }
  } else {
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(target, ctx)
  }
  if (is.null(fresh)) {
    # Finding 3.  The workers are gone, and the pool is dead until a later
    # recycle can rebuild it -- but a dead pool is still an owner.  This used to
    # return a bare stub with `dir = NULL` and no template, which discarded the
    # only handle on a live R process: `.emc_wpool_stop()` would then find
    # nothing to stop, and the template and its directory outlived the fit.  It
    # also removed the clean-backend recovery route, because a recycle needs the
    # template it had just thrown away.
    #
    # So keep everything and change only what is actually different: no workers,
    # not alive, and a new generation, so that a record written by the workers
    # that just died is recognisable on the completion channel this pool still
    # owns.
    if (!dead) .emc_wpool_degraded("could not re-fork the pool when recycling")
    fails <- if (is.null(pool$rebuild_fails)) 1L else pool$rebuild_fails + 1L
    pool$n <- target
    pool$jobs <- list()
    pool$wcs <- list()
    pool$rcs <- list()
    pool$alive <- FALSE
    pool$generation <- .emc_wpool_next_generation()
    .emc_wpool_retire(pool)
    pool$rebuild_fails <- fails
    pool$rebuild_skip <- min(fails, 8L)
    return(pool)
  }
  fresh
}

.emc_wpool_stop_workers <- function(pool) {
  if (is.null(pool)) return(invisible(NULL))
  # Send shutdown to every worker still there; surviving workers may still await
  # its token.  Writing into a pipe whose worker has already gone earns a
  # SIGPIPE, and R prints that to stderr from a signal handler where no
  # tryCatch() can reach it, so check first rather than explain the noise
  # afterwards.
  #
  # The token is the graceful half and it does work -- a worker that receives it
  # leaves its serve loop and stops answering.  What it does not do, and what
  # closing the pipes below does not do either, is end the process: R opens a
  # blocking read FIFO for reading and writing, so there is no end-of-stream to
  # see.  Ending the process is the backend-specific escalation below.
  for (w in seq_len(length(pool$wcs))) {
    pid <- if (w <= length(pool$jobs)) {
      .emc_wpool_job_pid(pool$jobs[[w]])
    } else NA_integer_
    if (is.na(pid) || .emc_wpool_pid_alive(pid)) {
      .emc_wpool_send(pool$wcs[[w]], NULL)
    }
  }
  for (cn in c(pool$wcs, pool$rcs)) try(close(cn), silent = TRUE)
  if (identical(pool$backend, "spawn")) {
    # What the token and the closed pipes actually achieve, measured rather
    # than assumed: the worker leaves its serve loop and stops answering, but
    # its *process* does not end.  A `mcparallel` child persists until it is
    # collected, and a spawned worker is the template's child, not ours -- only
    # the template can collect it.  Waiting for it here waits for something
    # this process cannot cause, and costs the full deadline every teardown.
    #
    # So: when the template is alive, leave it to do its job.  When it is not
    # -- a failed recycle, or a template already stopped -- nobody will, and a
    # worker with no reaper is a retained R process holding a whole context.
    # Signal it, and do not wait: `waitpid` is not ours to call, and once it is
    # killed and unparented init reaps it.  These are pids this pool holds
    # handles for; nothing here searches for a process by name.
    tpid <- if (is.null(pool$template)) NA_integer_ else {
      .emc_wpool_job_pid(pool$template$job)
    }
    orphaned <- !isTRUE(pool$template$alive) || !.emc_wpool_pid_alive(tpid)
    if (orphaned) {
      .emc_wpool_terminate_jobs(pool$jobs, wait = FALSE, terminate = TRUE)
    }
  } else {
    .emc_wpool_terminate_jobs(pool$jobs, wait = TRUE, terminate = TRUE)
  }
  invisible(NULL)
}

.emc_wpool_template_stop <- function(template) {
  if (is.null(template)) return(invisible(NULL))
  if (isTRUE(template$alive)) {
    .emc_wpool_send(template$wc, NULL)
  }
  for (cn in list(template$wc, template$rc)) try(close(cn), silent = TRUE)
  .emc_wpool_terminate_jobs(list(template$job), wait = TRUE, terminate = FALSE)
  invisible(NULL)
}

.emc_wpool_stop <- function(pool) {
  if (is.null(pool)) return(invisible(NULL))
  .emc_wpool_stop_workers(pool)
  # Recycling retains this channel while workers are replaced.
  if (!is.null(pool$done)) {
    for (cn in pool$done[c("rc", "wc")]) try(close(cn), silent = TRUE)
  }
  if (identical(pool$backend, "spawn")) {
    .emc_wpool_template_stop(pool$template)
  }
  unlink(pool$dir, recursive = TRUE)
  invisible(NULL)
}

# Fork workers may dirty block-context pages, so the default recycle period is
# conservative. The clean template avoids sampler-history sharing; retain an
# occasional reset for context pages without replacing workers every few steps.
.emc_wpool_recycle_default <- function(pool) {
  if (!is.null(pool) && identical(pool$backend, "spawn")) 50L else 10L
}

# --- worker -----------------------------------------------------------------

.emc_wpool_done_connect <- function(req) {
  tryCatch(fifo(file.path(dirname(req), "done"), "wb", blocking = TRUE),
           error = function(e) NULL)
}

# TRUE also when there is no channel: a worker with no way to announce itself
# still has a reply to write, and the master falls back to reading pipes in
# order.  FALSE means the channel exists and the write failed, which is a dead
# channel and the end of this worker's usefulness.
.emc_wpool_done_write <- function(con, msg, status) {
  if (is.null(con)) return(TRUE)
  tryCatch({
    writeBin(.emc_wpool_done_encode(msg$w, msg$generation, msg$request, status),
             con)
    flush(con)
    TRUE
  }, error = function(e) FALSE)
}

.emc_wpool_serve <- function(req, ans, ctx, inherited = NULL) {
  for (cn in inherited) try(close(cn), silent = TRUE)
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  done <- NULL
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE)
            if (!is.null(done)) try(close(done), silent = TRUE) }, add = TRUE)
  repeat {
    msg <- .emc_wpool_recv(rc)
    if (is.null(msg)) break                     # shutdown or master disconnect
    # Announce the pick-up before doing any of the work.  "Never collected the
    # request" and "collected it and is wedged inside the likelihood" are
    # different failures with different answers, and a live pid cannot tell
    # them apart.  Open the channel lazily so ordinary pool messages are
    # unchanged.
    if (isTRUE(msg$notify)) {
      if (is.null(done)) done <- .emc_wpool_done_connect(req)
      if (!.emc_wpool_done_write(done, msg, .EMC_WPOOL_ST_START)) break
    }
    out <- tryCatch(.emc_wpool_compute(msg, ctx),
                    error = function(e) list(failed = conditionMessage(e)))
    # Announce completion before writing the reply: a large reply can block on
    # the pipe, while the master waits for this record to know which reply to
    # read.  The status travels with it, so the master knows a failure is
    # coming before it commits to reading one.
    if (isTRUE(msg$notify)) {
      if (is.null(done)) done <- .emc_wpool_done_connect(req)
      if (!.emc_wpool_done_write(done, msg, if (is.null(out$failed)) {
        .EMC_WPOOL_ST_OK
      } else .EMC_WPOOL_ST_FAILED)) break
    }
    if (!.emc_wpool_send(wc, out)) break
    # One-shot wrapper jobs (currently the clean loo process) compute exactly
    # one request, return it, and exit so their allocator/native state cannot
    # accumulate work from later requests.
    if (isTRUE(ctx$oneshot)) break
  }
  invisible(NULL)
}

# Template entry point: load context once, then fork children from a process
# that has never held sampler history.
.emc_wpool_template <- function(boot) {
  ctx <- readRDS(boot$ctx)
  rc <- fifo(boot$req, "rb", blocking = TRUE)
  wc <- fifo(boot$ans, "wb", blocking = TRUE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)
  .emc_wpool_send(wc, list(ready = TRUE, pid = Sys.getpid()))
  jobs <- list()
  repeat {
    command <- .emc_wpool_recv(rc)
    if (is.null(command)) break
    reply <- tryCatch({
      if (!identical(command$command, "spawn")) stop("unknown template command")
      keys <- as.character(command$idx)
      old <- jobs[keys]
      old <- old[!vapply(old, is.null, logical(1))]
      if (length(old)) parallel::mccollect(old, wait = TRUE)
      jobs[keys] <- NULL
      # Finding 4.  Each child is tracked the moment it exists, not after the
      # loop: a fork that fails on the fourth of eight used to leave the first
      # three outside `jobs`, so the template's own exit handler -- the only
      # thing that reaps them -- did not know they were its.  Cleanup has to own
      # a child at the moment it is acquired, and a partial spawn has to roll
      # back rather than leave half a pool running.
      made <- tryCatch({
        for (w in seq_along(keys)) {
          jobs[[keys[w]]] <- local({
            i <- w
            parallel::mcparallel(
              .emc_wpool_serve(command$req[i], command$ans[i], ctx,
                               inherited = list(rc, wc)),
              detached = FALSE
            )
          })
        }
        jobs[keys]
      }, error = function(e) e)
      if (inherits(made, "error")) {
        born <- jobs[keys]
        born <- born[!vapply(born, is.null, logical(1))]
        if (length(born)) {
          for (job in born) try(tools::pskill(job$pid), silent = TRUE)
          try(parallel::mccollect(born, wait = TRUE), silent = TRUE)
        }
        jobs[keys] <- NULL
        stop(conditionMessage(made))
      }
      list(ok = TRUE, pids = vapply(made, function(x) x$pid, integer(1)))
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    if (!.emc_wpool_send(wc, reply)) break
  }
  if (length(jobs)) {
    for (job in jobs) try(tools::pskill(job$pid), silent = TRUE)
    try(parallel::mccollect(jobs, wait = TRUE), silent = TRUE)
  }
  invisible(NULL)
}

# Run `expr` while restoring the caller's RNG state. Workers and fallback paths
# may install subject streams in the master process.
.emc_with_preserved_rng <- function(expr) {
  if (exists(".Random.seed", envir = globalenv())) {
    seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  expr
}

# Shared unit of work for workers and the serial fallback, ensuring identical
# results after a pool failure.
.emc_wpool_compute <- function(msg, ctx, shared = NULL) {
  # "ll" computes likelihoods for proposal rows and performs no random draws.
  # Particle messages carry group state and use the RNG-preserving path below.
  if (identical(msg$kind, "ll")) {
    return(list(ll = calc_ll_manager(msg$proposals, dadm = ctx$data[[msg$s]],
                                     model = ctx$model, component = msg$component,
                                     r_cores = 1L)))
  }
  # loo::loo parallelises over columns with mclapply. Run it in a clean,
  # one-shot process so those inner workers inherit only this model's matrix,
  # never compare()'s accumulated matrices for earlier models.
  if (identical(msg$kind, "loo")) {
    loo_warnings <- character()
    value <- withCallingHandlers(
      loo::loo(ctx$ll_mat, cores = ctx$cores),
      warning = function(w) {
        loo_warnings <<- c(loo_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    return(list(value = value, warnings = loo_warnings))
  }
  .emc_wpool_compute_particle(msg, ctx, shared)
}

# Execute loo::loo in a clean one-shot process. Its own implementation forks
# once per pointwise column, so isolating that call is the important part: the
# inner forks see only the current model's matrix and disappear with the
# wrapper. A NULL return asks the caller to use a serial loo fallback.
.emc_wpool_run_loo <- function(ll_mat, cores) {
  cores <- suppressWarnings(as.integer(cores))
  if (is.na(cores) || cores <= 1L) return(NULL)
  ctx <- list(kind = "loo", ll_mat = ll_mat, cores = cores, oneshot = TRUE)
  pool <- .emc_wpool_start(1L, ctx, allow_one = TRUE)
  if (is.null(pool)) return(NULL)
  # A fork fallback would still inherit compare()'s accumulated matrices.
  # Only use the history-free spawned backend here; the caller will run loo
  # serially if a clean process cannot be created.
  if (!identical(pool$backend, "spawn")) {
    .emc_wpool_stop(pool)
    return(NULL)
  }
  on.exit(.emc_wpool_stop(pool), add = TRUE)

  sent <- .emc_wpool_request(pool, 1L, list(kind = "loo", notify = FALSE, w = 1L),
                             deadline = .emc_wpool_deadline(
                               .EMC_WPOOL_TIMEOUT_LOO, adapt = FALSE))
  if (!sent) return(NULL)
  # loo's own work is legitimately expensive, so this one gets a generous
  # deadline rather than the sampler's: a bound exists to stop a hang, not to
  # decide how long real work is allowed to take.
  res <- .emc_wpool_reply(pool, 1L,
                          deadline = .emc_wpool_deadline(.EMC_WPOOL_TIMEOUT_LOO,
                                                         adapt = FALSE))
  if (is.null(res) || !is.null(res$failed) || is.null(res$value)) return(NULL)
  if (length(res$warnings) && exists(".loo_warn_store")) {
    .loo_warn_store$msgs <- c(.loo_warn_store$msgs, res$warnings)
  }
  res$value
}

.emc_wpool_compute_particle <- function(msg, ctx, shared = NULL) .emc_with_preserved_rng({
  # Workers receive a path to the shared group draw; the master fallback already
  # has the decoded object.  Keep the raw-byte form as a compatibility path for
  # direct callers and older serialized messages.
  if (is.null(shared)) {
    shared <- if (!is.null(msg$shared_file)) {
      .emc_wpool_shared_read(msg$shared_file)
    } else {
      unserialize(msg$shared)
    }
  }
  # Only the covariance crosses the pipe; its factorisation is rebuilt here.
  # `.chol_factor()` returns the root *and* its inverse, so the cache is three
  # P x P matrices where the covariance alone is one, and the wire is the
  # scarce resource: a blocking FIFO moves ~16 MB/s, while chol() runs at
  # GFLOP/s.  The trade only improves with P (wire is O(P^2), chol O(P^3)/3,
  # and they do not cross until P is in the thousands).  Every worker starts
  # from the identical matrix and makes the identical LAPACK call, so the
  # factors are bit-for-bit the ones the master would have sent.
  group_chol <- shared$group_chol
  if (is.null(group_chol) && !is.null(shared$group_var)) {
    group_chol <- build_group_chol_cache(shared$group_var, shared$idx_list)
  }
  subs <- msg$subs
  props <- matrix(0, ctx$n_pars + 1L, length(subs))
  pm <- vector("list", length(subs))
  seeds <- vector("list", length(subs))
  times <- numeric(length(subs))
  # The master cannot see how long a worker's own request took, only how long
  # it waited for the reply, and those differ by the queue delay plus whatever
  # else the worker was descheduled for.  Report both ends from here.
  t_recv <- Sys.time()
  w_started <- proc.time()
  rejects_before <- .emc_reject_counts("particle")
  for (k in seq_along(subs)) {
    s <- subs[k]
    assign(".Random.seed", msg$seeds[[k]], envir = globalenv())
    # Use CPU time for partition costs so descheduling does not distort routing.
    t0 <- sum(proc.time()[c("user.self", "sys.self")])
    out <- safe_new_particle(
      s = s, data = ctx$data[[s]], pm_settings = msg$pm[[k]],
      eff_mu = ctx$eff_mu[[s]], eff_var = ctx$eff_var[[s]],
      chains_mu = ctx$chains_mu[[s]], chains_var = ctx$chains_var[[s]],
      prev_ll = msg$prev_ll[k], parameters = NULL, model = ctx$model,
      stage = ctx$stage, type = ctx$type, tune = ctx$tune,
      marginalise = ctx$marginalise, r_cores = ctx$r_cores,
      chol_cache = ctx$chol_caches[[s]], group_chol = group_chol,
      current_alpha = msg$alpha[, k],
      population_mu = msg$population_mu[, k],
      # The same matrix new_particle() will compare the cache's `ref` against,
      # so one copy serves both and the cache always hits.
      population_var = shared$group_var
    )
    times[k] <- sum(proc.time()[c("user.self", "sys.self")]) - t0
    props[, k] <- c(out$proposal, out$ll)
    pm[[k]] <- out$pm_settings
    seeds[[k]] <- get(".Random.seed", envir = globalenv())
  }
  w_spent <- proc.time() - w_started
  list(props = props, pm = pm, seeds = seeds, times = times,
       # `times` is per subject and excludes this function's own overhead;
       # `cpu`/`elapsed` cover the whole request, so the difference between
       # them is the worker's own housekeeping rather than the master's.
       t_recv = t_recv, t_done = Sys.time(),
       cpu = unname(w_spent[["user.self"]] + w_spent[["sys.self"]]),
       elapsed = unname(w_spent[["elapsed"]]),
       rejects = .emc_reject_delta(rejects_before, "particle"))
})

# With one subject, distribute proposal rows across persistent workers. Likelihood
# evaluation draws nothing, so row assignment and pool width cannot affect values.
#
# Return NULL when the pool cannot serve the call so the caller computes the
# vector itself.
.emc_wpool_ll <- function(proposals, s, component = NULL, dynamic = FALSE) {
  pool <- .emc_pool_state$ll_pool
  if (is.null(pool) || !isTRUE(pool$alive) || pool$n <= 1L) return(NULL)
  if (!is.matrix(proposals)) return(NULL)
  n <- nrow(proposals)
  # Avoid a round trip when there are too few rows to keep workers busy.
  if (n <= pool$n) return(NULL)
  # Record whether static or dynamic splitting actually ran; routing decisions
  # must distinguish a declined dynamic request from a completed queue.
  # No worker-count cap any more: the completion record carries a 32-bit worker
  # id, so the queue is bounded by what the machine will give it rather than by
  # what fits in a byte.
  .emc_pool_state$ll_dynamic <- isTRUE(dynamic) && !is.null(pool$done) &&
    n >= .EMC_LL_DYN_MIN_ROWS * pool$n
  if (.emc_pool_state$ll_dynamic) {
    return(.emc_wpool_ll_dynamic(pool, proposals, s, component, n))
  }
  idx <- .split_work_indices(n, pool$n)
  n_workers <- max(idx)

  sent <- logical(n_workers)
  for (w in seq_len(n_workers)) {
    sent[w] <- .emc_wpool_request(
      pool, w, list(kind = "ll", s = s, component = component,
                    proposals = proposals[idx == w, , drop = FALSE]))
    if (!sent[w]) {
      .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
      break
    }
  }

  # Drain every sent worker after a failed round so stale replies cannot answer
  # the next call.
  out <- numeric(n)
  ok <- all(sent)
  deadline <- .emc_wpool_deadline()
  for (w in which(sent)) {
    res <- .emc_wpool_reply(pool, w, deadline)
    if (is.null(res) || !is.null(res$failed) ||
        length(res$ll) != sum(idx == w)) {
      ok <- FALSE
      # A missing reply means the worker or transport is gone; a reported
      # failure means the worker remains available.
      if (is.null(res)) .emc_pool_state$ll_pool$alive <- FALSE
      next
    }
    out[idx == w] <- res$ll
  }
  if (!ok) {
    if (!all(sent)) .emc_pool_state$ll_pool$alive <- FALSE
    return(NULL)
  }
  out
}

# Dynamic work queue for heterogeneous likelihood costs. Static splitting
# balances row counts; the queue keeps workers fed when row costs are uneven.
#
# Chunk sizes take a fixed fraction of the remaining rows: early chunks amortise
# pipe overhead and later chunks limit tail imbalance. Single rows are the floor.
#
# Routing is measured on the live fit because queue overhead can outweigh its
# benefit for homogeneous or cheap likelihoods. The queue cannot beat the cost
# of the most expensive single row and does not know costs in advance.
.EMC_LL_DYN_MIN_ROWS <- 4L   # rows per worker below which the queue cannot help
# Chunk size is remaining / (GRAIN * n_workers); tune grain for tail balance
# versus message overhead.
.EMC_LL_DYN_GRAIN <- 4

.emc_wpool_ll_dynamic <- function(pool, proposals, s, component, n) {
  k <- pool$n
  out <- numeric(n)
  pending <- vector("list", k)
  next_row <- 1L

  # Larger grain gives smaller chunks and better balance at the cost of messages.
  grain <- max(1, as.numeric(getOption("emc2.ll_queue_grain", .EMC_LL_DYN_GRAIN)))
  take <- function() {
    left <- n - next_row + 1L
    if (left <= 0L) return(NULL)
    size <- min(left, max(1L, as.integer(ceiling(left / (grain * k)))))
    rows <- seq.int(next_row, length.out = size)
    next_row <<- next_row + size
    rows
  }
  generation <- .emc_wpool_one_int(pool$generation)
  # A fresh request id per chunk, not per call: the queue hands the same worker
  # several chunks, and a record naming a chunk that has already been collected
  # is a late answer rather than the current one.
  issued <- integer(k)
  send <- function(w, rows) {
    issued[w] <<- .emc_wpool_next_request()
    ok <- .emc_wpool_request(
      pool, w, list(kind = "ll", s = s, component = component,
                    notify = TRUE, w = w, generation = generation,
                    request = issued[w],
                    proposals = proposals[rows, , drop = FALSE]))
    if (!ok) {
      .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
    }
    ok
  }
  # Poll the shared completion channel so worker death cannot hang the master,
  # and give the wait a deadline: a worker that is alive and wedged used to hold
  # this queue open indefinitely, because liveness was the only thing being
  # checked.  Returns the worker index, or NA for anything that ends the queue.
  outstanding <- function() which(!vapply(pending, is.null, logical(1)))
  deadline <- .emc_wpool_deadline()
  await <- function() {
    repeat {
      rec <- .emc_wpool_await_record(pool, outstanding, deadline)
      if (!identical(rec$code, "record")) return(NA_integer_)
      # Progress rather than completion: the chunk was collected.
      if (identical(rec$status, .EMC_WPOOL_ST_START)) next
      w <- rec$w
      # A record for a chunk this worker has already been credited with is a
      # late answer from an earlier round of the same queue.
      if (is.na(w) || w < 1L || w > k) return(NA_integer_)
      if (!identical(rec$request, issued[w])) next
      return(w)
    }
  }

  fail <- function() {
    # Drain outstanding completion tokens and replies before abandoning the
    # queue, preventing stale answers from contaminating the next call.
    n_out <- sum(!vapply(pending, is.null, logical(1)))
    for (i in seq_len(n_out)) {
      w <- await()
      if (is.na(w)) break
      .emc_wpool_reply(pool, w, deadline)
    }
    .emc_pool_state$ll_pool$alive <- FALSE
    NULL
  }

  for (w in seq_len(k)) {
    rows <- take()
    if (is.null(rows)) break
    if (!send(w, rows)) return(fail())
    pending[[w]] <- rows
  }

  while (any(!vapply(pending, is.null, logical(1)))) {
    w <- await()
    # An unexpected token or a deleted pending slot means channel/reply
    # alignment is lost; no later result is trustworthy.
    if (is.na(w) || w < 1L || w > k || is.null(pending[[w]])) return(fail())
    rows <- pending[[w]]
    res <- .emc_wpool_reply(pool, w, deadline)
    if (is.null(res) || !is.null(res$failed) || length(res$ll) != length(rows)) {
      pending[w] <- list(NULL)
      return(fail())
    }
    out[rows] <- res$ll
    pending[w] <- list(NULL)
    more <- take()
    if (!is.null(more)) {
      if (!send(w, more)) return(fail())
      pending[[w]] <- more
    }
  }
  out
}

# Register the pool for the single-subject likelihood path without changing
# other `calc_ll_manager()` callers.
.emc_ll_pool_set <- function(pool, s, reset = TRUE) {
  .emc_pool_state$ll_pool <- pool
  .emc_pool_state$ll_subject <- s
  if (reset) .emc_ll_route_reset()
  invisible(NULL)
}

.emc_ll_pool_clear <- function() {
  .emc_pool_state$ll_pool <- NULL
  .emc_pool_state$ll_subject <- NULL
  invisible(NULL)
}

# Routing is empirical: compare serial and pooled likelihood cost per particle
# on the live fit, because proposal counts and model costs vary by iteration.
# Keep the cheaper arm; ties favour serial to avoid unnecessary pool overhead.
.EMC_LL_PROBE_N <- 6L      # calls per arm while probing
.EMC_LL_REPROBE <- 200L    # calls before the decision is re-examined
# A clearly expensive serial call can select the pool immediately; otherwise
# probe both arms before deciding.
.EMC_LL_OBVIOUS <- 0.5     # seconds in a single serial call

.emc_ll_route_reset <- function() {
  .emc_pool_state$ll_route <- list(mode = "probe", t_serial = numeric(0),
                                   t_pool = numeric(0), since = 0L,
                                   split = "probe", t_static = numeric(0),
                                   t_dynamic = numeric(0), since_split = 0L)
  invisible(NULL)
}

# Interleave serial and pooled probes so both arms see the same drift in particle
# costs rather than comparing separate clouds.
.emc_ll_route_use_pool <- function() {
  st <- .emc_pool_state$ll_route
  if (is.null(st)) { .emc_ll_route_reset(); st <- .emc_pool_state$ll_route }
  switch(st$mode,
         probe = length(st$t_pool) < length(st$t_serial),
         pool = TRUE,
         FALSE)
}

# Once the pool wins, independently choose static splitting or the dynamic queue.
# Probe this nested decision only while the pool arm is active, then periodically
# remeasure as particle costs and caches change.
.emc_ll_route_use_dynamic <- function() {
  if (!isTRUE(getOption("emc2.ll_queue", TRUE))) return(FALSE)
  st <- .emc_pool_state$ll_route
  if (is.null(st) || !identical(st$mode, "pool")) return(FALSE)
  switch(st$split,
         probe = length(st$t_dynamic) < length(st$t_static),
         dynamic = TRUE,
         FALSE)
}

.emc_ll_route_record <- function(pooled, seconds, n_particles, dynamic = FALSE) {
  st <- .emc_pool_state$ll_route
  if (is.null(st) || !is.finite(seconds) || n_particles <= 0L) return(invisible(NULL))
  per <- seconds / n_particles
  if (pooled) st$t_pool <- c(st$t_pool, per) else st$t_serial <- c(st$t_serial, per)
  if (pooled && identical(st$mode, "pool")) {
    if (dynamic) st$t_dynamic <- c(st$t_dynamic, per)
    else st$t_static <- c(st$t_static, per)
    if (identical(st$split, "probe")) {
      if (length(st$t_dynamic) >= .EMC_LL_PROBE_N &&
          length(st$t_static) >= .EMC_LL_PROBE_N) {
        # Ties go to static, which costs no polling and no extra messages.
        st$split <- if (stats::median(st$t_dynamic) < stats::median(st$t_static))
          "dynamic" else "static"
        st$since_split <- 0L
      }
    } else {
      st$since_split <- st$since_split + 1L
      if (st$since_split >= .EMC_LL_REPROBE) {
        st$split <- "probe"; st$t_static <- numeric(0)
        st$t_dynamic <- numeric(0); st$since_split <- 0L
      }
    }
  }
  if (st$mode == "probe") {
    if (!pooled && seconds >= .EMC_LL_OBVIOUS) {
      st$mode <- "pool"; st$since <- 0L
    } else if (length(st$t_pool) >= .EMC_LL_PROBE_N &&
               length(st$t_serial) >= .EMC_LL_PROBE_N) {
      # Median, not mean: one straggler call must not decide the arm.  A tie
      # goes to serial -- the pool costs processes and memory, so it has to
      # actually win something to be worth keeping in the path.
      st$mode <- if (stats::median(st$t_pool) < stats::median(st$t_serial))
        "pool" else "serial"
      st$since <- 0L
    }
  } else {
    st$since <- st$since + 1L
    if (st$since >= .EMC_LL_REPROBE) {
      # Re-time both arms as particle costs change; retain the split state so an
      # outer probe does not discard the nested routing decision.
      st <- list(mode = "probe", t_serial = numeric(0), t_pool = numeric(0),
                 since = 0L, split = st$split, t_static = st$t_static,
                 t_dynamic = st$t_dynamic, since_split = st$since_split)
    }
  }
  .emc_pool_state$ll_route <- st
  invisible(NULL)
}

# --- one iteration ----------------------------------------------------------

# Returns the same `proposals` matrix the mcmapply path builds, plus the
# updated per-subject settings, streams and timings.
.emc_wpool_iter <- function(pool, ctx, part, pars, group_chol, pm_settings,
                            prev_ll, seeds) {
  profile <- isTRUE(getOption("emc2.sampler_profile", FALSE))
  iter_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  n_subjects <- length(pm_settings)
  props <- matrix(0, ctx$n_pars + 1L, n_subjects)
  times <- numeric(n_subjects)

  # Only the population covariance and its factors are common to every worker.
  # The current random effects and subject-specific population means are sliced
  # by assignment, so each subject slice crosses only its assigned request pipe.
  # This retains per-iteration LPT balancing without its former O(workers * N)
  # subject-state broadcast.  Gibbs-only fields (tvinv, a_half, factor state,
  # etc.) never enter a particle message.
  population_mu <- vapply(seq_len(n_subjects), function(s) {
    as.numeric(get_group_level(pars, s, ctx$type)$mu)
  }, numeric(ctx$n_pars))
  rownames(population_mu) <- rownames(pars$alpha)
  alpha <- pars$alpha
  shared_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  # Publish the covariance, not its factorisation.  `build_group_chol_cache()`
  # holds the root and its inverse alongside the matrix itself, so the full
  # cache is three P x P matrices per worker per iteration; the covariance is
  # one.  The bytes are written once below and each worker receives only the
  # path; it rebuilds the rest in microseconds.  See
  # .emc_wpool_compute_particle().  The master keeps the cache it already built
  # so the recompute fallback below pays nothing for this.
  shared <- list(group_var = group_chol$ref, idx_list = group_chol$idx_list,
                 group_chol = group_chol)
  # A dead pool computes every share in the master from `shared` itself, so the
  # wire copy would be built and thrown away on every iteration of a fit that
  # has no workers at all.
  was_alive <- isTRUE(pool$alive)
  shared_raw <- if (was_alive) {
    serialize(shared[c("group_var", "idx_list")], NULL)
  } else raw(0)
  shared_file <- NULL
  # The file is removed only after every response or serial recomputation has
  # completed.  This also covers interrupts and errors while a worker is being
  # retired, so a surviving worker cannot observe a missing payload halfway
  # through the iteration.
  on.exit(if (!is.null(shared_file)) unlink(shared_file), add = TRUE)
  if (was_alive) {
    shared_file <- .emc_wpool_shared_write(pool$dir, shared_raw)
    if (is.null(shared_file)) {
      .emc_wpool_degraded(.emc_wpool_transport_error(
        "could not publish the shared broadcast"))
      pool$alive <- FALSE
    }
  }
  shared_elapsed <- if (profile) {
    proc.time()[["elapsed"]] - shared_started
  } else NA_real_
  # `notify`/`w` ask the worker to announce completion on the shared,
  # non-blocking `done` channel before it writes its (possibly large) reply.
  # That is what lets the master wait with a deadline instead of blocking on one
  # worker's pipe forever; see .emc_wpool_await_done().
  watch <- !is.null(pool$done)
  generation <- .emc_wpool_one_int(pool$generation)
  # One request id for the whole iteration: every worker is answering the same
  # question, and a record naming a different one is answering a question that
  # has already been settled.
  request <- .emc_wpool_next_request()
  msgs <- lapply(seq_along(part), function(w) {
    subs <- part[[w]]
    list(subs = subs, shared_file = shared_file,
         alpha = alpha[, subs, drop = FALSE],
         population_mu = population_mu[, subs, drop = FALSE],
         pm = pm_settings[subs], prev_ll = prev_ll[subs], seeds = seeds[subs],
         notify = watch, w = w, generation = generation, request = request)
  })
  sent <- rep(FALSE, pool$n)
  # Per-worker telemetry.  `dispatch_at`/`finish_at` are master elapsed seconds
  # relative to the start of this iteration; `dispatch_wall` is wall clock,
  # because the queue delay is the gap between two *different* processes and
  # `proc.time()` is measured from each process's own start.
  dispatch_at <- finish_at <- rep(NA_real_, pool$n)
  dispatch_wall <- rep(Sys.time()[NA], pool$n)
  worker_cpu <- worker_elapsed <- queue_delay <- rep(NA_real_, pool$n)
  fallback_workers <- 0L
  fallback_subjects <- 0L
  rejects <- stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                             .EMC_REJECT_CLASSES)
  send_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  sendable <- was_alive && !is.null(shared_file)
  if (sendable) {
    for (w in seq_len(pool$n)) {
      if (!length(part[[w]])) next
      sent[w] <- .emc_wpool_request(pool, w, msgs[[w]])
      if (profile && sent[w]) {
        dispatch_at[w] <- proc.time()[["elapsed"]] - iter_started
        dispatch_wall[w] <- Sys.time()
      }
      if (!sent[w]) {
        .emc_wpool_degraded(
          .emc_wpool_transport_error("send to a worker failed"))
        pool$alive <- FALSE
      }
    }
  }
  send_elapsed <- if (profile) {
    proc.time()[["elapsed"]] - send_started
  } else NA_real_
  receive_started <- proc.time()[["elapsed"]]
  # Sends are issued to every worker before any receive, and a send blocks until
  # that worker drains its pipe, so worker n starts a little after workers
  # 1..n-1 have each taken their message.  At current message sizes that skew is
  # far below the LPT partition's granularity; noted here so it is not
  # rediscovered later as a mystery imbalance.
  store <- function(w, res, source = "worker") {
    subs <- part[[w]]
    props[, subs] <<- res$props
    times[subs] <<- res$times
    pm_settings[subs] <<- res$pm
    seeds[subs] <<- res$seeds
    if (!is.null(res$rejects)) rejects <<- rejects + res$rejects
    # A pool that was already dead on entry (a single worker, or Windows) runs
    # every share in the master by design; that is the serial route, not a
    # fallback, and counting it would report a healthy fit as failing.
    if (identical(source, "master") && was_alive) {
      fallback_workers <<- fallback_workers + 1L
      fallback_subjects <<- fallback_subjects + length(subs)
    }
    if (!profile) return(invisible(NULL))
    finish_at[w] <<- proc.time()[["elapsed"]] - iter_started
    if (!is.null(res$cpu)) worker_cpu[w] <<- res$cpu
    if (!is.null(res$elapsed)) worker_elapsed[w] <<- res$elapsed
    # Only a genuine worker reply measures a queue: a share recomputed in the
    # master was never queued, and recording its "delay" would report the
    # failure-handling path as though the pool had been slow to pick work up.
    if (identical(source, "worker") && !is.null(res$t_recv) &&
        !is.na(dispatch_wall[w])) {
      queue_delay[w] <<- as.numeric(difftime(res$t_recv, dispatch_wall[w],
                                             units = "secs"))
    }
    invisible(NULL)
  }
  # A `failed` reply is not a dead worker: it caught an error in the work and is
  # still listening.  Recompute that share here -- same function, same streams,
  # so the result is identical -- but do not condemn the whole block over one
  # subject's bad iteration; if the error is deterministic the master's own call
  # raises it properly.
  recompute <- function(w) .emc_wpool_compute(msgs[[w]], ctx, shared = shared)

  pending <- sent & vapply(part, length, integer(1)) > 0L
  # Which requests were collected by a worker, as opposed to merely sent to one
  # that is still alive.  A worker that never picked its request up and a worker
  # wedged inside the likelihood look identical from a pid.
  picked_up <- rep(FALSE, pool$n)
  if (watch && any(pending)) {
    deadline <- .emc_wpool_deadline()
    while (any(pending)) {
      rec <- .emc_wpool_await_record(pool, function() which(pending), deadline)
      w <- if (identical(rec$code, "record")) rec$w else NA_integer_
      # A record for a request that has already been settled belongs to an
      # earlier iteration of this same pool -- a worker that answered after its
      # share had been recomputed here.  Acting on it would credit this
      # iteration with last iteration's work.
      if (identical(rec$code, "record") &&
          !identical(rec$request, request)) {
        .emc_pool_state$stale_records <-
          (if (is.null(.emc_pool_state$stale_records)) 0L
           else .emc_pool_state$stale_records) + 1L
        next
      }
      # Progress, not completion: the worker has the request and is working.
      if (identical(rec$code, "record") &&
          identical(rec$status, .EMC_WPOOL_ST_START) &&
          !is.na(w) && w >= 1L && w <= pool$n) {
        picked_up[w] <- TRUE
        next
      }
      if (is.na(w) || w < 1L || w > pool$n || !pending[w]) {
        if (identical(rec$code, "timeout")) {
          # Every worker we are still waiting on is alive and has produced
          # nothing for `deadline` seconds: one of them is wedged.  Name it and
          # the subjects it holds -- that is the diagnosis a frozen progress bar
          # cannot give -- then kill it so the recomputation below raises the
          # underlying error properly.  Saying whether it ever collected the
          # request separates a wedged likelihood from a worker that never
          # started, which need different fixes.
          stuck <- which(pending)
          idle <- stuck[!picked_up[stuck]]
          .emc_wpool_degraded(sprintf(
            "worker%s %s produced no reply within %.0f s (subject%s %s)%s",
            if (length(stuck) > 1L) "s" else "", paste(stuck, collapse = ", "),
            deadline, if (length(unlist(part[stuck])) > 1L) "s" else "",
            paste(unlist(part[stuck]), collapse = ", "),
            if (length(idle)) {
              sprintf("; worker%s %s never collected the request",
                      if (length(idle) > 1L) "s" else "",
                      paste(idle, collapse = ", "))
            } else ""))
          # Kill *and reap*: a forked worker that is signalled and never
          # collected is a zombie for the rest of the block.  Ownership is the
          # point -- these are jobs this pool holds handles for, not pids found
          # by searching for a process name.
          .emc_wpool_terminate_jobs(pool$jobs[stuck], wait = TRUE,
                                    terminate = TRUE)
        } else if (identical(rec$code, "dead") || is.na(w)) {
          .emc_wpool_degraded("worker gave no reply")
        } else {
          .emc_wpool_degraded("lost alignment on the worker completion channel")
        }
        # Nothing further on these pipes can be trusted or waited on.  Retire
        # the generation before recomputing rather than after: a record that
        # arrives while the master is redoing this work belongs to an
        # assignment that has already been abandoned, and leaving it in the
        # channel only defers the confusion to whoever reads next.
        pool$alive <- FALSE
        .emc_wpool_retire(pool)
        for (k in which(pending)) store(k, recompute(k), source = "master")
        pending[] <- FALSE
        break
      }
      # The record said a reply is on its way, so this read is owed bytes.  It
      # still gets a deadline: the announcement and the payload are separate
      # I/O, and a worker can die between them.
      res <- .emc_wpool_reply(pool, w, deadline)
      source <- "worker"
      if (is.null(res) || !is.null(res$failed)) {
        if (is.null(res)) {
          .emc_wpool_degraded(.emc_wpool_transport_error("worker gave no reply"))
          pool$alive <- FALSE
        }
        res <- recompute(w)
        source <- "master"
      }
      store(w, res, source = source)
      pending[w] <- FALSE
    }
    # Shares that were never sent (a failed send already retired the pool).
    for (w in seq_len(pool$n)) {
      if (length(part[[w]]) && !sent[w]) store(w, recompute(w), source = "master")
    }
  } else {
    # No completion channel: replies are read in order.  They are still bounded,
    # which is the difference between a slow block and a hung chain.
    deadline <- .emc_wpool_deadline()
    for (w in seq_len(pool$n)) {
      subs <- part[[w]]
      if (!length(subs)) next
      res <- if (sent[w]) .emc_wpool_reply(pool, w, deadline) else NULL
      source <- "worker"
      if (is.null(res) || !is.null(res$failed)) {
        if (sent[w] && is.null(res)) {
          # No reply means the transport or the worker itself is gone, so there
          # is nothing left to send to and the rest of the block is serial.
          .emc_wpool_degraded("worker gave no reply")
          pool$alive <- FALSE
        }
        res <- recompute(w)
        source <- "master"
      }
      store(w, res, source = source)
    }
  }
  receive_elapsed <- proc.time()[["elapsed"]] - receive_started
  if (was_alive) .emc_wpool_record_receive(receive_elapsed)
  # Sizing a message means serialising it a second time purely to measure it,
  # so it stays behind the expensive switch rather than being paid for by every
  # profiled run.
  measure_wire <- .emc_profile_expensive()
  private_bytes <- if (measure_wire) {
    sum(vapply(msgs, function(msg) {
      msg$shared <- NULL
      msg$shared_file <- NULL
      length(serialize(msg, NULL))
    }, integer(1)))
  } else NA_real_
  shared_path_bytes <- if (measure_wire && !is.null(shared_file)) {
    active <- vapply(part, length, integer(1)) > 0L
    sum(vapply(msgs[active], function(msg) {
      length(serialize(list(shared_file = msg$shared_file), NULL))
    }, integer(1)))
  } else if (measure_wire) 0 else NA_real_
  # The workload actually assigned to each worker, which is what a scheduler
  # change has to move.  `max(times)` is the slowest single *subject*: with
  # several subjects per worker the two differ by roughly that ratio, and using
  # the per-subject figure as "the slowest worker" made every subject after the
  # first in a partition look like transport overhead.
  worker_load <- vapply(part, function(subs) sum(times[subs]), numeric(1))
  # `rejects` is outside the profile block on purpose: a failure that was not
  # numerical rejection has to be visible in an ordinary run, not only in a
  # profiled one, and it is four integers.
  list(props = props, pm_settings = pm_settings, seeds = seeds, times = times,
       alive = pool$alive, rejects = rejects,
       profile = if (profile) list(
         elapsed = proc.time()[["elapsed"]] - iter_started,
         shared_serialize = shared_elapsed,
         send = send_elapsed,
         receive = receive_elapsed,
         subject_cpu_max = if (length(times)) max(times) else 0,
         subject_cpu_sum = sum(times),
         worker_cpu_max = if (length(worker_load)) max(worker_load) else 0,
         worker_cpu_sum = sum(worker_cpu, na.rm = TRUE),
         worker_elapsed_max = if (all(is.na(worker_elapsed))) NA_real_ else {
           max(worker_elapsed, na.rm = TRUE)
         },
         worker_elapsed_sum = sum(worker_elapsed, na.rm = TRUE),
         workers_active = sum(vapply(part, length, integer(1)) > 0L),
         dispatch_first = .emc_range_or_na(dispatch_at, min),
         dispatch_last = .emc_range_or_na(dispatch_at, max),
         finish_first = .emc_range_or_na(finish_at, min),
         finish_last = .emc_range_or_na(finish_at, max),
         queue_delay_max = .emc_range_or_na(queue_delay, max),
         queue_delay_mean = .emc_range_or_na(queue_delay, mean),
         fallback_workers = fallback_workers,
         fallback_subjects = fallback_subjects,
         degraded = as.integer(was_alive && !isTRUE(pool$alive)),
         rejects = rejects,
         shared_bytes = if (measure_wire) length(shared_raw) else NA_real_,
         private_bytes = private_bytes,
         wire_bytes = if (!measure_wire) NA_real_ else if (!is.null(shared_file)) {
           length(shared_raw) + shared_path_bytes + private_bytes
         } else private_bytes
       ) else NULL)
}

# `max(numeric(0))` warns and returns -Inf, and `max(all-NA, na.rm = TRUE)` does
# the same; a profile column should read NA there instead of an infinity that
# then propagates through every mean taken over it.
.emc_range_or_na <- function(x, f) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else as.numeric(f(x))
}

# --- partitioning -----------------------------------------------------------

# Longest-processing-time-first: repeatedly give the next most expensive
# subject to the worker with the least work so far.  Classic 4/3-approximation
# to the optimal makespan, and unlike mc.preschedule's contiguous equal-count
# split it actually looks at cost.
.emc_lpt_partition <- function(cost, n_workers) {
  n <- length(cost)
  if (n_workers <= 1L) return(list(seq_len(n)))
  out <- vector("list", n_workers)
  load <- numeric(n_workers)
  for (s in order(cost, decreasing = TRUE)) {
    w <- which.min(load)
    out[[w]] <- c(out[[w]], s)
    load[w] <- load[w] + cost[s]
  }
  lapply(out, function(x) if (is.null(x)) integer(0) else x)
}

# First-iteration cost proxy: what the likelihood actually loops over.
.emc_subject_cost <- function(data) {
  vapply(data, function(d) {
    if (is.data.frame(d)) nrow(d) else sum(vapply(d, function(x)
      if (is.data.frame(x)) nrow(x) else 0L, numeric(1)))
  }, numeric(1))
}

# One independent stream per subject, so the draws do not depend on how the
# subjects were partitioned or on how many workers there are.
#
# Streams have to come from L'Ecuyer-CMRG, but fitting must not leave a user's
# session switched over to it -- that would quietly change what every later
# `set.seed()` of theirs produces.  Where the caller is already on L'Ecuyer
# (the normal path: `make_emc()` puts it there) the generator is used and
# advanced directly.  Where it is not, a temporary one is seeded from a single
# draw off the caller's own generator, so the streams still follow from the
# caller's seed, and the caller's generator is handed back exactly as found
# apart from that one draw.
.emc_subject_streams <- function(n_subjects) {
  old_kind <- RNGkind()
  if (!identical(old_kind[1], "L'Ecuyer-CMRG")) {
    s <- sample.int(.Machine$integer.max, 1L)
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      RNGkind(old_kind[1], old_kind[2], old_kind[3])
      assign(".Random.seed", old_seed, envir = globalenv())
    }, add = TRUE)
    set.seed(s, kind = "L'Ecuyer-CMRG")
  }
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1)
  seed <- get(".Random.seed", envir = globalenv())
  streams <- vector("list", n_subjects)
  for (s in seq_len(n_subjects)) {
    seed <- parallel::nextRNGStream(seed)
    streams[[s]] <- seed
  }
  assign(".Random.seed", parallel::nextRNGStream(seed), envir = globalenv())
  streams
}
