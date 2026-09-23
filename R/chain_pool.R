
.emc_pool_state <- new.env(parent = emptyenv())

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


.EMC_WPOOL_READ_CHUNK <- 65536L
.EMC_WPOOL_WRITE_CHUNK <- 4096L
.EMC_WPOOL_MAX_FRAME <- 2^31
.EMC_WPOOL_MAX_SEGMENT <- 2^30

.emc_wpool_transport_fail <- function(msg) {
  .emc_pool_state$last_transport_error <- msg
  FALSE
}

.emc_wpool_transport_error <- function(default) {
  msg <- .emc_pool_state$last_transport_error
  if (is.null(msg) || !nzchar(msg)) default else msg
}

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

.emc_wpool_request <- function(pool, w, msg, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_send(pool$wcs[[w]], msg, deadline = deadline)
}

.emc_wpool_reply <- function(pool, w, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_recv(pool$rcs[[w]], deadline = deadline,
                  alive = .emc_wpool_alive_fn(pool, w))
}

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
    if (!is.finite(n) || n < 0 || n > .EMC_WPOOL_MAX_FRAME) {
      .emc_wpool_transport_fail("implausible message length on the pipe")
      return(NULL)
    }
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

.EMC_WPOOL_REAP_WAIT <- 2

.emc_wpool_reap_wait <- function() {
  opt <- suppressWarnings(as.numeric(
    getOption("emc2.worker_reap_wait", .EMC_WPOOL_REAP_WAIT)))
  if (length(opt) != 1L || is.na(opt) || opt < 0) .EMC_WPOOL_REAP_WAIT else opt
}

.EMC_WPOOL_TIMEOUT_FLOOR <- 120
.EMC_WPOOL_TIMEOUT_COLD <- 600
.EMC_WPOOL_TIMEOUT_TEMPLATE <- 300
.EMC_WPOOL_TIMEOUT_LOO <- 24 * 3600

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

.emc_wpool_await_record <- function(pool, pending, deadline = NULL,
                                    on_start = NULL) {
  started <- proc.time()[["elapsed"]]
  idle <- 0L
  wait <- 0
  repeat {
    rec <- .emc_wpool_done_take(pool$done)
    if (!is.null(rec)) {
      if (identical(rec$generation, .emc_wpool_one_int(pool$generation))) {
        if (is.function(on_start) &&
            identical(rec$status, .EMC_WPOOL_ST_START)) {
          on_start(rec)
          next
        }
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
  if (!dir.exists("/proc")) return(TRUE)
  stat <- suppressWarnings(
    tryCatch(readLines(sprintf("/proc/%d/stat", as.integer(pid)), n = 1L),
             error = function(e) character()))
  if (!length(stat)) return(FALSE)
  tail <- sub("^.*\\) ", "", stat)
  state <- strsplit(tail, " ", fixed = TRUE)[[1L]][1L]
  !(state %in% c("Z", "X", "x"))
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
    suppressWarnings(try(parallel::mccollect(forked, wait = wait), silent = TRUE))
  }
  if (length(spawned)) {
    pids <- vapply(spawned, .emc_wpool_job_pid, integer(1))
    if (wait) {
      deadline <- proc.time()[["elapsed"]] + .emc_wpool_reap_wait()
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
      unlink(c(req, ans))
    }
  }, add = TRUE)

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
      f <- NULL
      deadline <- proc.time()[["elapsed"]] + 20
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
      blocking <- fifo(req[w], "wb", blocking = TRUE)
      close(f)
      wcs[[w]] <- blocking
    }
    for (w in seq_len(n)) {
      rcs[[w]] <- fifo(ans[w], "rb", blocking = FALSE)
    }
    success <- TRUE
  }, error = function(e) NULL)

  if (!success) return(NULL)
  list(jobs = jobs, wcs = wcs, rcs = rcs)
}

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

.emc_wpool_start <- function(n_workers, ctx, allow_one = FALSE) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L ||
      (!isTRUE(allow_one) && n_workers <= 1L)) return(NULL)
  if (Sys.info()[["sysname"]] == "Windows") return(NULL)
  if (!nzchar(Sys.which("mkfifo"))) return(NULL)
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

  .emc_pool_state$warned <- FALSE
  spawned <- .emc_wpool_spawn(dir, seq_len(n_workers), ctx, backend,
                              template = template)
  if (is.null(spawned)) {
    if (!is.null(template)) .emc_wpool_template_stop(template)
    unlink(dir, recursive = TRUE)
    return(NULL)
  }

  list(n = n_workers, dir = dir, jobs = spawned$jobs, wcs = spawned$wcs,
       rcs = spawned$rcs, alive = TRUE, backend = backend,
       generation = .emc_wpool_next_generation(),
       ctx_file = ctx_file, lib = lib, template = template,
       done = .emc_wpool_done_open(dir))
}

.EMC_WPOOL_DONE_BYTES <- 16L
.EMC_WPOOL_ST_START <- 1L
.EMC_WPOOL_ST_OK <- 2L
.EMC_WPOOL_ST_FAILED <- 3L

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

.emc_wpool_next_generation <- function() {
  n <- .emc_pool_state$generation_seq
  n <- if (is.null(n)) 1L else n + 1L
  if (n >= .Machine$integer.max) n <- 1L
  .emc_pool_state$generation_seq <- n
  n
}

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

.emc_wpool_done_take <- function(done) {
  buf <- done$buf
  have <- length(buf$bytes)
  if (have < .EMC_WPOOL_DONE_BYTES) {
    chunk <- tryCatch(readBin(done$rc, "raw", .EMC_WPOOL_DONE_BYTES - have),
                      error = function(e) raw(0))
    if (!length(chunk)) return(NULL)
    buf$bytes <- if (have) c(buf$bytes, chunk) else chunk
    if (length(buf$bytes) < .EMC_WPOOL_DONE_BYTES) return(NULL)
  }
  bytes <- buf$bytes
  buf$bytes <- raw(0)
  .emc_wpool_done_decode(bytes)
}

.emc_wpool_done_open <- function(dir) {
  path <- file.path(dir, "done")
  made <- tryCatch(
    system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE),
    warning = function(e) 1L, error = function(e) 1L
  )
  if (!identical(as.integer(made), 0L)) return(NULL)
  rc <- tryCatch(fifo(path, "rb", blocking = FALSE), error = function(e) NULL)
  if (is.null(rc)) return(NULL)
  wc <- tryCatch(fifo(path, "wb", blocking = FALSE), error = function(e) NULL)
  if (is.null(wc)) { try(close(rc), silent = TRUE); return(NULL) }
  buf <- new.env(parent = emptyenv())
  buf$bytes <- raw(0)
  list(rc = rc, wc = wc, buf = buf)
}


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

.EMC_WPOOL_CONN_PER_WORKER <- 2L
.EMC_WPOOL_CONN_OVERHEAD <- 6L

.emc_wpool_conn_limit <- function() {
  limit <- suppressWarnings(as.integer(
    Sys.getenv("R_MAX_NUM_CONNECTIONS", "128")))
  if (is.na(limit) || limit <= 0L) 128L else limit
}

.emc_wpool_feasible_workers <- function(n_workers) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L) return(n_workers)
  open <- tryCatch(length(getAllConnections()),
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

.emc_wpool_grow <- function(pool, n_total, ctx) {
  n_total <- as.integer(n_total)
  if (is.null(pool) || !isTRUE(pool$alive) || is.na(n_total)) return(pool)
  n_total <- min(n_total, pool$n + .emc_wpool_feasible_workers(n_total - pool$n))
  n_new <- n_total - pool$n
  if (n_new <= 0L) return(pool)

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
    pool$grow_skip <- bitwShiftL(1L, min(fails, 6L))
    return(pool)
  }

  .emc_wpool_retire(pool)
  list(n = n_total, dir = pool$dir, jobs = c(pool$jobs, spawned$jobs),
       wcs = c(pool$wcs, spawned$wcs), rcs = c(pool$rcs, spawned$rcs),
       alive = TRUE, grow_fails = 0L, grow_skip = 0L,
       generation = .emc_wpool_next_generation(),
       backend = pool$backend, ctx_file = pool$ctx_file, lib = pool$lib,
       template = pool$template, done = pool$done)
}

.emc_wpool_recycle <- function(pool, i, every, ctx) {
  if (is.null(pool)) return(pool)
  if (is.na(every) || every <= 0L) return(pool)
  if (i <= 1L || pool$n <= 1L) return(pool)
  if (((i - 1L) %% every) != 0L) return(pool)

  target <- pool$n
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
      .emc_wpool_retire(pool)
      pool
    }
  } else {
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(target, ctx)
  }
  if (is.null(fresh)) {
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
  if (!is.null(pool$done)) {
    for (cn in pool$done[c("rc", "wc")]) try(close(cn), silent = TRUE)
  }
  if (identical(pool$backend, "spawn")) {
    .emc_wpool_template_stop(pool$template)
  }
  unlink(pool$dir, recursive = TRUE)
  invisible(NULL)
}

.emc_wpool_recycle_default <- function(pool) {
  if (!is.null(pool) && identical(pool$backend, "spawn")) 50L else 10L
}


.emc_wpool_done_connect <- function(req) {
  tryCatch(fifo(file.path(dirname(req), "done"), "wb", blocking = TRUE),
           error = function(e) NULL)
}

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
  if (identical(tryCatch(emc_arm_parent_death(), error = function(e) 0L), -1L)) {
    return(invisible(NULL))
  }
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  done <- NULL
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE)
            if (!is.null(done)) try(close(done), silent = TRUE) }, add = TRUE)
  repeat {
    msg <- .emc_wpool_recv(rc)
    if (is.null(msg)) break
    if (isTRUE(msg$notify)) {
      if (is.null(done)) done <- .emc_wpool_done_connect(req)
      if (!.emc_wpool_done_write(done, msg, .EMC_WPOOL_ST_START)) break
    }
    out <- tryCatch(.emc_wpool_compute(msg, ctx),
                    error = function(e) list(failed = conditionMessage(e)))
    if (isTRUE(msg$notify)) {
      if (is.null(done)) done <- .emc_wpool_done_connect(req)
      if (!.emc_wpool_done_write(done, msg, if (is.null(out$failed)) {
        .EMC_WPOOL_ST_OK
      } else .EMC_WPOOL_ST_FAILED)) break
    }
    if (!.emc_wpool_send(wc, out)) break
    if (isTRUE(ctx$oneshot)) break
  }
  invisible(NULL)
}

.emc_wpool_template <- function(boot) {
  ctx <- readRDS(boot$ctx)
  rc <- fifo(boot$req, "rb", blocking = TRUE)
  wc <- fifo(boot$ans, "wb", blocking = TRUE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)
  .emc_wpool_send(wc, list(ready = TRUE, pid = Sys.getpid()))
  jobs <- list()
  orphaned <- FALSE
  repeat {
    command <- .emc_wpool_recv(rc)
    if (is.null(command)) { orphaned <- TRUE; break }
    reply <- tryCatch({
      if (!identical(command$command, "spawn")) stop("unknown template command")
      keys <- as.character(command$idx)
      old <- jobs[keys]
      old <- old[!vapply(old, is.null, logical(1))]
      if (length(old)) parallel::mccollect(old, wait = TRUE)
      jobs[keys] <- NULL
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
    if (!.emc_wpool_send(wc, reply)) { orphaned <- TRUE; break }
  }
  if (length(jobs)) {
    for (job in jobs) try(tools::pskill(job$pid), silent = TRUE)
    try(parallel::mccollect(jobs, wait = TRUE), silent = TRUE)
  }
  if (orphaned) {
    dir <- dirname(boot$req)
    if (is.character(dir) && length(dir) == 1L && nzchar(dir)) {
      try(unlink(dir, recursive = TRUE), silent = TRUE)
    }
  }
  invisible(NULL)
}

.emc_with_preserved_rng <- function(expr) {
  if (exists(".Random.seed", envir = globalenv())) {
    seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  expr
}

.emc_wpool_compute <- function(msg, ctx, shared = NULL) {
  if (identical(msg$kind, "ll")) {
    return(list(ll = calc_ll_manager(msg$proposals, dadm = ctx$data[[msg$s]],
                                     model = ctx$model, component = msg$component,
                                     r_cores = 1L, varying = msg$varying)))
  }
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

.emc_wpool_run_loo <- function(ll_mat, cores) {
  cores <- suppressWarnings(as.integer(cores))
  if (is.na(cores) || cores <= 1L) return(NULL)
  ctx <- list(kind = "loo", ll_mat = ll_mat, cores = cores, oneshot = TRUE)
  pool <- .emc_wpool_start(1L, ctx, allow_one = TRUE)
  if (is.null(pool)) return(NULL)
  if (!identical(pool$backend, "spawn")) {
    .emc_wpool_stop(pool)
    return(NULL)
  }
  on.exit(.emc_wpool_stop(pool), add = TRUE)

  sent <- .emc_wpool_request(pool, 1L, list(kind = "loo", notify = FALSE, w = 1L),
                             deadline = .emc_wpool_deadline(
                               .EMC_WPOOL_TIMEOUT_LOO, adapt = FALSE))
  if (!sent) return(NULL)
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
  if (is.null(shared)) {
    shared <- if (!is.null(msg$shared_file)) {
      .emc_wpool_shared_read(msg$shared_file)
    } else {
      unserialize(msg$shared)
    }
  }
  group_chol <- shared$group_chol
  if (is.null(group_chol) && !is.null(shared$group_var)) {
    group_chol <- build_group_chol_cache(shared$group_var, shared$idx_list,
                                         shared$marginal_idx,
                                         include_full = if (is.null(shared$include_full))
                                           TRUE else shared$include_full)
  }
  subs <- msg$subs
  props <- matrix(0, ctx$n_pars + 1L, length(subs))
  pm <- vector("list", length(subs))
  seeds <- vector("list", length(subs))
  times <- numeric(length(subs))
  t_recv <- Sys.time()
  w_started <- proc.time()
  rejects_before <- .emc_reject_counts("particle")
  for (k in seq_along(subs)) {
    s <- subs[k]
    assign(".Random.seed", msg$seeds[[k]], envir = globalenv())
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
      population_var = shared$group_var
    )
    times[k] <- sum(proc.time()[c("user.self", "sys.self")]) - t0
    props[, k] <- c(out$proposal, out$ll)
    pm[[k]] <- out$pm_settings
    seeds[[k]] <- get(".Random.seed", envir = globalenv())
  }
  w_spent <- proc.time() - w_started
  list(props = props, pm = pm, seeds = seeds, times = times,
       t_recv = t_recv, t_done = Sys.time(),
       cpu = unname(w_spent[["user.self"]] + w_spent[["sys.self"]]),
       elapsed = unname(w_spent[["elapsed"]]),
       rejects = .emc_reject_delta(rejects_before, "particle"))
})

.emc_wpool_ll <- function(proposals, s, component = NULL, dynamic = FALSE,
                          varying = NULL) {
  pool <- .emc_pool_state$ll_pool
  if (is.null(pool) || !isTRUE(pool$alive) || pool$n <= 1L) return(NULL)
  if (!is.matrix(proposals)) return(NULL)
  n <- nrow(proposals)
  if (n <= pool$n) return(NULL)
  .emc_pool_state$ll_dynamic <- isTRUE(dynamic) && !is.null(pool$done) &&
    n >= .EMC_LL_DYN_MIN_ROWS * pool$n
  if (.emc_pool_state$ll_dynamic) {
    return(.emc_wpool_ll_dynamic(pool, proposals, s, component, n, varying))
  }
  idx <- .split_work_indices(n, pool$n)
  n_workers <- max(idx)

  sent <- logical(n_workers)
  for (w in seq_len(n_workers)) {
    sent[w] <- .emc_wpool_request(
      pool, w, list(kind = "ll", s = s, component = component,
                    varying = varying,
                    proposals = proposals[idx == w, , drop = FALSE]))
    if (!sent[w]) {
      .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
      break
    }
  }

  out <- numeric(n)
  ok <- all(sent)
  deadline <- .emc_wpool_deadline()
  for (w in which(sent)) {
    res <- .emc_wpool_reply(pool, w, deadline)
    if (is.null(res) || !is.null(res$failed) ||
        length(res$ll) != sum(idx == w)) {
      ok <- FALSE
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

.EMC_LL_DYN_MIN_ROWS <- 4L
.EMC_LL_DYN_GRAIN <- 4

.emc_wpool_ll_dynamic <- function(pool, proposals, s, component, n,
                                  varying = NULL) {
  k <- pool$n
  out <- numeric(n)
  pending <- vector("list", k)
  next_row <- 1L

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
  issued <- integer(k)
  send <- function(w, rows) {
    issued[w] <<- .emc_wpool_next_request()
    ok <- .emc_wpool_request(
      pool, w, list(kind = "ll", s = s, component = component,
                    notify = TRUE, w = w, generation = generation,
                    request = issued[w], varying = varying,
                    proposals = proposals[rows, , drop = FALSE]))
    if (!ok) {
      .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
    }
    ok
  }
  outstanding <- function() which(!vapply(pending, is.null, logical(1)))
  deadline <- .emc_wpool_deadline()
  await <- function() {
    repeat {
      rec <- .emc_wpool_await_record(pool, outstanding, deadline,
                                     on_start = function(rec) NULL)
      if (!identical(rec$code, "record")) return(NA_integer_)
      if (identical(rec$status, .EMC_WPOOL_ST_START)) next
      w <- rec$w
      if (is.na(w) || w < 1L || w > k) return(NA_integer_)
      if (!identical(rec$request, issued[w])) next
      return(w)
    }
  }

  fail <- function() {
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

.emc_ll_pool_set <- function(pool, s, reset = TRUE) {
  .emc_pool_state$ll_pool <- pool
  .emc_pool_state$ll_subject <- s
  if (reset) .emc_ll_route_reset()
  invisible(NULL)
}

.emc_ll_pool_clear <- function() {
  .emc_pool_state$ll_pool <- NULL
  .emc_pool_state$ll_subject <- NULL
  .emc_pool_state$ll_route <- NULL
  invisible(NULL)
}

.EMC_LL_PROBE_N <- 6L
.EMC_LL_REPROBE <- 200L
.EMC_LL_OBVIOUS <- 0.5

.emc_ll_route_reset <- function() {
  .emc_pool_state$ll_route <- list(mode = "probe", t_serial = numeric(0),
                                   t_pool = numeric(0), since = 0L,
                                   split = "probe", t_static = numeric(0),
                                   t_dynamic = numeric(0), since_split = 0L,
                                   last_route = NULL)
  invisible(NULL)
}

.emc_ll_route_use_pool <- function() {
  st <- .emc_pool_state$ll_route
  if (is.null(st)) { .emc_ll_route_reset(); st <- .emc_pool_state$ll_route }
  switch(st$mode,
         probe = length(st$t_pool) < length(st$t_serial),
         pool = TRUE,
         FALSE)
}

.emc_ll_route_use_dynamic <- function() {
  if (!isTRUE(getOption("emc2.ll_queue", TRUE))) return(FALSE)
  st <- .emc_pool_state$ll_route
  if (is.null(st) || !identical(st$mode, "pool")) return(FALSE)
  switch(st$split,
         probe = length(st$t_dynamic) < length(st$t_static),
         dynamic = TRUE,
         FALSE)
}

.emc_ll_route_record <- function(pooled, seconds, n_particles, dynamic = FALSE,
                                  route = NULL) {
  st <- .emc_pool_state$ll_route
  if (is.null(st) || !is.finite(seconds) || n_particles <= 0L) return(invisible(NULL))
  per <- seconds / n_particles
  if (is.null(route)) route <- if (pooled) "persistent_pool" else "serial"
  st$last_route <- route
  if (pooled) st$t_pool <- c(st$t_pool, per) else st$t_serial <- c(st$t_serial, per)
  if (pooled && identical(st$mode, "pool")) {
    if (dynamic) st$t_dynamic <- c(st$t_dynamic, per)
    else st$t_static <- c(st$t_static, per)
    if (identical(st$split, "probe")) {
      if (length(st$t_dynamic) >= .EMC_LL_PROBE_N &&
          length(st$t_static) >= .EMC_LL_PROBE_N) {
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
      st$mode <- if (stats::median(st$t_pool) < stats::median(st$t_serial))
        "pool" else "serial"
      st$since <- 0L
    }
  } else {
    st$since <- st$since + 1L
    if (st$since >= .EMC_LL_REPROBE) {
      st <- list(mode = "probe", t_serial = numeric(0), t_pool = numeric(0),
                 since = 0L, split = st$split, t_static = st$t_static,
                 t_dynamic = st$t_dynamic, since_split = st$since_split,
                 last_route = st$last_route)
    }
  }
  .emc_pool_state$ll_route <- st
  invisible(NULL)
}


.emc_wpool_iter <- function(pool, ctx, part, pars, group_chol, pm_settings,
                            prev_ll, seeds) {
  profile <- isTRUE(getOption("emc2.sampler_profile", FALSE))
  iter_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  n_subjects <- length(pm_settings)
  props <- matrix(0, ctx$n_pars + 1L, n_subjects)
  times <- numeric(n_subjects)

  population_mu <- vapply(seq_len(n_subjects), function(s) {
    as.numeric(get_group_level(pars, s, ctx$type)$mu)
  }, numeric(ctx$n_pars))
  rownames(population_mu) <- rownames(pars$alpha)
  alpha <- pars$alpha
  shared_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  shared <- list(group_var = group_chol$ref, idx_list = group_chol$idx_list,
                 marginal_idx = group_chol$marginal_idx,
                 include_full = group_chol$include_full,
                 group_chol = group_chol)
  was_alive <- isTRUE(pool$alive)
  shared_raw <- if (was_alive) {
    serialize(shared[c("group_var", "idx_list", "marginal_idx", "include_full")], NULL)
  } else raw(0)
  shared_file <- NULL
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
  watch <- !is.null(pool$done)
  generation <- .emc_wpool_one_int(pool$generation)
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
    send_deadline <- .emc_wpool_deadline()
    for (w in seq_len(pool$n)) {
      if (!length(part[[w]])) next
      sent[w] <- .emc_wpool_request(pool, w, msgs[[w]], send_deadline)
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
  store <- function(w, res, source = "worker") {
    subs <- part[[w]]
    props[, subs] <<- res$props
    times[subs] <<- res$times
    pm_settings[subs] <<- res$pm
    seeds[subs] <<- res$seeds
    if (!is.null(res$rejects)) rejects <<- rejects + res$rejects
    if (identical(source, "master") && was_alive) {
      fallback_workers <<- fallback_workers + 1L
      fallback_subjects <<- fallback_subjects + length(subs)
    }
    if (!profile) return(invisible(NULL))
    finish_at[w] <<- proc.time()[["elapsed"]] - iter_started
    if (!is.null(res$cpu)) worker_cpu[w] <<- res$cpu
    if (!is.null(res$elapsed)) worker_elapsed[w] <<- res$elapsed
    if (identical(source, "worker") && !is.null(res$t_recv) &&
        !is.na(dispatch_wall[w])) {
      queue_delay[w] <<- as.numeric(difftime(res$t_recv, dispatch_wall[w],
                                             units = "secs"))
    }
    invisible(NULL)
  }
  recompute <- function(w) .emc_wpool_compute(msgs[[w]], ctx, shared = shared)

  pending <- sent & vapply(part, length, integer(1)) > 0L
  picked_up <- rep(FALSE, pool$n)
  if (watch && any(pending)) {
    deadline <- .emc_wpool_deadline()
    note_start <- function(rec) {
      w <- rec$w
      if (!identical(rec$request, request)) {
        .emc_pool_state$stale_records <-
          (if (is.null(.emc_pool_state$stale_records)) 0L
           else .emc_pool_state$stale_records) + 1L
      } else if (!is.na(w) && w >= 1L && w <= pool$n) {
        picked_up[w] <<- TRUE
      }
      invisible(NULL)
    }
    while (any(pending)) {
      rec <- .emc_wpool_await_record(pool, function() which(pending), deadline,
                                     on_start = note_start)
      w <- if (identical(rec$code, "record")) rec$w else NA_integer_
      if (identical(rec$code, "record") &&
          !identical(rec$request, request)) {
        .emc_pool_state$stale_records <-
          (if (is.null(.emc_pool_state$stale_records)) 0L
           else .emc_pool_state$stale_records) + 1L
        next
      }
      if (identical(rec$code, "record") &&
          identical(rec$status, .EMC_WPOOL_ST_START) &&
          !is.na(w) && w >= 1L && w <= pool$n) {
        picked_up[w] <- TRUE
        next
      }
      if (is.na(w) || w < 1L || w > pool$n || !pending[w]) {
        if (identical(rec$code, "timeout")) {
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
          .emc_wpool_terminate_jobs(pool$jobs[stuck], wait = TRUE,
                                    terminate = TRUE)
        } else if (identical(rec$code, "dead") || is.na(w)) {
          .emc_wpool_degraded("worker gave no reply")
        } else {
          .emc_wpool_degraded("lost alignment on the worker completion channel")
        }
        pool$alive <- FALSE
        .emc_wpool_retire(pool)
        for (k in which(pending)) store(k, recompute(k), source = "master")
        pending[] <- FALSE
        break
      }
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
    for (w in seq_len(pool$n)) {
      if (length(part[[w]]) && !sent[w]) store(w, recompute(w), source = "master")
    }
  } else {
    deadline <- .emc_wpool_deadline()
    for (w in seq_len(pool$n)) {
      subs <- part[[w]]
      if (!length(subs)) next
      res <- if (sent[w]) .emc_wpool_reply(pool, w, deadline) else NULL
      source <- "worker"
      if (is.null(res) || !is.null(res$failed)) {
        if (sent[w] && is.null(res)) {
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
  worker_load <- vapply(part, function(subs) sum(times[subs]), numeric(1))
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

.emc_range_or_na <- function(x, f) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else as.numeric(f(x))
}


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

.emc_subject_cost <- function(data) {
  vapply(data, function(d) {
    if (is.data.frame(d)) nrow(d) else sum(vapply(d, function(x)
      if (is.data.frame(x)) nrow(x) else 0L, numeric(1)))
  }, numeric(1))
}

.emc_l_ecuyer_seed <- function() {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  caller_seed <- old_seed
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", caller_seed, envir = globalenv())
    else suppressWarnings(rm(".Random.seed", envir = globalenv()))
  }, add = TRUE)
  if (!had_seed) stats::runif(1)
  s <- sample.int(.Machine$integer.max, 1L)
  caller_seed <- get(".Random.seed", envir = globalenv())
  set.seed(s, kind = "L'Ecuyer-CMRG")
  get(".Random.seed", envir = globalenv())
}

.emc_subject_streams <- function(n_subjects, seed = NULL) {
  if (is.null(seed)) seed <- .emc_l_ecuyer_seed()
  if (!is.numeric(seed) || length(seed) != 7L || seed[[1L]] != 10407L) {
    stop("seed must be an L'Ecuyer-CMRG RNG state", call. = FALSE)
  }
  streams <- vector("list", n_subjects)
  for (s in seq_len(n_subjects)) {
    seed <- parallel::nextRNGSubStream(seed)
    streams[[s]] <- seed
  }
  streams
}

.emc_chain_rng <- function(chain_seed, n_subjects) {
  gibbs <- parallel::nextRNGSubStream(chain_seed)
  subjects <- .emc_subject_streams(n_subjects, seed = gibbs)
  list(chain = chain_seed, gibbs = gibbs, subjects = subjects)
}

.emc_rng_valid <- function(x, n_subjects = NULL) {
  if (!is.list(x) || is.null(x$chain) || is.null(x$gibbs) ||
      is.null(x$subjects) || !is.list(x$subjects)) return(FALSE)
  state_ok <- function(s) is.numeric(s) && length(s) == 7L &&
    identical(as.integer(s[[1L]]), 10407L)
  state_ok(x$chain) && state_ok(x$gibbs) &&
    (is.null(n_subjects) || length(x$subjects) == n_subjects) &&
    all(vapply(x$subjects, state_ok, logical(1)))
}

.emc_prepare_chain_rng <- function(emc) {
  if (!length(emc)) return(emc)
  n_subjects <- suppressWarnings(as.integer(emc[[1L]]$n_subjects))
  if (length(n_subjects) != 1L || is.na(n_subjects) || n_subjects < 0L) {
    return(emc)
  }
  if (all(vapply(emc, function(x) .emc_rng_valid(x$rng, n_subjects), logical(1)))) {
    return(emc)
  }
  base <- .emc_l_ecuyer_seed()
  for (i in seq_along(emc)) {
    base <- parallel::nextRNGStream(base)
    emc[[i]]$rng <- .emc_chain_rng(base, n_subjects)
  }
  emc
}

.emc_prepare_sampler_rng <- function(sampler) {
  n_subjects <- suppressWarnings(as.integer(sampler$n_subjects))
  if (length(n_subjects) != 1L || is.na(n_subjects) || n_subjects < 0L) {
    return(sampler)
  }
  if (!.emc_rng_valid(sampler$rng, n_subjects)) {
    base <- .emc_l_ecuyer_seed()
    sampler$rng <- .emc_chain_rng(base, n_subjects)
  }
  sampler
}

.emc_set_sampler_rng <- function(sampler) {
  if (!is.null(sampler$rng$gibbs)) {
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  }
  sampler
}
