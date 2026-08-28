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
  stat <- tryCatch(readLines(sprintf("/proc/%d/stat", as.integer(pid)), n = 1L),
                   error = function(e) character())
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
      deadline <- Sys.time() + 2
      while (any(vapply(pids, .emc_wpool_pid_alive, logical(1))) &&
             Sys.time() < deadline) Sys.sleep(0.005)
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
      serialize(list(command = "spawn", idx = idx, req = req, ans = ans),
                template$wc)
      flush(template$wc)
      reply <- unserialize(template$rc)
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
      deadline <- Sys.time() + 20  # startup timeout ceiling
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
        if (Sys.time() > deadline) break
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
      rcs[[w]] <- fifo(ans[w], "rb", blocking = TRUE)
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
  deadline <- Sys.time() + 20
  wait <- 0.0002
  repeat {
    pid <- .emc_wpool_job_pid(job)
    if (!is.na(pid)) {
      job$pid <- pid
      if (!.emc_wpool_pid_alive(pid)) break
    }
    f <- tryCatch(suppressWarnings(fifo(req, "wb", blocking = FALSE)),
                  error = function(e) NULL)
    if (!is.null(f) || Sys.time() > deadline) break
    Sys.sleep(wait)
    wait <- min(wait * 2, 0.01)
  }
  if (is.null(f)) return(NULL)
  wc <- fifo(req, "wb", blocking = TRUE)
  close(f)
  rc <- fifo(ans, "rb", blocking = TRUE)
  ready <- tryCatch(unserialize(rc), error = function(e) NULL)
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

  list(n = n_workers, dir = dir, jobs = spawned$jobs, wcs = spawned$wcs,
       rcs = spawned$rcs, alive = TRUE, backend = backend,
       ctx_file = ctx_file, lib = lib, template = template,
       done = .emc_wpool_done_open(dir))
}

# Completion FIFO shared by workers in the dynamic queue. Each worker writes its
# one-byte index before its reply, allowing the master to dispatch whichever
# worker finished rather than blocking on a particular pipe.
#
# A single-byte token cannot be split across short non-blocking reads, so the
# stream remains aligned. The reader is non-blocking to detect worker death;
# keeping a master writer open prevents spurious EOF. NULL disables the queue
# and selects the static split.
.EMC_WPOOL_MAX_DYN_WORKERS <- 255L

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
  list(rc = rc, wc = wc)
}

# Grow the pool into cores released by completed chains. Per-subject streams make
# changing worker count and assignment deterministic.
.emc_wpool_grow <- function(pool, n_total, ctx) {
  n_total <- as.integer(n_total)
  if (is.null(pool) || !isTRUE(pool$alive) || is.na(n_total)) return(pool)
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

  # Keep the existing completion FIFO so added workers join the same queue.
  list(n = n_total, dir = pool$dir, jobs = c(pool$jobs, spawned$jobs),
       wcs = c(pool$wcs, spawned$wcs), rcs = c(pool$rcs, spawned$rcs),
       alive = TRUE, grow_fails = 0L, grow_skip = 0L,
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

  if (identical(pool$backend, "spawn") && !is.null(pool$template)) {
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
      pool
    }
  } else {
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(target, ctx)
  }
  if (is.null(fresh)) {
    # Workers are already gone; return a dead pool so the master recomputes this
    # work serially with identical streams. Close the completion FIFO before
    # discarding it.
    if (!is.null(pool$done)) {
      for (cn in pool$done) try(close(cn), silent = TRUE)
    }
    if (!dead) .emc_wpool_degraded("could not re-fork the pool when recycling")
    fails <- if (is.null(pool$rebuild_fails)) 1L else pool$rebuild_fails + 1L
    return(list(n = target, dir = NULL, jobs = list(), wcs = list(),
                rcs = list(), alive = FALSE,
                rebuild_fails = fails, rebuild_skip = min(fails, 8L)))
  }
  fresh
}

.emc_wpool_stop_workers <- function(pool) {
  if (is.null(pool)) return(invisible(NULL))
  # Send shutdown to every worker; surviving workers may still await its token.
  for (w in seq_len(length(pool$wcs))) {
    try({ serialize(NULL, pool$wcs[[w]]); flush(pool$wcs[[w]]) }, silent = TRUE)
  }
  for (cn in c(pool$wcs, pool$rcs)) try(close(cn), silent = TRUE)
  if (!identical(pool$backend, "spawn")) {
    .emc_wpool_terminate_jobs(pool$jobs, wait = TRUE, terminate = TRUE)
  }
  invisible(NULL)
}

.emc_wpool_template_stop <- function(template) {
  if (is.null(template)) return(invisible(NULL))
  if (isTRUE(template$alive)) {
    try({ serialize(NULL, template$wc); flush(template$wc) }, silent = TRUE)
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
    for (cn in pool$done) try(close(cn), silent = TRUE)
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

.emc_wpool_serve <- function(req, ans, ctx, inherited = NULL) {
  for (cn in inherited) try(close(cn), silent = TRUE)
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  done <- NULL
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE)
            if (!is.null(done)) try(close(done), silent = TRUE) }, add = TRUE)
  repeat {
    msg <- tryCatch(unserialize(rc), error = function(e) NULL)
    if (is.null(msg)) break                     # shutdown or master disconnect
    out <- tryCatch(.emc_wpool_compute(msg, ctx),
                    error = function(e) list(failed = conditionMessage(e)))
    # Announce completion before writing the reply: a large reply can block on
    # the pipe, while the master waits for this token to know which reply to read.
    # Open the channel lazily so ordinary pool messages remain unchanged.
    if (isTRUE(msg$notify)) {
      if (is.null(done)) {
        done <- tryCatch(fifo(file.path(dirname(req), "done"), "wb",
                              blocking = TRUE), error = function(e) NULL)
      }
      if (!is.null(done)) {
        if (!tryCatch({ writeBin(as.raw(msg$w), done); flush(done); TRUE },
                      error = function(e) FALSE)) break
      }
    }
    if (!tryCatch({ serialize(out, wc); flush(wc); TRUE },
                  error = function(e) FALSE)) break
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
  serialize(list(ready = TRUE, pid = Sys.getpid()), wc)
  flush(wc)
  jobs <- list()
  repeat {
    command <- tryCatch(unserialize(rc), error = function(e) NULL)
    if (is.null(command)) break
    reply <- tryCatch({
      if (!identical(command$command, "spawn")) stop("unknown template command")
      keys <- as.character(command$idx)
      old <- jobs[keys]
      old <- old[!vapply(old, is.null, logical(1))]
      if (length(old)) parallel::mccollect(old, wait = TRUE)
      jobs[keys] <- NULL
      made <- vector("list", length(keys))
      for (w in seq_along(keys)) {
        made[[w]] <- local({
          i <- w
          parallel::mcparallel(
            .emc_wpool_serve(command$req[i], command$ans[i], ctx,
                             inherited = list(rc, wc)),
            detached = FALSE
          )
        })
      }
      names(made) <- keys
      jobs[keys] <- made
      list(ok = TRUE, pids = vapply(made, function(x) x$pid, integer(1)))
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    if (!tryCatch({ serialize(reply, wc); flush(wc); TRUE },
                  error = function(e) FALSE)) break
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

  sent <- tryCatch({
    serialize(list(kind = "loo", notify = FALSE, w = 1L), pool$wcs[[1L]])
    flush(pool$wcs[[1L]])
    TRUE
  }, error = function(e) FALSE)
  if (!sent) return(NULL)
  res <- tryCatch(unserialize(pool$rcs[[1L]]), error = function(e) NULL)
  if (is.null(res) || !is.null(res$failed) || is.null(res$value)) return(NULL)
  if (length(res$warnings) && exists(".loo_warn_store")) {
    .loo_warn_store$msgs <- c(.loo_warn_store$msgs, res$warnings)
  }
  res$value
}

.emc_wpool_compute_particle <- function(msg, ctx, shared = NULL) .emc_with_preserved_rng({
  # Workers receive the shared group draw as serialised bytes; the master
  # fallback already has the decoded object.
  if (is.null(shared)) shared <- unserialize(msg$shared)
  subs <- msg$subs
  props <- matrix(0, ctx$n_pars + 1L, length(subs))
  pm <- vector("list", length(subs))
  seeds <- vector("list", length(subs))
  times <- numeric(length(subs))
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
      chol_cache = ctx$chol_caches[[s]], group_chol = shared$group_chol,
      current_alpha = msg$alpha[, k],
      population_mu = msg$population_mu[, k],
      # build_group_chol_cache() retains the exact covariance in `ref` for its
      # cache-validity check.  Reuse it here instead of sending the same P x P
      # matrix a second time in the shared payload.
      population_var = shared$group_chol$ref
    )
    times[k] <- sum(proc.time()[c("user.self", "sys.self")]) - t0
    props[, k] <- c(out$proposal, out$ll)
    pm[[k]] <- out$pm_settings
    seeds[[k]] <- get(".Random.seed", envir = globalenv())
  }
  list(props = props, pm = pm, seeds = seeds, times = times)
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
  .emc_pool_state$ll_dynamic <- isTRUE(dynamic) && !is.null(pool$done) &&
    pool$n <= .EMC_WPOOL_MAX_DYN_WORKERS &&
    n >= .EMC_LL_DYN_MIN_ROWS * pool$n
  if (.emc_pool_state$ll_dynamic) {
    return(.emc_wpool_ll_dynamic(pool, proposals, s, component, n))
  }
  idx <- .split_work_indices(n, pool$n)
  n_workers <- max(idx)

  sent <- logical(n_workers)
  for (w in seq_len(n_workers)) {
    sent[w] <- tryCatch({
      serialize(list(kind = "ll", s = s, component = component,
                     proposals = proposals[idx == w, , drop = FALSE]),
                pool$wcs[[w]])
      flush(pool$wcs[[w]])
      TRUE
    }, error = function(e) { .emc_wpool_degraded(conditionMessage(e)); FALSE })
    if (!sent[w]) break
  }

  # Drain every sent worker after a failed round so stale replies cannot answer
  # the next call.
  out <- numeric(n)
  ok <- all(sent)
  for (w in which(sent)) {
    res <- tryCatch(unserialize(pool$rcs[[w]]), error = function(e) NULL)
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
  send <- function(w, rows) {
    tryCatch({
      serialize(list(kind = "ll", s = s, component = component,
                     notify = TRUE, w = w,
                     proposals = proposals[rows, , drop = FALSE]),
                pool$wcs[[w]])
      flush(pool$wcs[[w]])
      TRUE
    }, error = function(e) { .emc_wpool_degraded(conditionMessage(e)); FALSE })
  }
  # Poll the non-blocking completion channel so worker death cannot hang the
  # master. Spin briefly, then back off while checking worker liveness.
  await <- function() {
    idle <- 0L
    wait <- 0
    repeat {
      tok <- tryCatch(readBin(pool$done$rc, "raw", 1L),
                      error = function(e) raw(0))
      if (length(tok)) return(as.integer(tok[1L]))
      idle <- idle + 1L
      if ((idle %% 64L) == 0L) {
        for (w in which(!vapply(pending, is.null, logical(1)))) {
          pid <- .emc_wpool_job_pid(pool$jobs[[w]])
          if (!is.null(pid) && !.emc_wpool_pid_alive(pid)) return(NA_integer_)
        }
      }
      if (idle > 16L) { Sys.sleep(wait); wait <- min(max(wait * 2, 5e-5), 2e-3) }
    }
  }

  fail <- function() {
    # Drain outstanding completion tokens and replies before abandoning the
    # queue, preventing stale answers from contaminating the next call.
    n_out <- sum(!vapply(pending, is.null, logical(1)))
    for (i in seq_len(n_out)) {
      w <- await()
      if (is.na(w)) break
      tryCatch(unserialize(pool$rcs[[w]]), error = function(e) NULL)
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
    res <- tryCatch(unserialize(pool$rcs[[w]]), error = function(e) NULL)
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
  # by assignment, so collectively they cross one pipe rather than every pipe.
  # This retains per-iteration LPT balancing without its former O(workers * N)
  # subject-state broadcast.  Gibbs-only fields (tvinv, a_half, factor state,
  # etc.) never enter a particle message.
  population_mu <- vapply(seq_len(n_subjects), function(s) {
    as.numeric(get_group_level(pars, s, ctx$type)$mu)
  }, numeric(ctx$n_pars))
  rownames(population_mu) <- rownames(pars$alpha)
  alpha <- pars$alpha
  shared_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  shared <- list(group_chol = group_chol)
  shared_raw <- serialize(shared, NULL)
  shared_elapsed <- if (profile) {
    proc.time()[["elapsed"]] - shared_started
  } else NA_real_
  msgs <- lapply(part, function(subs) {
    list(subs = subs, shared = shared_raw,
         alpha = alpha[, subs, drop = FALSE],
         population_mu = population_mu[, subs, drop = FALSE],
         pm = pm_settings[subs], prev_ll = prev_ll[subs], seeds = seeds[subs])
  })
  was_alive <- isTRUE(pool$alive)
  sent <- rep(FALSE, pool$n)
  send_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  if (was_alive) {
    for (w in seq_len(pool$n)) {
      if (!length(part[[w]])) next
      sent[w] <- tryCatch({
        serialize(msgs[[w]], pool$wcs[[w]]); flush(pool$wcs[[w]]); TRUE
      }, error = function(e) {
        .emc_wpool_degraded(conditionMessage(e)); FALSE
      })
      if (!sent[w]) pool$alive <- FALSE
    }
  }
  send_elapsed <- if (profile) {
    proc.time()[["elapsed"]] - send_started
  } else NA_real_
  receive_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  for (w in seq_len(pool$n)) {
    subs <- part[[w]]
    if (!length(subs)) next
    res <- if (sent[w]) {
      tryCatch(unserialize(pool$rcs[[w]]), error = function(e) NULL)
    } else NULL
    # A dead or erroring worker must not lose an iteration: recompute its share
    # here.  Same function, same streams, so the result is identical.
    if (is.null(res) || !is.null(res$failed)) {
      if (sent[w] && is.null(res)) {
        # No reply means the transport or the worker itself is gone, so there
        # is nothing left to send to and the rest of the block is serial.
        .emc_wpool_degraded("worker gave no reply")
        pool$alive <- FALSE
      }
      # A reported `failed` is different: the worker caught an error in the
      # work and is still listening.  Recompute this share, but do not condemn
      # the whole block over one subject's bad iteration -- if it is
      # deterministic the master's own call raises it properly.
      res <- .emc_wpool_compute(msgs[[w]], ctx, shared = shared)
    }
    props[, subs] <- res$props
    times[subs] <- res$times
    pm_settings[subs] <- res$pm
    seeds[subs] <- res$seeds
  }
  receive_elapsed <- if (profile) {
    proc.time()[["elapsed"]] - receive_started
  } else NA_real_
  private_bytes <- if (profile) {
    sum(vapply(msgs, function(msg) {
      msg$shared <- NULL
      length(serialize(msg, NULL))
    }, integer(1)))
  } else NA_real_
  active_workers <- sum(vapply(part, length, integer(1)) > 0L)
  list(props = props, pm_settings = pm_settings, seeds = seeds, times = times,
       alive = pool$alive,
       profile = if (profile) list(
         elapsed = proc.time()[["elapsed"]] - iter_started,
         shared_serialize = shared_elapsed,
         send = send_elapsed,
         receive = receive_elapsed,
         worker_max = if (length(times)) max(times) else 0,
         worker_sum = sum(times),
         shared_bytes = length(shared_raw),
         private_bytes = private_bytes,
         wire_bytes = active_workers * length(shared_raw) + private_bytes
       ) else NULL)
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
