# ===========================================================================
# Persistent worker pool for the particle step
# ===========================================================================
#
# `run_stage()` forks a fresh set of `cores_per_chain` workers on *every* MCMC
# iteration, via `mcmapply`.  A block of 40 iterations therefore pays 40
# fork/join cycles per chain, and each cycle has two costs:
#
#   * The fork itself is serial in the chain process and scales with its heap.
#     Measured on a 330 MB heap, `mclapply` over 8 workers costs 33 ms of pure
#     fork/join with no work at all.  A real fit's heap is larger.
#   * `mc.preschedule` splits the subjects into contiguous, equal-*count*
#     chunks.  Subjects do not cost the same -- trial counts differ, and the
#     adaptive particle number differs -- so the iteration ends when the
#     unluckiest chunk ends, and every other core idles until it does.
#
# Together those are what makes core usage sag between iterations rather than
# staying flat.  They are a per-iteration cost, so unlike inter-chain straggle
# they cannot be recovered by handing cores around at block boundaries.
#
# Here the workers are started once per block and kept alive across the whole
# block, so process startup is paid once instead of `iter` times, and the
# subject partition is rebalanced every iteration from measured times.  On a
# real 112-subject, 17-parameter RDMSWTN fit (3 chains x 8 cores, 40-iteration
# blocks) this takes a block from 37.3 s to 32.4 s.
#
# Transport.  A socket cluster is the obvious way to keep workers alive and it
# does not work: `makeForkCluster` costs 44 ms per dispatch round on the same
# probe -- worse than the fork it was meant to replace, because R's socket
# connections stall on delayed ACK.  Named pipes cost 1.0 ms for the same
# round.  One clean R template reads the block-constant state (the per-subject
# data, the design, the Cholesky caches) once at startup and forks the workers
# from that history-free heap.  The only thing crossing the pipe each iteration
# is the group-level draw and per-subject bookkeeping -- about 5 kB on a small
# model.  That grows
# with the square of the parameter count, though, since the group covariance
# travels with it: ~35 kB at 24 parameters, ~90 kB at 42.  See
# `.emc_wpool_spawn()` for why the writer has to block rather than fail once
# that passes the 64 kB pipe buffer.
#
#     mclapply(8) fork+join   33.5 ms
#     makeForkCluster round   44.1 ms
#     named-pipe pool round    1.0 ms
#
# Scheduling.  Each worker gets one message per iteration containing a set of
# subjects, so dispatch stays at one round trip regardless of subject count.
# The partition is longest-processing-time-first over the *previous*
# iteration's measured per-subject times (falling back to dadm row counts for
# the first iteration).  Per-subject cost is stable from one iteration to the
# next, so a one-iteration lag costs nothing and this tracks the adaptive
# particle counts as they drift.
#
# RNG.  Every subject carries its own L'Ecuyer stream, held in the master and
# advanced by whichever worker handled it.  Draws are therefore independent of
# the number of workers and of how subjects happen to be partitioned.  The
# `mcmapply` fallback below does not have that property -- it derives its
# children's streams from `mc.cores` -- so a pooled fit is reproducible in a
# way that one is not.
#
# The clean template is important.  A worker forked from the chain inherits the
# chain's complete sample history.  R's GC subsequently writes to object
# headers throughout that inherited heap, turning shared pages into a private
# copy in every worker.  The template has never held those samples, so history
# cannot enter any worker's address space; its workers still share the large
# immutable context initially.  The old chain-fork backend remains as a
# fallback for source-loaded packages and custom external pointers.
#
# The pool is the normal path.  `.emc_wpool_start()` returns NULL where it
# cannot work (Windows, no mkfifo, a single worker), and `run_stage()` then
# falls back to the original `mcmapply` call.

# Where a failed worker's message is kept, so a pool that has quietly fallen
# back to serial recomputation can be diagnosed after the fact.
.emc_pool_state <- new.env(parent = emptyenv())

# Losing the pool costs the block its parallelism and nothing else -- the
# results are identical, just computed one subject at a time in the master.
# That is exactly why it has to be said out loud: a silent 8x slowdown looks
# like a slow model, not like a bug.  `warning()` alone would not do it, since
# a chain normally runs inside an `mclapply` fork whose deferred warnings are
# discarded when it exits; `immediate. = TRUE` writes to the inherited stderr
# there and then.
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

# A clean Rscript worker can load only an installed namespace.  `pkgload`
# namespaces point at the source tree and have no Meta/package.rds; keep the
# fork backend for that development case.  It is also the escape hatch for a
# custom trend kernel: an external pointer cannot survive serialisation into a
# fresh process, whereas it remains valid across fork().
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

# Design constructors return a tiny function whose enclosing frame can retain
# data and intermediate objects that the likelihood never uses.  That frame is
# harmless under fork but needlessly inflates the clean template's context
# file (15.8 MB for Annika).  Reclose each model around only its realised model
# specification.  The likelihood-facing contract remains the same: call the
# function and receive the model list.
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
  length(pid) == 1L && !is.na(pid) &&
    isTRUE(tryCatch(tools::pskill(pid, 0L), error = function(e) FALSE))
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

# Fork one batch of workers on `idx` and connect to them.  Shared by the
# initial start and by a later grow, which differ only in which indices they
# claim.  Returns NULL if any worker fails to come up, having cleaned up after
# itself.
#
# The request connection is opened twice on purpose.  Opening the write end of
# a FIFO blocks until a reader appears, so the handshake has to be
# non-blocking -- but a *non-blocking* write fails outright once the payload
# exceeds the 64 kB pipe buffer, which a group covariance does at about 35
# parameters.  A second, blocking writer on the same FIFO does not signal EOF
# to the reader, so it can take over from the handshake handle and then wait
# for the worker to drain rather than failing.  The protocol reads every reply
# before the next send, so no worker is ever mid-compute when the master
# writes and a blocking write cannot deadlock.
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
      # Leave no half-built FIFO behind, or a later grow onto the same indices
      # would fail on mkfifo and the pool could never take those cores.
      unlink(c(req, ans))
    }
  }, add = TRUE)

  # fork()/process creation fails with EAGAIN once the process limit is reached,
  # and this design
  # makes that reachable: donation lets one chain claim the whole budget, and a
  # recycle forks a fresh set while the old one is still being reaped.  Every
  # other failure in this file returns NULL and lets the caller fall back, so
  # this one must too.  Letting the error escape would kill the chain outright
  # -- and from .emc_wpool_recycle(), where the old workers have already been
  # stopped, it would take the whole block's work with it.
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
      # Poll until the worker has opened its end.  A worker is normally there
      # within a millisecond, so back off from well under that rather than
      # sleeping a flat 10 ms: this loop runs once per worker, serially, and is
      # paid again on every grow and every recycle, so at eight workers the flat
      # version could spend ~80 ms per spawn waiting on processes that were
      # already up.  (It is not the whole of a spawn -- forking the processes
      # and reaping the old ones dominate -- but it is the part that was pure
      # sleeping.)
      f <- NULL
      deadline <- Sys.time() + 20  # same ceiling as before
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
      # The reader is attached now, so this returns at once; then drop the
      # handshake handle and keep the blocking one for all traffic.
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

# Start one pristine R process for the chain.  It loads the package and context
# once, then fork()s all persistent workers from that small template.  The
# workers therefore share immutable data and compiled package state, but the
# chain's sampler history has never existed in their ancestry.  Keeping the
# template idle also gives recycling a clean parent: replacing eight workers
# does not reload or reserialise the context.
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
  # The spawned template (and every chain worker forked from it) inherits this
  # command line, so lead with a self-identifying no-op assignment: `ps` and
  # `top` then show `EMC2-worker-pool[<dir>|master=<pid>]` instead of an
  # anonymous `readRDS(...)`.  The tag is sanitised because it is pasted into
  # the expression as a literal.
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

# `ctx` is everything constant for the whole block.  The legacy backend
# inherits it directly from the chain.  The clean backend serialises it once to
# the template, whose workers then inherit the loaded pages through fork().
.emc_wpool_start <- function(n_workers, ctx) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers <= 1L) return(NULL)
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
       ctx_file = ctx_file, lib = lib, template = template)
}

# Take over cores released by a chain that has finished its block.  Chains do
# not finish a block together -- the slowest takes ~30% longer than the fastest
# on real fits -- so without this the finished chains' cores idle until the
# block ends.  Growing the pool is safe precisely because every subject owns
# its stream: who computes it, and how many workers there are, cannot change a
# single draw.  (Feeding a timing-dependent core count to `mc.cores` would not
# be safe, since that seeds the children.)
.emc_wpool_grow <- function(pool, n_total, ctx) {
  n_total <- as.integer(n_total)
  if (is.null(pool) || !isTRUE(pool$alive) || is.na(n_total)) return(pool)
  n_new <- n_total - pool$n
  if (n_new <= 0L) return(pool)

  # A grow that cannot spawn must not be retried on every iteration.  The cheap
  # failure (mkfifo refuses) costs nothing to repeat, but the expensive one --
  # the workers fork and then never complete the handshake -- costs up to 20 s
  # *per worker*, which would dwarf the iteration the extra cores were meant to
  # speed up.  Back off exponentially instead and let the block run at its
  # current width; the cores were a donation, not the budget this chain needs.
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
    pool$grow_skip <- bitwShiftL(1L, min(fails, 6L))  # 2, 4, 8 ... 64 iterations
    return(pool)
  }

  list(n = n_total, dir = pool$dir, jobs = c(pool$jobs, spawned$jobs),
       wcs = c(pool$wcs, spawned$wcs), rcs = c(pool$rcs, spawned$rcs),
       alive = TRUE, grow_fails = 0L, grow_skip = 0L,
       backend = pool$backend, ctx_file = pool$ctx_file, lib = pool$lib,
       template = pool$template)
}

# Replace the workers with a fresh fork every `every` iterations.
#
# R's garbage collector writes to the header of each object it marks, and each
# write un-shares a page for good.  Under the legacy backend a worker therefore
# converges on a private copy of the entire chain heap, including history.
# Measured on a 250 MB chain heap with 8 workers, total physical memory (PSS,
# so sharing is already accounted for) against block length:
#
#     no pool (the old per-iteration mcmapply forks)   0.28 GB
#      10 iterations                                   0.75 GB
#      25 iterations                                   1.10 GB
#      50 iterations                                   2.36 GB
#     100 iterations                                   3.54 GB
#     200 iterations                                   3.83 GB  (~8 x the heap)
#
# The old path never showed this because its workers were forked per iteration
# and died before they had collected anything -- the cost arrived with their
# longer life, not with the pool itself.  Left alone it saturates at roughly
# n_workers x heap, which on several chains of a large fit is enough to reach
# the OOM killer overnight, and the donation path makes it worst exactly when
# one chain is left holding every core.
#
# The clean-template backend removes history from that ceiling entirely.  Its
# workers can still dirty their copies of the block context, so periodic
# replacement remains useful; replacement forks from the already-loaded idle
# template and does not repeat R startup, context serialisation, or context
# loading.  At 10 iterations most sharing is still intact.  The original
# chain-fork measurements over 100 iterations were:
#
#                       peak PSS   elapsed
#     250 MB heap, off     3.58 GB    8.5 s
#     250 MB heap, 10      1.01 GB    9.9 s   -72% memory, +17% time
#     small heap,   off    0.92 GB    8.2 s
#     small heap,   10     0.74 GB   10.1 s   -19% memory, +23% time
#
# Re-forking costs a roughly fixed ~0.19 s -- eight processes, the FIFO
# handshake, and waiting for the old workers to exit -- near enough independent
# of the heap (and no more for the clean template).  The period converts that
# into a share of the block on its own:
# an iteration here is 85 ms, so ten of them are 0.85 s and a recycle is +22%,
# whereas on a real fit at ~800 ms an iteration those ten are 8 s and the same
# recycle is +2%.  The overhead is therefore large as a fraction only where the
# whole block is already seconds long, and it shrinks exactly as the fits get
# big enough for the memory to matter.  Raise `emc2.worker_recycle` to trade
# more memory for less of it; 0 or less turns recycling off.
#
# A wall-clock floor was tried here instead and is the wrong control: divergence
# gets *faster* with a bigger heap, so a 250 MB fit was already fully diverged
# eight seconds in, and any floor long enough to matter skipped the case that
# needed it most.
#
# Safe for the same reason growing is: every subject's RNG stream lives in the
# master and travels in the message, so which process computes a subject, and
# how many processes there are, cannot change a draw.
.emc_wpool_recycle <- function(pool, i, every, ctx) {
  if (is.null(pool)) return(pool)
  if (is.na(every) || every <= 0L) return(pool)
  # Nothing has run yet on iteration 1, and pool$n <= 1 is the degenerate pool
  # that .emc_wpool_start() refuses anyway.
  if (i <= 1L || pool$n <= 1L) return(pool)
  if (((i - 1L) %% every) != 0L) return(pool)

  target <- pool$n
  # A pool that lost a worker mid-block used to stay lost for the rest of the
  # block: `alive = FALSE` was permanent, so the surviving workers idled and
  # every subject was recomputed in the master -- the silent 8x slowdown this
  # file warns about, paid for the remaining iterations even when the loss was
  # one transient write error.  A recycle point is where a fresh set is forked
  # anyway, so it is also the natural place to recover.  Back off, though: a
  # machine that genuinely cannot fork should not be asked every period.
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
    # The old workers are already gone, so there is nothing to keep running on.
    # Hand back a dead pool: the master then computes every subject itself,
    # which is slow but produces exactly the same numbers.
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
  # Try every connection even when one worker marked the aggregate pool dead;
  # the surviving workers are still blocked waiting for their shutdown token.
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
  if (identical(pool$backend, "spawn")) {
    .emc_wpool_template_stop(pool$template)
  }
  unlink(pool$dir, recursive = TRUE)
  invisible(NULL)
}

# Workers forked from the chain can eventually privatise the complete sampler
# heap, so the conservative historical period remains appropriate there.  The
# clean template has never held sampler history; measured on 140-subject blocks
# a ten-iteration period retained essentially no additional PSS but imposed a
# material fixed process-replacement cost.  Keep an occasional reset for dirty
# block-context pages without paying it four times in a 40-iteration block.
.emc_wpool_recycle_default <- function(pool) {
  if (!is.null(pool) && identical(pool$backend, "spawn")) 50L else 10L
}

# --- worker -----------------------------------------------------------------

.emc_wpool_serve <- function(req, ans, ctx, inherited = NULL) {
  for (cn in inherited) try(close(cn), silent = TRUE)
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)
  repeat {
    msg <- tryCatch(unserialize(rc), error = function(e) NULL)
    if (is.null(msg)) break                     # shutdown, or master went away
    out <- tryCatch(.emc_wpool_compute(msg, ctx),
                    error = function(e) list(failed = conditionMessage(e)))
    if (!tryCatch({ serialize(out, wc); flush(wc); TRUE },
                  error = function(e) FALSE)) break
  }
  invisible(NULL)
}

# Entry point for the clean, idle template process.  Only this process reads the
# context file.  Its children inherit the resulting pages from a heap that has
# never held sampler history; periodic replacement forks from the same pristine
# template again.
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

# Run `expr` with the caller's generator handed back exactly as it was found.
#
# Per-subject streams are installed on the *global* seed, which is private and
# disposable inside a forked worker but not in the master.  Both of the places
# that install them can end up running in the master -- the pool's fallback
# path, and `mclapply` at one core, which does not fork at all -- and either
# would otherwise leave the master sitting on some subject's stream, colliding
# with that subject on every later draw.  That is also what makes a one-core
# fit differ from a many-core one.
.emc_with_preserved_rng <- function(expr) {
  if (exists(".Random.seed", envir = globalenv())) {
    seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  expr
}

# The unit of work, shared by the workers and by the master's fallback path so
# that a broken pool computes exactly what a working one would have.
.emc_wpool_compute <- function(msg, ctx, shared = NULL) .emc_with_preserved_rng({
  # `shared` is the group-level draw, which reaches a worker as the raw bytes
  # .emc_wpool_iter() serialised once for every worker.  The master's fallback
  # path already holds the decoded object and passes it directly rather than
  # paying a pointless round trip through serialize/unserialize.
  if (is.null(shared)) shared <- unserialize(msg$shared)
  subs <- msg$subs
  props <- matrix(0, ctx$n_pars + 1L, length(subs))
  pm <- vector("list", length(subs))
  seeds <- vector("list", length(subs))
  times <- numeric(length(subs))
  for (k in seq_along(subs)) {
    s <- subs[k]
    assign(".Random.seed", msg$seeds[[k]], envir = globalenv())
    # CPU, not wall: on a shared machine a descheduled worker would otherwise
    # look expensive, and that error would feed straight back into the next
    # iteration's partition.
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
