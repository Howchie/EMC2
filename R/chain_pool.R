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
# Here the workers are forked once per block and kept alive across the whole
# block, so the fork is paid once instead of `iter` times, and the subject
# partition is rebalanced every iteration from measured times.  On a real
# 112-subject, 17-parameter RDMSWTN fit (3 chains x 8 cores, 40-iteration
# blocks) this takes a block from 37.3 s to 32.4 s.
#
# Transport.  A socket cluster is the obvious way to keep workers alive and it
# does not work: `makeForkCluster` costs 44 ms per dispatch round on the same
# probe -- worse than the fork it was meant to replace, because R's socket
# connections stall on delayed ACK.  Named pipes cost 1.0 ms for the same
# round, and forked workers inherit the block-constant state (the per-subject
# data, the design, the Cholesky caches) rather than receiving a copy of it,
# so the only thing crossing the pipe each iteration is the group-level draw
# and the per-subject bookkeeping -- about 5 kB on a small model.  That grows
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
.emc_wpool_spawn <- function(dir, idx, ctx) {
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
      for (job in live) try(tools::pskill(job$pid), silent = TRUE)
      try(parallel::mccollect(live, wait = TRUE), silent = TRUE)
      # Leave no half-built FIFO behind, or a later grow onto the same indices
      # would fail on mkfifo and the pool could never take those cores.
      unlink(c(req, ans))
    }
  }, add = TRUE)

  for (w in seq_len(n)) {
    jobs[[w]] <- local({
      i <- w
      parallel::mcparallel(.emc_wpool_serve(req[i], ans[i], ctx), detached = FALSE)
    })
  }

  tryCatch({
    for (w in seq_len(n)) {
      f <- NULL
      for (poll in 1:2000) { # 20 seconds max
        res <- parallel::mccollect(jobs[[w]], wait = FALSE, timeout = 0)
        if (!is.null(res)) stop("Worker process died before pool initialization")
        f <- tryCatch(suppressWarnings(fifo(req[w], "wb", blocking = FALSE)),
                      error = function(e) NULL)
        if (!is.null(f)) break
        Sys.sleep(0.01)
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
  }, error = function(e) NULL, interrupt = function(e) NULL)

  if (!success) return(NULL)
  list(jobs = jobs, wcs = wcs, rcs = rcs)
}

# `ctx` is everything constant for the whole block.  It is captured by the
# child expression, so on Unix the workers inherit it through the fork; it is
# never serialised.
.emc_wpool_start <- function(n_workers, ctx) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers <= 1L) return(NULL)
  if (Sys.info()[["sysname"]] == "Windows") return(NULL)
  if (!nzchar(Sys.which("mkfifo"))) return(NULL)

  dir <- tempfile("emc_wpool_")
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) return(NULL)

  # Per pool, not per session: a later block that degrades must say so too.
  .emc_pool_state$warned <- FALSE
  spawned <- .emc_wpool_spawn(dir, seq_len(n_workers), ctx)
  if (is.null(spawned)) {
    unlink(dir, recursive = TRUE)
    return(NULL)
  }

  list(n = n_workers, dir = dir, jobs = spawned$jobs, wcs = spawned$wcs,
       rcs = spawned$rcs, alive = TRUE)
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

  spawned <- .emc_wpool_spawn(pool$dir, pool$n + seq_len(n_new), ctx)
  if (is.null(spawned)) return(pool)

  list(n = n_total, dir = pool$dir, jobs = c(pool$jobs, spawned$jobs),
       wcs = c(pool$wcs, spawned$wcs), rcs = c(pool$rcs, spawned$rcs),
       alive = TRUE)
}

.emc_wpool_stop <- function(pool) {
  if (is.null(pool)) return(invisible(NULL))
  if (isTRUE(pool$alive)) {
    for (w in seq_len(pool$n)) {
      try({ serialize(NULL, pool$wcs[[w]]); flush(pool$wcs[[w]]) }, silent = TRUE)
    }
  }
  for (cn in c(pool$wcs, pool$rcs)) try(close(cn), silent = TRUE)
  try(parallel::mccollect(pool$jobs, wait = TRUE), silent = TRUE)
  unlink(pool$dir, recursive = TRUE)
  invisible(NULL)
}

# --- worker -----------------------------------------------------------------

.emc_wpool_serve <- function(req, ans, ctx) {
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
.emc_wpool_compute <- function(msg, ctx) .emc_with_preserved_rng({
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
      prev_ll = msg$prev_ll[k], parameters = msg$pars, model = ctx$model,
      stage = ctx$stage, type = ctx$type, tune = ctx$tune,
      marginalise = ctx$marginalise, r_cores = ctx$r_cores,
      chol_cache = ctx$chol_caches[[s]], group_chol = msg$group_chol
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
  n_subjects <- length(pm_settings)
  props <- matrix(0, ctx$n_pars + 1L, n_subjects)
  times <- numeric(n_subjects)

  msgs <- lapply(part, function(subs) {
    list(subs = subs, pars = pars, group_chol = group_chol,
         pm = pm_settings[subs], prev_ll = prev_ll[subs], seeds = seeds[subs])
  })
  was_alive <- isTRUE(pool$alive)
  sent <- rep(FALSE, pool$n)
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
      res <- .emc_wpool_compute(msgs[[w]], ctx)
    }
    props[, subs] <- res$props
    times[subs] <- res$times
    pm_settings[subs] <- res$pm
    seeds[subs] <- res$seeds
  }
  list(props = props, pm_settings = pm_settings, seeds = seeds, times = times,
       alive = pool$alive)
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
