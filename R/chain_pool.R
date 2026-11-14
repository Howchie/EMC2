
.emc_pool_state <- new.env(parent = emptyenv())

# A persistent chain keeps these large inputs from its initial fork.  They are
# immutable during a run_emc() stage; the rest of the sampler is the bounded
# continuation state exchanged at block boundaries.
.EMC_CHAIN_RETAINED_FIELDS <- c("data", "model", "prior")

.emc_chain_state_pack <- function(sampler) {
  fields <- setdiff(names(sampler), .EMC_CHAIN_RETAINED_FIELDS)
  list(version = 1L,
       names = names(sampler),
       attributes = attributes(sampler),
       fields = sampler[fields])
}

.emc_chain_state_apply <- function(base, state) {
  if (!is.list(state) || !identical(state$version, 1L) ||
      !is.character(state$names) || anyDuplicated(state$names) ||
      !is.list(state$fields) || is.null(names(state$fields)) ||
      !is.list(state$attributes) ||
      !all(setdiff(state$names, .EMC_CHAIN_RETAINED_FIELDS) %in%
             names(state$fields)) ||
      !all(intersect(state$names, .EMC_CHAIN_RETAINED_FIELDS) %in%
             names(base))) {
    stop("invalid persistent-chain state message", call. = FALSE)
  }
  out <- vector("list", length(state$names))
  names(out) <- state$names
  retained <- intersect(state$names, .EMC_CHAIN_RETAINED_FIELDS)
  dynamic <- setdiff(state$names, .EMC_CHAIN_RETAINED_FIELDS)
  out[retained] <- base[retained]
  out[dynamic] <- state$fields[dynamic]
  attributes(out) <- state$attributes
  out
}

.emc_chain_process_serve <- function(req, ans, sampler, runner, slot,
                                     other_slots) {
  if (identical(tryCatch(emc_arm_parent_death(), error = function(e) 0L), -1L))
    return(invisible(NULL))
  for (other in other_slots) {
    if (is.null(other$template) || identical(other, slot)) next
    for (cn in list(other$template$wc, other$template$rc))
      if (!is.null(cn)) try(close(cn), silent = TRUE)
  }
  .emc_pool_state$session_slot <- slot
  .emc_pool_state$session_in_use <- FALSE
  .emc_pool_state$chain_persistent <- TRUE
  .emc_pool_state$worker_stage_blocks <- 0L
  .emc_pool_state$subject_pool <- NULL
  .emc_pool_state$subject_ll_pool <- NULL
  on.exit({
    subject_pool <- .emc_pool_state$subject_pool
    subject_ll_pool <- .emc_pool_state$subject_ll_pool
    .emc_pool_state$subject_pool <- NULL
    .emc_pool_state$subject_ll_pool <- NULL
    .emc_pool_state$chain_persistent <- FALSE
    .emc_wpool_stop(subject_pool)
    .emc_wpool_stop(subject_ll_pool)
  }, add = TRUE)
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE) },
          add = TRUE)
  .emc_wpool_send(wc, list(ready = TRUE, pid = Sys.getpid()))
  repeat {
    msg <- .emc_wpool_recv(rc)
    if (is.null(msg)) break
    if (identical(msg$command, "stop")) {
      .emc_wpool_send(wc, list(id = msg$id, stopped = TRUE))
      break
    }
    if (identical(msg$command, "rebuild")) {
      result <- tryCatch(.emc_chain_process_rebuild(msg$tasks, sampler),
                         error = identity)
      reply <- if (inherits(result, "error")) {
        list(id = msg$id, error = result)
      } else list(id = msg$id, results = result)
      if (!is.null(msg$done)) try(file.create(msg$done), silent = TRUE)
      sent <- .emc_wpool_send(wc, reply)
      if (!sent || inherits(result, "error")) break
      next
    }
    warnings <- character()
    result <- tryCatch(withCallingHandlers({
      sampler <- .emc_chain_state_apply(sampler, msg$state)
      sampler <- do.call(runner, c(list(sampler = sampler), msg$args))
      list(state = .emc_chain_state_pack(sampler))
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
    reply <- if (inherits(result, "error")) {
      list(id = msg$id, error = result, warnings = warnings)
    } else {
      list(id = msg$id, state = result$state, warnings = warnings)
    }
    # The parent waits on this marker before reading the response. Publish it
    # before writing: an error condition can serialize to more than a FIFO's
    # capacity, otherwise leaving the child blocked in write while the parent
    # waits for the marker that follows the write.
    if (!is.null(msg$done)) {
      try(file.create(msg$done), silent = TRUE)
    }
    sent <- .emc_wpool_send(wc, reply)
    if (!sent || inherits(result, "error")) break
  }
  invisible(NULL)
}

.emc_chain_run_stage <- function(sampler, ...) {
  sampler <- .emc_set_sampler_rng(sampler)
  run_stages(sampler, ...)
}

.emc_chain_process_start <- function(sampler, runner = .emc_chain_run_stage,
                                     slot = NULL, other_slots = list()) {
  if (Sys.info()[["sysname"]] == "Windows" || !nzchar(Sys.which("mkfifo")))
    return(NULL)
  dir <- tempfile("emc_chain_")
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) return(NULL)
  req <- file.path(dir, "req")
  ans <- file.path(dir, "ans")
  made <- tryCatch(system2("mkfifo", shQuote(c(req, ans)),
                           stdout = FALSE, stderr = FALSE),
                   warning = function(e) 1L, error = function(e) 1L)
  if (!identical(as.integer(made), 0L)) {
    unlink(dir, recursive = TRUE)
    return(NULL)
  }
  wc <- rc <- NULL
  job <- NULL
  success <- FALSE
  on.exit({
    if (!success) {
      for (cn in list(wc, rc)) if (!is.null(cn)) try(close(cn), silent = TRUE)
      if (!is.null(job)) .emc_wpool_terminate_jobs(list(job), wait = TRUE)
      unlink(dir, recursive = TRUE)
    }
  }, add = TRUE)

  old_slot <- .emc_pool_state$session_slot
  .emc_pool_state$session_slot <- slot
  on.exit(.emc_pool_state$session_slot <- old_slot, add = TRUE)
  serve <- .emc_armed(function() .emc_chain_process_serve(
    req, ans, sampler, runner, slot, other_slots))
  job <- tryCatch(parallel::mcparallel(serve(), detached = FALSE),
                  error = function(e) NULL)
  .emc_pool_state$session_slot <- old_slot
  if (is.null(job)) return(NULL)

  f <- NULL
  deadline <- proc.time()[["elapsed"]] + 20
  wait <- 0.0002
  repeat {
    if (!.emc_wpool_pid_alive(job$pid)) break
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
  alive <- local({ j <- job; function() .emc_wpool_pid_alive(j$pid) })
  ready <- .emc_wpool_recv(rc, deadline = 20, alive = alive)
  if (!isTRUE(ready$ready) || !identical(as.integer(ready$pid), job$pid))
    return(NULL)
  success <- TRUE
  slot_pid <- if (is.null(slot$template$job)) NA_integer_ else
    .emc_wpool_job_pid(slot$template$job)
  list(dir = dir, job = job, wc = wc, rc = rc, next_id = 0L,
       slot_pid = slot_pid)
}

.emc_chain_slot_pid <- function(slot) {
  if (is.null(slot$template$job)) NA_integer_ else
    .emc_wpool_job_pid(slot$template$job)
}

.emc_chain_process_stop <- function(worker) {
  if (is.null(worker)) return(invisible(NULL))
  pid <- .emc_wpool_job_pid(worker$job)
  stopped <- FALSE
  if (!is.na(pid) && .emc_wpool_pid_alive(pid)) {
    id <- worker$next_id + 1L
    if (.emc_wpool_send(worker$wc, list(command = "stop", id = id))) {
      reply <- .emc_wpool_recv(worker$rc, deadline = 2,
                               alive = function() .emc_wpool_pid_alive(pid))
      stopped <- isTRUE(reply$stopped)
    }
  }
  for (cn in list(worker$wc, worker$rc)) try(close(cn), silent = TRUE)
  .emc_wpool_terminate_jobs(list(worker$job), wait = TRUE,
                            terminate = !stopped)
  unlink(worker$dir, recursive = TRUE)
  invisible(stopped)
}

.emc_chain_process_pool_start <- function(samplers, runner = .emc_chain_run_stage,
                                          slots = NULL) {
  if (is.null(slots)) slots <- vector("list", length(samplers))
  if (length(slots) != length(samplers))
    stop("one process slot is required per chain", call. = FALSE)
  if (!is.null(getOption("emc2.chain_process_memory_file", NULL)))
    .emc_pool_state$chain_process_blocks <- 0L
  workers <- vector("list", length(samplers))
  ok <- FALSE
  on.exit(if (!ok) for (w in workers) .emc_chain_process_stop(w), add = TRUE)
  for (i in seq_along(samplers)) {
    workers[[i]] <- .emc_chain_process_start(
      samplers[[i]], runner, slot = slots[[i]], other_slots = slots[-i])
    if (is.null(workers[[i]])) return(NULL)
  }
  ok <- TRUE
  pool <- new.env(parent = emptyenv())
  pool$workers <- workers
  pool$next_id <- 0L
  pool$blocks <- integer(length(workers))
  pool$runner <- runner
  pool$slots <- slots
  pool
}

.emc_chain_process_pool_stop <- function(pool) {
  if (is.null(pool)) return(invisible(NULL))
  stopped <- vapply(pool$workers, function(w)
    isTRUE(.emc_chain_process_stop(w)), logical(1))
  invisible(all(stopped))
}

.emc_chain_process_pool_refresh_slots <- function(pool, samplers, slots) {
  if (is.null(pool) || length(pool$workers) != length(samplers) ||
      length(slots) != length(samplers))
    stop("persistent chain slots do not match the active chains", call. = FALSE)
  for (i in seq_along(samplers)) {
    new_pid <- .emc_chain_slot_pid(slots[[i]])
    if (identical(pool$workers[[i]]$slot_pid, new_pid)) next
    .emc_chain_process_stop(pool$workers[[i]])
    replacement <- .emc_chain_process_start(
      samplers[[i]], pool$runner, slot = slots[[i]], other_slots = slots[-i])
    if (is.null(replacement))
      stop(sprintf("could not refresh persistent chain process %d", i),
           call. = FALSE)
    pool$workers[[i]] <- replacement
    pool$slots[[i]] <- slots[[i]]
    pool$blocks[[i]] <- 0L
  }
  invisible(pool)
}

.emc_chain_process_memory_record <- function(pool) {
  path <- getOption("emc2.chain_process_memory_file", NULL)
  if (!is.character(path) || length(path) != 1L || !nzchar(path) ||
      !dir.exists(dirname(path))) return(invisible(NULL))
  block <- if (is.null(.emc_pool_state$chain_process_blocks)) 1L else
    as.integer(.emc_pool_state$chain_process_blocks) + 1L
  .emc_pool_state$chain_process_blocks <- block
  rows <- lapply(seq_along(pool$workers), function(i) {
    pid <- .emc_wpool_job_pid(pool$workers[[i]]$job)
    file <- file.path("/proc", pid, "smaps_rollup")
    lines <- if (file.exists(file))
      tryCatch(suppressWarnings(readLines(file, warn = FALSE)),
               error = function(e) character()) else character()
    value <- function(key) {
      hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
      if (!length(hit)) return(NA_real_)
      suppressWarnings(as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", hit[1L])))
    }
    data.frame(block = block, chain = i, pid = pid,
               rss_kb = value("Rss"), pss_kb = value("Pss"),
               private_kb = value("Private_Clean") + value("Private_Dirty"))
  })
  write.table(do.call(rbind, rows), path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path))
  invisible(NULL)
}

.emc_chain_process_pool_run <- function(pool, samplers, args = list(),
                                        max_active = length(samplers)) {
  if (is.null(pool) || length(pool$workers) != length(samplers))
    stop("persistent chain process pool is unavailable", call. = FALSE)
  max_active <- max(1L, min(length(samplers), as.integer(max_active)))
  out <- vector("list", length(samplers))
  pending <- seq_along(samplers)
  active <- integer()
  ids <- integer(length(samplers))
  submit <- function(i) {
    worker <- pool$workers[[i]]
    pool$next_id <- pool$next_id + 1L
    id <- pool$next_id
    marker <- file.path(worker$dir, paste0("done_", id))
    unlink(marker)
    ok <- .emc_wpool_send(worker$wc, list(
      command = "run", id = id, done = marker,
      state = .emc_chain_state_pack(samplers[[i]]), args = args))
    if (!ok) stop(.emc_wpool_transport_error("chain request send failed"),
                   call. = FALSE)
    ids[i] <<- id
    active <<- c(active, i)
  }
  while (length(pending) && length(active) < max_active) {
    i <- pending[[1L]]; pending <- pending[-1L]; submit(i)
  }
  # The done markers are cheap to test; the liveness check reads /proc, so it
  # runs at most every 0.2 s, and the wait backs off to 5 ms between tests.
  alive_checked <- proc.time()[["elapsed"]]
  nap <- 0.001
  while (length(active)) {
    finished <- active[vapply(active, function(i) {
      file.exists(file.path(pool$workers[[i]]$dir,
                            paste0("done_", ids[i])))
    }, logical(1))]
    if (!length(finished)) {
      now <- proc.time()[["elapsed"]]
      if (now - alive_checked >= 0.2) {
        alive_checked <- now
        dead <- active[!vapply(active, function(i)
          .emc_wpool_pid_alive(pool$workers[[i]]$job$pid), logical(1))]
        if (length(dead))
          stop("persistent chain process exited without a reply", call. = FALSE)
      }
      Sys.sleep(nap)
      nap <- min(0.005, nap * 1.5)
      next
    }
    nap <- 0.001
    for (i in finished) {
      worker <- pool$workers[[i]]
      reply <- .emc_wpool_recv(
        worker$rc, deadline = .emc_wpool_deadline(),
        alive = local({ pid <- worker$job$pid
          function() .emc_wpool_pid_alive(pid) }))
      if (is.null(reply) || !identical(as.integer(reply$id), ids[i]))
        stop(.emc_wpool_transport_error("invalid persistent-chain reply"),
             call. = FALSE)
      for (msg in reply$warnings) warning(msg, call. = FALSE)
      if (inherits(reply$error, "error")) stop(reply$error)
      if (is.null(reply$state))
        stop("persistent chain returned no continuation state", call. = FALSE)
      out[[i]] <- .emc_chain_state_apply(samplers[[i]], reply$state)
      active <- active[active != i]
      if (length(pending)) {
        next_i <- pending[[1L]]; pending <- pending[-1L]; submit(next_i)
      }
    }
  }
  .emc_chain_process_memory_record(pool)
  # A chain process's heap only grows: it starts with a copy of the session's
  # and ratchets up with every block's largest temporaries. A fresh fork costs
  # about as much as the block boundary it replaces.
  restart_blocks <- getOption("emc2.persistent_chain_restart_blocks", 2L)
  if (identical(restart_blocks, FALSE) || is.null(restart_blocks)) {
    restart_blocks <- Inf
  } else {
    restart_blocks <- suppressWarnings(as.numeric(restart_blocks)[1L])
    if (!length(restart_blocks) || !is.finite(restart_blocks) ||
        restart_blocks < 1) restart_blocks <- Inf
  }
  if (is.finite(restart_blocks)) {
    for (i in seq_along(samplers)) {
      pool$blocks[i] <- pool$blocks[i] + 1L
      if (pool$blocks[i] < restart_blocks) next
      slot <- pool$slots[[i]]
      .emc_chain_process_stop(pool$workers[[i]])
      replacement <- .emc_chain_process_start(
        out[[i]], pool$runner, slot = slot, other_slots = pool$slots[-i])
      if (is.null(replacement))
        stop(sprintf("could not restart persistent chain process %d", i),
             call. = FALSE)
      pool$workers[[i]] <- replacement
      pool$blocks[i] <- 0L
    }
  }
  out
}

# Rebuild jobs run on the same workers retained by this chain's sampling
# process. A lost worker pool falls back to the owning chain process; each task
# restores its original subject RNG seed before it computes.
.emc_chain_process_rebuild <- function(tasks, sampler) {
  usable <- function(p) !is.null(p) && (isTRUE(p$alive) || isTRUE(p$parked))
  pool_name <- if (usable(.emc_pool_state$subject_pool)) {
    "subject_pool"
  } else "subject_ll_pool"
  pool <- .emc_pool_state[[pool_name]]
  if (isTRUE(pool$parked)) {
    unparked <- .emc_wpool_unpark(pool, pool$ctx, pool$n)
    if (!is.null(unparked)) {
      pool <- unparked
      .emc_pool_state[[pool_name]] <- pool
      on.exit(.emc_pool_state[[pool_name]] <-
                .emc_wpool_park(.emc_pool_state[[pool_name]]), add = TRUE)
    }
  }
  ctx <- if (!is.null(pool) && !is.null(pool$ctx)) pool$ctx else
    list(data = sampler$data, model = sampler$model,
         par_names = rownames(sampler$samples$alpha), n_pars = sampler$n_pars)
  ctx$data <- sampler$data
  run <- function(task) .emc_laplace_rebuild_compute(task, ctx)
  if (is.null(pool) || !isTRUE(pool$alive) || !length(pool$jobs))
    return(lapply(tasks, run))

  costs <- .emc_subject_cost(ctx$data)
  task_cost <- vapply(tasks, function(task) costs[[task$subject]], numeric(1))
  partition <- .emc_lpt_partition(task_cost, pool$n)
  messages <- lapply(partition, function(idx) list(
    kind = "laplace_rebuild_batch", tasks = tasks[idx]))
  sent <- vapply(seq_len(pool$n), function(w) {
    if (!length(partition[[w]])) return(TRUE)
    .emc_wpool_request(pool, w, messages[[w]])
  }, logical(1))
  replies <- vector("list", pool$n)
  if (all(sent)) replies <- lapply(seq_len(pool$n), function(w) {
    if (!length(partition[[w]])) return(list(results = list()))
    .emc_wpool_reply(pool, w)
  })
  valid <- all(vapply(seq_len(pool$n), function(w)
    !length(partition[[w]]) ||
      (is.list(replies[[w]]) && is.list(replies[[w]]$results) &&
       length(replies[[w]]$results) == length(partition[[w]]) &&
       all(vapply(replies[[w]]$results, function(x)
         is.list(x) && !is.null(x$seed), logical(1)))), logical(1)))
  if (!valid) {
    pool$alive <- FALSE
    .emc_pool_state[[pool_name]] <- pool
    detail <- if (!all(sent)) "request send failed" else {
      failed <- vapply(replies, function(reply) {
        if (!is.list(reply)) return("worker returned no reply")
        if (!is.null(reply$failed)) return(as.character(reply$failed)[1L])
        if (!is.list(reply$results)) return("worker reply has no results")
        "worker result omitted its subject seed"
      }, character(1))
      paste(unique(failed[nzchar(failed)]), collapse = "; ")
    }
    .emc_wpool_degraded(paste0(
      "retained rebuild worker failed (", detail,
      "); retrying in its chain"))
    return(lapply(tasks, run))
  }
  out <- vector("list", length(tasks))
  for (w in seq_len(pool$n)) if (length(partition[[w]]))
    out[partition[[w]]] <- replies[[w]]$results
  out
}

# Fan one chain's rebuild tasks out over its retained worker pool. Rebuild
# tasks are independent and carry explicit subject streams, so an LPT batch per
# worker keeps scheduling overhead low while balancing uneven subjects.
.emc_chain_process_pool_rebuild <- function(pool, samplers, tasks_by_chain,
                                           contexts,
                                           max_active = length(samplers)) {
  if (is.null(pool) || length(pool$workers) != length(samplers) ||
      length(tasks_by_chain) != length(samplers) ||
      length(contexts) != length(samplers)) return(NULL)
  max_active <- max(1L, min(length(samplers), as.integer(max_active)))
  out <- vector("list", length(samplers))
  pending <- which(lengths(tasks_by_chain) > 0L)
  active <- integer()
  ids <- integer(length(samplers))
  recover <- function(i) {
    worker <- pool$workers[[i]]
    for (cn in list(worker$wc, worker$rc)) try(close(cn), silent = TRUE)
    .emc_wpool_terminate_jobs(list(worker$job), wait = TRUE, terminate = TRUE)
    unlink(worker$dir, recursive = TRUE)
    replacement <- .emc_chain_process_start(
      samplers[[i]], pool$runner, slot = pool$slots[[i]],
      other_slots = pool$slots[-i])
    if (is.null(replacement))
      stop(sprintf("could not restart persistent chain process %d after rebuild failure", i),
           call. = FALSE)
    pool$workers[[i]] <- replacement
    lapply(tasks_by_chain[[i]], function(task)
      .emc_laplace_rebuild_compute(task, contexts[[i]]))
  }
  submit <- function(i) {
    worker <- pool$workers[[i]]
    pool$next_id <- pool$next_id + 1L
    id <- pool$next_id
    marker <- file.path(worker$dir, paste0("done_", id))
    unlink(marker)
    ok <- .emc_wpool_send(worker$wc, list(
      command = "rebuild", id = id, done = marker,
      tasks = tasks_by_chain[[i]]))
    if (!ok) {
      out[[i]] <<- recover(i)
      return(FALSE)
    }
    ids[i] <<- id
    active <<- c(active, i)
    TRUE
  }
  while (length(pending) && length(active) < max_active) {
    i <- pending[[1L]]; pending <- pending[-1L]; submit(i)
  }
  alive_checked <- proc.time()[["elapsed"]]
  nap <- 0.001
  while (length(active)) {
    finished <- active[vapply(active, function(i) {
      file.exists(file.path(pool$workers[[i]]$dir,
                            paste0("done_", ids[i])))
    }, logical(1))]
    if (!length(finished)) {
      now <- proc.time()[["elapsed"]]
      dead <- integer()
      if (now - alive_checked >= 0.2) {
        alive_checked <- now
        dead <- active[!vapply(active, function(i)
          .emc_wpool_pid_alive(pool$workers[[i]]$job$pid), logical(1))]
      }
      if (length(dead)) finished <- dead else {
        Sys.sleep(nap)
        nap <- min(0.005, nap * 1.5)
        next
      }
    }
    nap <- 0.001
    for (i in finished) {
      worker <- pool$workers[[i]]
      reply <- tryCatch(.emc_wpool_recv(
        worker$rc, deadline = .emc_wpool_deadline(),
        alive = local({ pid <- worker$job$pid
          function() .emc_wpool_pid_alive(pid) })), error = function(e) NULL)
      if (is.list(reply) && identical(as.integer(reply$id), ids[i]) &&
          is.list(reply$results) &&
          length(reply$results) == length(tasks_by_chain[[i]])) {
        out[[i]] <- reply$results
      } else {
        out[[i]] <- recover(i)
      }
      active <- active[active != i]
      if (length(pending)) {
        next_i <- pending[[1L]]; pending <- pending[-1L]; submit(next_i)
      }
    }
  }
  out
}

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

.emc_wpool_request <- function(pool, w, msg, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_send(pool$wcs[[w]], msg, deadline = deadline)
}

.emc_wpool_reply <- function(pool, w, deadline = NULL) {
  if (is.null(deadline)) deadline <- .emc_wpool_deadline()
  .emc_wpool_recv(pool$rcs[[w]], deadline = deadline,
                  alive = .emc_wpool_alive_fn(pool, w),
                  busy = .emc_wpool_busy_fn(pool, w))
}

.emc_wpool_busy_fn <- function(pool, w) {
  job <- if (!is.null(pool$jobs) && w <= length(pool$jobs)) pool$jobs[[w]] else NULL
  if (is.null(job)) return(NULL)
  function() .emc_wpool_busy(.emc_wpool_job_pid(job))
}

# CPU ticks (user + system) a process has used, from /proc; NA where /proc is
# unavailable or the process is gone.
.emc_wpool_cpu_ticks <- function(pid) {
  if (length(pid) != 1L || is.na(pid) || !dir.exists("/proc")) return(NA_real_)
  stat <- suppressWarnings(
    tryCatch(readLines(sprintf("/proc/%d/stat", as.integer(pid)), n = 1L),
             error = function(e) character()))
  if (!length(stat)) return(NA_real_)
  # Fields after "(comm) " start at the state (field 3): utime and stime are
  # fields 14 and 15.
  f <- strsplit(sub("^.*\\) ", "", stat), " ", fixed = TRUE)[[1L]]
  if (length(f) < 13L) return(NA_real_)
  sum(suppressWarnings(as.numeric(f[12:13])))
}

# A worker past the reply deadline whose CPU time still advances over a short
# probe is computing (slow, not lost); one that is alive but idle is stuck.
.emc_wpool_busy <- function(pids, probe = 1) {
  pids <- pids[!is.na(pids)]
  if (!length(pids)) return(FALSE)
  t0 <- vapply(pids, .emc_wpool_cpu_ticks, numeric(1))
  if (!any(is.finite(t0))) return(FALSE)
  Sys.sleep(probe)
  t1 <- vapply(pids, .emc_wpool_cpu_ticks, numeric(1))
  busy <- any(is.finite(t0) & is.finite(t1) & t1 > t0)
  if (busy) {
    .emc_pool_state$slow_waits <-
      (if (is.null(.emc_pool_state$slow_waits)) 0L
       else .emc_pool_state$slow_waits) + 1L
  }
  busy
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
  tryCatch({
    sent <- emc_wpool_send_cpp(con, payload, deadline, .EMC_WPOOL_MAX_SEGMENT,
                               .EMC_WPOOL_WRITE_CHUNK)
    if (!isTRUE(sent)) stop(sent, call. = FALSE)
    flush(con)
    TRUE
  }, error = function(e) .emc_wpool_transport_fail(conditionMessage(e)))
}

.emc_wpool_recv <- function(con, deadline = NULL, alive = NULL, busy = NULL) {
  payload <- tryCatch(
    emc_wpool_recv_cpp(con, deadline, alive, busy, .EMC_WPOOL_READ_CHUNK,
                       .EMC_WPOOL_MAX_FRAME),
    error = function(e) conditionMessage(e))
  if (!is.raw(payload)) {
    .emc_wpool_transport_fail(payload)
    return(NULL)
  }
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
  payload <- readBin(path, "raw", as.integer(size))
  if (length(payload) != as.integer(size)) {
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
        pids <- vapply(pending(), function(w) {
          pid <- .emc_wpool_job_pid(pool$jobs[[w]])
          if (is.null(pid)) NA_integer_ else as.integer(pid)
        }, integer(1))
        if (!.emc_wpool_busy(pids)) return(list(code = "timeout"))
        started <- proc.time()[["elapsed"]]
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
  if (!installed || !file.exists(rscript) ||
      (.emc_wpool_has_custom_ptr(ctx) && !.emc_wpool_custom_restorable(ctx$model))) {
    return("fork")
  }
  "spawn"
}

# A spawned template rebuilds custom trend kernels from their registration
# (.emc_wpool_template()), so only a kernel without one forces the fork backend.
.emc_wpool_custom_restorable <- function(model) {
  spec <- if (is.function(model)) tryCatch(model(), error = function(e) NULL)
          else model
  if (!is.list(spec)) return(TRUE)
  if (!is.null(spec$trend)) {
    return(all(vapply(spec$trend, function(entry) {
      if (!identical(entry$kernel, "custom")) return(TRUE)
      reg <- attr(entry, "custom_registration")
      !is.null(reg) && !is.null(reg$file) && file.exists(reg$file)
    }, logical(1))))
  }
  if (any(c("p_types", "rfun", "Ttransform") %in% names(spec))) return(TRUE)
  all(vapply(spec, .emc_wpool_custom_restorable, logical(1)))
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
  fields <- strsplit(sub("^.*\\) ", "", stat), " ", fixed = TRUE)[[1L]]
  if (fields[1L] %in% c("Z", "X", "x")) return(FALSE)
  # A process past exit() still runs while the kernel tears it down (PF_EXITING
  # in its flags); it will serve nothing more.
  flags <- suppressWarnings(as.numeric(fields[7L]))
  !(is.finite(flags) && bitwAnd(as.integer(flags %% 8), 4L) != 0L)
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


.emc_wpool_spawn <- function(dir, idx, ctx, backend = "fork", template = NULL,
                             ctx_update = NULL, context_epoch = 0L) {
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
      id <- .emc_wpool_next_request()
      if (!.emc_wpool_send(template$wc, list(command = "spawn", idx = idx,
                                             req = req, ans = ans,
                                             ctx = ctx_update, id = id,
                                             context_epoch = context_epoch))) {
        stop(.emc_wpool_transport_error("template stopped accepting requests"))
      }
      # A template kept across blocks may still hold the reply to a command
      # whose sender died before reading it.
      repeat {
        reply <- .emc_wpool_recv(
          template$rc, deadline = .EMC_WPOOL_TIMEOUT_TEMPLATE,
          alive = local({ j <- template$job; function() {
            pid <- .emc_wpool_job_pid(j)
            is.null(pid) || is.na(pid) || .emc_wpool_pid_alive(pid)
          } }))
        if (is.null(reply) || identical(reply$id, id)) break
      }
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
          parallel::mcparallel(.emc_wpool_serve(
            req[i], ans[i], ctx, context_epoch = context_epoch),
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

# One Rscript serving several template slots. It loads the package and the
# base context once and forks a server per slot, each running the template loop
# on that slot's FIFOs; the slots then share the booted heap instead of each
# holding a copy. The root exits once every slot has. The slots are detached
# and watched by pid: collecting them with mccollect() returns while they are
# still serving, and the root's exit cleanup then kills them.
.emc_wpool_template_root <- function(boot) {
  base <- readRDS(boot$ctx)
  base$model <- .restore_custom_model_container(base$model, quiet = TRUE)
  pids <- vapply(boot$slots, function(slot) parallel::mcparallel(
    .emc_wpool_template(slot, base = base), detached = TRUE)$pid, integer(1))
  while (any(vapply(pids, .emc_wpool_pid_alive, logical(1)))) Sys.sleep(0.25)
  invisible(NULL)
}

# Launch the slots of a session from one root (.emc_wpool_template_root()).
# Returns one launch record per directory, each connected like a single
# template's; NULL entries for slots whose FIFOs could not be made.
.emc_wpool_template_launch_shared <- function(dirs, ctx_file, lib) {
  slots <- lapply(dirs, function(dir) {
    req <- file.path(dir, "template_req")
    ans <- file.path(dir, "template_ans")
    made <- tryCatch(
      system2("mkfifo", shQuote(c(req, ans)), stdout = FALSE, stderr = FALSE),
      warning = function(e) 1L, error = function(e) 1L
    )
    if (!identical(as.integer(made), 0L)) return(NULL)
    list(req = req, ans = ans, pidfile = file.path(dir, "template.pid"))
  })
  ok <- !vapply(slots, is.null, logical(1))
  if (!any(ok)) return(vector("list", length(dirs)))
  root_dir <- dirs[which(ok)[1L]]
  boot <- file.path(root_dir, "template_root.rds")
  log <- file.path(root_dir, "template.log")
  saveRDS(list(slots = slots[ok], ctx = ctx_file, lib = lib), boot,
          compress = FALSE)
  tag <- sprintf("EMC2-worker-pool[%s|master=%d]",
                 gsub("[^A-Za-z0-9._-]", "_", basename(root_dir)), Sys.getpid())
  expr <- paste0(
    "EMC2.worker<-'", tag, "';",
    "b<-readRDS(commandArgs(TRUE)[1L]);",
    ".libPaths(unique(c(b$lib,.libPaths())));",
    "ns<-loadNamespace('EMC2');",
    "get('.emc_wpool_template_root',ns)(b)"
  )
  status <- tryCatch(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(expr), shQuote(boot)),
    stdout = log, stderr = log, wait = FALSE,
    env = paste0("EMC2_CPP_CACHE=", shQuote(.emc_cpp_cache_dir()))
  ), error = function(e) 1L)
  if (!identical(as.integer(status), 0L)) {
    for (slot in slots[ok]) unlink(c(slot$req, slot$ans))
    unlink(c(boot, log))
    return(vector("list", length(dirs)))
  }
  lapply(slots, function(slot) if (is.null(slot)) NULL else
    list(req = slot$req, ans = slot$ans, boot = boot, pidfile = slot$pidfile,
         log = log))
}

.emc_wpool_template_start <- function(dir, ctx_file, lib) {
  .emc_wpool_template_connect(.emc_wpool_template_launch(dir, ctx_file, lib))
}

# Start the template's Rscript without waiting for it to load, so several can
# boot at once; .emc_wpool_template_connect() then waits for each.
.emc_wpool_template_launch <- function(dir, ctx_file, lib) {
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
    stdout = log, stderr = log, wait = FALSE,
    env = paste0("EMC2_CPP_CACHE=", shQuote(.emc_cpp_cache_dir()))
  ), error = function(e) 1L)
  if (!identical(as.integer(status), 0L)) {
    unlink(c(req, ans, boot, pidfile, log))
    return(NULL)
  }
  list(req = req, ans = ans, boot = boot, pidfile = pidfile, log = log)
}

.emc_wpool_template_connect <- function(launch) {
  if (is.null(launch)) return(NULL)
  req <- launch$req
  ans <- launch$ans
  boot <- launch$boot
  pidfile <- launch$pidfile
  log <- launch$log
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

# `persistent = TRUE` forks the workers from this chain's template kept by
# run_emc() across blocks (.emc_wpool_session_start()) when one is available,
# instead of booting a template for this pool alone.
# `retain_across_blocks = TRUE` instead keeps fork workers owned by one
# persistent chain process, refreshing their mutable context between blocks;
# they too are forked from the chain's session template when it has one.
.emc_wpool_start <- function(n_workers, ctx, allow_one = FALSE,
                             persistent = FALSE, backend = NULL,
                             retain_across_blocks = FALSE) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L ||
      (!isTRUE(allow_one) && n_workers <= 1L)) return(NULL)
  if (Sys.info()[["sysname"]] == "Windows") return(NULL)
  if (!nzchar(Sys.which("mkfifo"))) return(NULL)
  n_workers <- .emc_wpool_feasible_workers(n_workers)
  if (n_workers < 1L || (!isTRUE(allow_one) && n_workers <= 1L)) return(NULL)

  dir <- tempfile("emc_wpool_")
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) return(NULL)

  if ((isTRUE(persistent) || isTRUE(retain_across_blocks)) && is.null(backend)) {
    pool <- .emc_wpool_start_kept(dir, n_workers, ctx)
    if (!is.null(pool)) {
      if (isTRUE(retain_across_blocks)) {
        pool$ctx <- ctx
        pool$retained <- TRUE
        pool$iterations <- 0L
        pool$base_n <- n_workers
        pool$context_epoch <- 0L
      }
      return(pool)
    }
  }

  if (is.null(backend)) backend <- .emc_wpool_backend(ctx)
  backend <- match.arg(backend, c("spawn", "fork"))
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
       template_data = if (is.null(template)) NULL else ctx$data,
       done = .emc_wpool_done_open(dir), ctx = ctx,
       retained = isTRUE(retain_across_blocks), iterations = 0L,
       base_n = n_workers, context_epoch = 0L)
}

.emc_wpool_retained_idle <- function(ctx, iterations = 0L, base_n = 1L) {
  base_n <- max(1L, as.integer(base_n))
  list(n = base_n, dir = NULL, jobs = list(), wcs = list(), rcs = list(),
       alive = FALSE, backend = "fork", template = NULL, done = NULL,
       ctx = ctx, retained = TRUE, iterations = as.integer(iterations),
       base_n = as.integer(base_n), context_epoch = 0L,
       generation = .emc_wpool_next_generation())
}

.emc_wpool_trace_retained <- function(pool, role, block, stage) {
  path <- getOption("emc2.worker_trace_file", NULL)
  if (!is.character(path) || length(path) != 1L || !nzchar(path) ||
      is.null(pool) || !isTRUE(pool$retained)) return(invisible(NULL))
  pids <- if (length(pool$jobs)) {
    vapply(pool$jobs, .emc_wpool_job_pid, integer(1))
  } else integer(0)
  memory <- lapply(pids, function(pid) {
    lines <- tryCatch(readLines(sprintf("/proc/%d/smaps_rollup", pid),
                                warn = FALSE), error = function(e) character())
    value <- function(key) {
      hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
      if (!length(hit)) return(NA_real_)
      suppressWarnings(as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", hit[1L])))
    }
    c(rss = value("Rss"), pss = value("Pss"),
      private = value("Private_Clean") + value("Private_Dirty"))
  })
  mem <- if (length(memory)) {
    do.call(rbind, memory)
  } else matrix(numeric(), nrow = 0L, ncol = 3L,
                dimnames = list(NULL, c("rss", "pss", "private")))
  mem_sum <- function(name) {
    x <- mem[, name]
    if (!length(x) || all(!is.finite(x))) NA_real_ else sum(x, na.rm = TRUE)
  }
  file <- sprintf("%s-%d.tsv", path, Sys.getpid())
  row <- data.frame(owner = Sys.getpid(), block = as.integer(block),
                    stage = as.character(stage), role = as.character(role),
                    workers = as.integer(pool$n), alive = isTRUE(pool$alive),
                    epoch = as.integer(pool$context_epoch),
                    worker_rss_kb = mem_sum("rss"),
                    worker_pss_kb = mem_sum("pss"),
                    worker_private_kb = mem_sum("private"),
                    pids = paste(pids, collapse = ";"),
                    stringsAsFactors = FALSE)
  utils::write.table(row, file = file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = !file.exists(file),
                     append = file.exists(file))
  invisible(NULL)
}

.emc_wpool_retained_static_equal <- function(a, b) {
  strip_pooled <- function(x) {
    if (is.list(x)) for (i in seq_along(x)) x[[i]] <- strip_pooled(x[[i]])
    attr(x, "pooled_values") <- NULL
    x
  }
  identical(strip_pooled(a$data), strip_pooled(b$data)) &&
    isTRUE(all.equal(a$model, b$model, tolerance = 0)) &&
    identical(a$type, b$type) && identical(a$n_pars, b$n_pars)
}

.emc_wpool_context_update <- function(pool, ctx) {
  if (is.null(pool) || !isTRUE(pool$retained) ||
      !isTRUE(pool$alive) || !length(pool$jobs) ||
      !.emc_wpool_retained_static_equal(pool$ctx, ctx)) return(NULL)
  values <- ctx[setdiff(names(ctx), c("data", "model"))]
  epoch <- as.integer(pool$context_epoch) + 1L
  sent <- vapply(seq_len(pool$n), function(w) {
    .emc_wpool_send(pool$wcs[[w]],
                    list(control = "context_update", epoch = epoch,
                         values = values))
  }, logical(1))
  if (!all(sent)) return(NULL)
  replies <- lapply(seq_len(pool$n), function(w)
    .emc_wpool_reply(pool, w))
  ok <- vapply(replies, function(x) is.list(x) && isTRUE(x$ok) &&
                 identical(as.integer(x$epoch), epoch), logical(1))
  if (!all(ok)) return(NULL)
  pool$ctx <- ctx
  pool$context_epoch <- epoch
  pool
}

.emc_wpool_shrink <- function(pool, n_total) {
  if (is.null(pool) || !isTRUE(pool$retained)) return(pool)
  n_total <- max(1L, as.integer(n_total))
  if (n_total >= pool$n) return(pool)
  iterations <- pool$iterations
  base_n <- pool$base_n
  ctx <- pool$ctx
  if (n_total <= 1L) {
    .emc_wpool_stop(pool)
    return(.emc_wpool_retained_idle(ctx, iterations, base_n))
  }
  if (!isTRUE(pool$alive) || is.null(pool$dir)) {
    .emc_wpool_stop(pool)
    return(.emc_wpool_retained_idle(ctx, iterations, n_total))
  }
  drop_idx <- seq.int(n_total + 1L, pool$n)
  dropped <- list(n = length(drop_idx), dir = pool$dir,
                  jobs = pool$jobs[drop_idx], wcs = pool$wcs[drop_idx],
                  rcs = pool$rcs[drop_idx], backend = pool$backend)
  .emc_wpool_stop_workers(dropped)
  unlink(file.path(pool$dir,
                   c(sprintf("req%d", drop_idx), sprintf("ans%d", drop_idx))))
  pool$n <- n_total
  pool$jobs <- pool$jobs[seq_len(n_total)]
  pool$wcs <- pool$wcs[seq_len(n_total)]
  pool$rcs <- pool$rcs[seq_len(n_total)]
  pool$generation <- .emc_wpool_next_generation()
  .emc_wpool_retire(pool)
  pool
}

.emc_wpool_reconfigure <- function(pool, ctx, n_workers) {
  n_workers <- max(1L, as.integer(n_workers))
  iterations <- if (is.null(pool$iterations)) 0L else pool$iterations
  compatible <- !is.null(pool) && isTRUE(pool$retained) &&
    .emc_wpool_retained_static_equal(pool$ctx, ctx)
  if (n_workers <= 1L) {
    if (compatible) {
      pool$base_n <- 1L
      return(.emc_wpool_shrink(pool, 1L))
    }
    if (!is.null(pool)) .emc_wpool_stop(pool)
    return(.emc_wpool_retained_idle(ctx, iterations, 1L))
  }
  if (compatible && isTRUE(pool$parked)) {
    unparked <- .emc_wpool_unpark(pool, ctx, n_workers)
    if (!is.null(unparked)) {
      unparked$iterations <- iterations
      unparked$base_n <- n_workers
      return(unparked)
    }
  }
  if (!compatible || !isTRUE(pool$alive)) {
    if (!is.null(pool)) .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(n_workers, ctx, retain_across_blocks = TRUE)
    if (is.null(fresh)) return(.emc_wpool_retained_idle(ctx, iterations,
                                                       n_workers))
    fresh$iterations <- iterations
    fresh$base_n <- n_workers
    return(fresh)
  }
  refreshed <- .emc_wpool_context_update(pool, ctx)
  if (is.null(refreshed)) {
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(n_workers, ctx, retain_across_blocks = TRUE)
    if (is.null(fresh)) return(.emc_wpool_retained_idle(ctx, iterations,
                                                       n_workers))
    fresh$iterations <- iterations
    fresh$base_n <- n_workers
    return(fresh)
  }
  pool <- refreshed
  pool$base_n <- n_workers
  if (pool$n > n_workers) pool <- .emc_wpool_shrink(pool, n_workers)
  if (pool$n < n_workers) pool <- .emc_wpool_grow(pool, n_workers, ctx)
  pool
}

.emc_wpool_reclaim_persistent <- function() {
  if (!isTRUE(.emc_pool_state$chain_persistent)) return(invisible(NULL))
  for (name in c("subject_pool", "subject_ll_pool")) {
    pool <- .emc_pool_state[[name]]
    if (is.null(pool) || !isTRUE(pool$retained)) next
    .emc_pool_state[[name]] <- .emc_wpool_park(
      .emc_wpool_shrink(pool, pool$base_n))
  }
  invisible(NULL)
}

.emc_wpool_template_live <- function(pool) {
  identical(pool$backend, "spawn") && !is.null(pool$template) &&
    isTRUE(pool$template$alive) &&
    .emc_wpool_pid_alive(.emc_wpool_job_pid(pool$template$job))
}

# Between blocks a persistent chain's workers would sit idle while the session
# forks its own jobs, so a pool with a live template stops its workers and
# keeps the template; the workers are forked again, with the then-current
# context, when the pool is next needed. Workers are re-forked from the
# template every few dozen iterations anyway, so this adds little.
.emc_wpool_park <- function(pool) {
  if (is.null(pool) || !isTRUE(pool$retained) || !isTRUE(pool$alive) ||
      !length(pool$jobs) || !.emc_wpool_template_live(pool)) return(pool)
  .emc_wpool_stop_workers(pool)
  unlink(file.path(pool$dir, c(sprintf("req%d", seq_len(pool$n)),
                               sprintf("ans%d", seq_len(pool$n)))))
  pool$jobs <- list()
  pool$wcs <- list()
  pool$rcs <- list()
  pool$alive <- FALSE
  pool$parked <- TRUE
  pool$generation <- .emc_wpool_next_generation()
  .emc_wpool_retire(pool)
  pool
}

# Writes ctx as an update over the context a template booted with.  The data
# are left out when the template already holds the same data: a template
# reads every update into its own heap, so resending them grows it.
.emc_wpool_save_update <- function(ctx, template_data, path) {
  upd <- ctx
  if (!is.null(template_data) && identical(ctx$data, template_data))
    upd$data <- NULL
  tryCatch({ saveRDS(upd, path, compress = FALSE); TRUE },
           error = function(e) FALSE)
}

# NULL when the template can no longer fork the workers.
.emc_wpool_unpark <- function(pool, ctx, n_workers) {
  if (!isTRUE(pool$parked) || !.emc_wpool_template_live(pool)) return(NULL)
  n_workers <- max(1L, as.integer(n_workers))
  ctx_update <- file.path(pool$dir, "context_current.rds")
  if (!.emc_wpool_save_update(ctx, pool$template_data, ctx_update)) return(NULL)
  spawned <- .emc_wpool_spawn(pool$dir, seq_len(n_workers), ctx,
                              backend = "spawn", template = pool$template,
                              ctx_update = ctx_update,
                              context_epoch = pool$context_epoch)
  unlink(ctx_update)
  if (is.null(spawned)) return(NULL)
  pool$n <- n_workers
  pool$jobs <- spawned$jobs
  pool$wcs <- spawned$wcs
  pool$rcs <- spawned$rcs
  pool$ctx <- ctx
  pool$alive <- TRUE
  pool$parked <- NULL
  pool$generation <- .emc_wpool_next_generation()
  .emc_wpool_retire(pool)
  pool
}

# Fork a pool from the session template for this chain, sending it this
# block's context (without the data when the template already holds it).
.emc_wpool_start_kept <- function(dir, n_workers, ctx) {
  slot <- .emc_pool_state$session_slot
  template <- slot$template
  if (is.null(template) || !isTRUE(template$alive) ||
      isTRUE(.emc_pool_state$session_in_use) ||
      !.emc_wpool_pid_alive(.emc_wpool_job_pid(template$job)) ||
      !identical(.emc_wpool_backend(ctx), "spawn")) return(NULL)
  ctx_update <- file.path(dir, "context.rds")
  if (!.emc_wpool_save_update(ctx, slot$data, ctx_update)) return(NULL)
  .emc_pool_state$warned <- FALSE
  spawned <- .emc_wpool_spawn(dir, seq_len(n_workers), ctx, "spawn",
                              template = template, ctx_update = ctx_update)
  if (is.null(spawned)) {
    # The template is no use to this chain any more; later pools boot their own.
    .emc_pool_state$session_slot <- NULL
    unlink(ctx_update)
    return(NULL)
  }
  .emc_pool_state$session_in_use <- TRUE
  list(n = n_workers, dir = dir, jobs = spawned$jobs, wcs = spawned$wcs,
       rcs = spawned$rcs, alive = TRUE, backend = "spawn",
       generation = .emc_wpool_next_generation(),
       ctx_file = NULL, lib = NULL, template = template, kept = TRUE,
       template_data = slot$data,
       done = .emc_wpool_done_open(dir))
}

# Templates kept for the length of a run_emc() call, one per chain slot.
# Booting a template (Rscript, the package, the data) is most of the cost of
# starting a pool, and run_emc() forks fresh chain processes every block.  The
# workers are still forked afresh every block and recycled within one, so
# their memory does not accumulate; only the template, which runs no
# likelihoods, lives on.  A template exits when its control FIFO closes, so it
# cannot outlive this process.
.emc_wpool_session_start <- function(base_ctx, n_slots) {
  n_slots <- as.integer(n_slots)
  if (is.na(n_slots) || n_slots < 1L) return(NULL)
  if (Sys.info()[["sysname"]] == "Windows") return(NULL)
  if (!nzchar(Sys.which("mkfifo"))) return(NULL)
  if (!identical(.emc_wpool_backend(base_ctx), "spawn")) return(NULL)
  ns_path <- tryCatch(getNamespaceInfo(asNamespace("EMC2"), "path"),
                      error = function(e) "")
  if (!nzchar(ns_path)) return(NULL)
  lib <- dirname(ns_path)
  dirs <- character(0)
  launches <- vector("list", n_slots)
  slot_dirs <- rep(NA_character_, n_slots)
  for (k in seq_len(n_slots)) {
    dir <- tempfile("emc_wtemplate_")
    if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) next
    dirs <- c(dirs, dir)
    slot_dirs[k] <- dir
  }
  made <- !is.na(slot_dirs)
  if (any(made)) {
    ctx_file <- file.path(slot_dirs[which(made)[1L]], "context.rds")
    saved <- tryCatch({ saveRDS(base_ctx, ctx_file, compress = FALSE); TRUE },
                      error = function(e) FALSE)
    if (saved) launches[made] <- .emc_wpool_template_launch_shared(
      slot_dirs[made], ctx_file, lib)
  }
  templates <- lapply(launches, function(l) {
    if (is.null(l)) return(NULL)
    t <- .emc_wpool_template_connect(l)
    if (!is.null(t)) t$dir <- dirname(l$req)
    t
  })
  if (all(vapply(templates, is.null, logical(1)))) {
    unlink(dirs, recursive = TRUE)
    return(NULL)
  }
  list(templates = templates, data = base_ctx$data, base_ctx = base_ctx,
       lib = lib)
}

# Replace any template that has died since the last block.
.emc_wpool_session_check <- function(session) {
  if (is.null(session)) return(session)
  for (k in seq_along(session$templates)) {
    t <- session$templates[[k]]
    if (!is.null(t) && .emc_wpool_pid_alive(.emc_wpool_job_pid(t$job))) next
    if (!is.null(t)) {
      .emc_wpool_template_stop(t)
      unlink(t$dir, recursive = TRUE)
    }
    dir <- tempfile("emc_wtemplate_")
    session$templates[k] <- list(NULL)
    if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) next
    ctx_file <- file.path(dir, "context.rds")
    saved <- tryCatch({ saveRDS(session$base_ctx, ctx_file, compress = FALSE); TRUE },
                      error = function(e) FALSE)
    fresh <- if (saved) .emc_wpool_template_start(dir, ctx_file, session$lib) else NULL
    if (is.null(fresh)) {
      unlink(dir, recursive = TRUE)
      next
    }
    fresh$dir <- dir
    session$templates[[k]] <- fresh
  }
  session
}

.emc_wpool_session_stop <- function(session) {
  if (is.null(session)) return(invisible(NULL))
  for (t in session$templates) {
    if (is.null(t)) next
    .emc_wpool_template_stop(t)
    unlink(t$dir, recursive = TRUE)
  }
  invisible(NULL)
}

# What a chain process needs to use its slot's template.
.emc_wpool_session_slot <- function(session, k) {
  if (is.null(session) || k > length(session$templates) ||
      is.null(session$templates[[k]])) return(NULL)
  list(template = session$templates[[k]], data = session$data)
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
  if (is.null(pool) || is.na(n_total)) return(pool)
  if (isTRUE(pool$retained) && n_total > pool$n &&
      (!isTRUE(pool$alive) || is.null(pool$dir))) {
    iterations <- pool$iterations
    base_n <- pool$base_n
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(n_total, ctx, retain_across_blocks = TRUE)
    if (is.null(fresh)) return(.emc_wpool_retained_idle(ctx, iterations,
                                                       base_n))
    fresh$iterations <- iterations
    fresh$base_n <- base_n
    return(fresh)
  }
  if (!isTRUE(pool$alive)) return(pool)
  n_total <- min(n_total, pool$n + .emc_wpool_feasible_workers(n_total - pool$n))
  n_new <- n_total - pool$n
  if (n_new <= 0L) return(pool)

  if (isTRUE(pool$grow_skip > 0L)) {
    pool$grow_skip <- pool$grow_skip - 1L
    return(pool)
  }

  ctx_update <- NULL
  if (isTRUE(pool$retained) && identical(pool$backend, "spawn") &&
      !is.null(pool$template)) {
    ctx_update <- file.path(pool$dir, "context_current.rds")
    if (!.emc_wpool_save_update(ctx, pool$template_data, ctx_update))
      ctx_update <- NULL
  }
  spawned <- .emc_wpool_spawn(pool$dir, pool$n + seq_len(n_new), ctx,
                              backend = if (is.null(pool$backend)) "fork" else pool$backend,
                              template = pool$template, ctx_update = ctx_update,
                              context_epoch = pool$context_epoch)
  if (!is.null(ctx_update)) unlink(ctx_update)
  if (is.null(spawned)) {
    fails <- if (is.null(pool$grow_fails)) 1L else pool$grow_fails + 1L
    pool$grow_fails <- fails
    pool$grow_skip <- bitwShiftL(1L, min(fails, 6L))
    return(pool)
  }

  .emc_wpool_retire(pool)
  pool$n <- n_total
  pool$jobs <- c(pool$jobs, spawned$jobs)
  pool$wcs <- c(pool$wcs, spawned$wcs)
  pool$rcs <- c(pool$rcs, spawned$rcs)
  pool$ctx <- ctx
  pool$alive <- TRUE
  pool$grow_fails <- 0L
  pool$grow_skip <- 0L
  pool$generation <- .emc_wpool_next_generation()
  pool
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

  # A retained pool with a live template re-forks its workers from the
  # template, which takes the current context as an update; other retained
  # pools are restarted.
  spawn_alive <- identical(pool$backend, "spawn") && !is.null(pool$template) &&
    isTRUE(pool$template$alive)
  if (isTRUE(pool$retained) && !spawn_alive) {
    iterations <- pool$iterations
    base_n <- pool$base_n
    .emc_wpool_stop(pool)
    fresh <- .emc_wpool_start(target, ctx, retain_across_blocks = TRUE)
    if (is.null(fresh)) {
      if (!dead) .emc_wpool_degraded("could not re-fork the retained pool when recycling")
      return(.emc_wpool_retained_idle(ctx, iterations, base_n))
    }
    fresh$iterations <- iterations
    fresh$base_n <- base_n
    return(fresh)
  }

  if (spawn_alive) {
    .emc_wpool_stop_workers(pool)
    unlink(file.path(pool$dir,
                     c(sprintf("req%d", seq_len(target)),
                       sprintf("ans%d", seq_len(target)))))
    ctx_update <- NULL
    if (isTRUE(pool$retained)) {
      ctx_update <- file.path(pool$dir, "context_current.rds")
      if (!.emc_wpool_save_update(ctx, pool$template_data, ctx_update))
        ctx_update <- NULL
    }
    spawned <- .emc_wpool_spawn(pool$dir, seq_len(target), ctx,
                                backend = "spawn", template = pool$template,
                                ctx_update = ctx_update,
                                context_epoch = pool$context_epoch)
    if (!is.null(ctx_update)) unlink(ctx_update)
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
  if (isTRUE(pool$kept)) {
    .emc_pool_state$session_in_use <- FALSE
  } else if (identical(pool$backend, "spawn")) {
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

.emc_wpool_serve <- function(req, ans, ctx, inherited = NULL,
                             context_epoch = 0L) {
  for (cn in inherited) try(close(cn), silent = TRUE)
  if (identical(tryCatch(emc_arm_parent_death(), error = function(e) 0L), -1L)) {
    return(invisible(NULL))
  }
  rc <- fifo(req, "rb", blocking = TRUE)
  wc <- fifo(ans, "wb", blocking = TRUE)
  done <- NULL
  context_epoch <- as.integer(context_epoch)
  on.exit({ try(close(rc), silent = TRUE); try(close(wc), silent = TRUE)
            if (!is.null(done)) try(close(done), silent = TRUE) }, add = TRUE)
  repeat {
    msg <- .emc_wpool_recv(rc)
    if (is.null(msg)) break
    if (identical(msg$control, "context_update")) {
      values <- msg$values
      epoch <- suppressWarnings(as.integer(msg$epoch))
      ok <- is.list(values) && !is.null(names(values)) &&
        length(epoch) == 1L && !is.na(epoch) &&
        identical(epoch, context_epoch + 1L) &&
        !any(names(values) %in% c("data", "model"))
      if (ok) {
        ctx[names(values)] <- values
        context_epoch <- epoch
      }
      if (!.emc_wpool_send(wc, list(ok = ok, epoch = epoch))) break
      next
    }
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

.emc_wpool_template <- function(boot, base = NULL) {
  if (is.null(base)) {
    base <- readRDS(boot$ctx)
    base$model <- .restore_custom_model_container(base$model, quiet = TRUE)
  } else {
    writeLines(as.character(Sys.getpid()), boot$pidfile)
  }
  ctx <- base
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
      # A template kept across blocks takes each block's context as an update
      # over the one it booted with; recycling within a block keeps it.
      if (!is.null(command$ctx)) {
        upd <- readRDS(command$ctx)
        if (!is.null(upd$model))
          upd$model <- .restore_custom_model_container(upd$model, quiet = TRUE)
        ctx <- base
        ctx[names(upd)] <- upd
      }
      keys <- as.character(command$idx)
      old <- jobs[keys]
      old <- old[!vapply(old, is.null, logical(1))]
      if (length(old)) .emc_wpool_wait_gone(old)
      jobs[keys] <- NULL
      made <- tryCatch({
        for (w in seq_along(keys)) {
          jobs[[keys[w]]] <- local({
            i <- w
            # Detached: a worker exits as soon as its pool stops it.  An
            # attached one would wait to be collected, and the template
            # collects only when it next forks a worker in the same slot.
            parallel::mcparallel(
              .emc_wpool_serve(command$req[i], command$ans[i], ctx,
                               inherited = list(rc, wc),
                               context_epoch = if (is.null(command$context_epoch))
                                 0L else command$context_epoch),
              detached = TRUE
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
          .emc_wpool_wait_gone(born)
        }
        jobs[keys] <- NULL
        stop(conditionMessage(made))
      }
      list(ok = TRUE, pids = vapply(made, function(x) x$pid, integer(1)),
           id = command$id)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e),
                                id = command$id))
    if (!.emc_wpool_send(wc, reply)) { orphaned <- TRUE; break }
  }
  if (length(jobs)) {
    for (job in jobs) try(tools::pskill(job$pid), silent = TRUE)
    .emc_wpool_wait_gone(jobs)
  }
  if (orphaned) {
    dir <- dirname(boot$req)
    if (is.character(dir) && length(dir) == 1L && nzchar(dir)) {
      try(unlink(dir, recursive = TRUE), silent = TRUE)
    }
  }
  invisible(NULL)
}

# Wait for detached workers to exit; parallel reaps them.
.emc_wpool_wait_gone <- function(jobs) {
  pids <- vapply(jobs, function(x) as.integer(x$pid), integer(1))
  wait <- 0.0005
  while (any(vapply(pids, .emc_wpool_pid_alive, logical(1)))) {
    Sys.sleep(wait)
    wait <- min(wait * 2, 0.02)
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

# Valid log-likelihoods of each subject's rows in one native call
# (src/ll_fast_cache.cpp); NULL for a subject leaves it to the checked path.
.emc_ll_subjects_fast <- function(rows, data, model, pooled, ctx) {
  n <- length(rows)
  if (!n || !isTRUE(ctx$ll_cache_validated) || !is.null(pooled) ||
      !is.null(ctx$marginalise))
    return(vector("list", n))
  emc2_ll_cached_subjects(rows, data, model)
}

.emc_wpool_compute <- function(msg, ctx, shared = NULL) {
  if (identical(msg$kind, "subject_ll")) {
    return(.emc_with_preserved_rng({
      started <- proc.time()
      failures <- list()
      fast <- .emc_ll_subjects_fast(msg$rows, ctx$data[msg$subs], ctx$model,
                                    msg$pooled, ctx)
      ll <- lapply(seq_along(msg$subs), function(k) {
        if (!is.null(fast[[k]])) return(fast[[k]])
        tryCatch(.emc_ensemble_ll_checked(msg$rows[[k]],
                                         .emc_pooled_set_one(ctx$data[[msg$subs[k]]], msg$pooled),
                                         ctx$model,
                                         ll_cache_validated = isTRUE(ctx$ll_cache_validated)),
          error = function(e) {
            failures[[length(failures) + 1L]] <<- e
            rep(-Inf, nrow(msg$rows[[k]]))
          })
      })
      spent <- proc.time() - started
      list(ll = ll, failures = failures,
           cpu = unname(spent[["user.self"]] + spent[["sys.self"]]))
    }))
  }
  if (identical(msg$kind, "laplace_rebuild_batch")) {
    return(list(results = lapply(
      msg$tasks, .emc_laplace_rebuild_compute, ctx = ctx)))
  }
  if (identical(msg$kind, "subject_ll_claim")) {
    return(.emc_with_preserved_rng({
      started <- proc.time()
      map <- .emc_wpool_claim_attach(msg$claim_path)
      if (is.null(map)) stop("cannot map the subject queue counter")
      b <- .emc_wpool_shared_read(msg$shared_file)
      failures <- list()
      subs <- integer(0)
      ll <- list()
      repeat {
        i <- emc_claim_next(map$ptr, msg$key, msg$n)
        if (i < 0L) break
        s <- b$order[[i + 1L]]
        subs <- c(subs, s)
        ll[[length(ll) + 1L]] <- tryCatch(
          .emc_ensemble_ll_checked(b$rows[[s]],
                                   .emc_pooled_set_one(ctx$data[[s]], b$pooled),
                                   ctx$model,
                                   ll_cache_validated = isTRUE(ctx$ll_cache_validated)),
          error = function(e) {
            failures[[length(failures) + 1L]] <<- e
            rep(-Inf, nrow(b$rows[[s]]))
          })
      }
      spent <- proc.time() - started
      list(subs = subs, ll = ll, failures = failures,
           cpu = unname(spent[["user.self"]] + spent[["sys.self"]]))
    }))
  }
  if (identical(msg$kind, "hmc")) return(.emc_hmc_compute(msg, ctx))
  if (identical(msg$kind, "claim")) return(.emc_wpool_compute_claim(msg, ctx))
  if (identical(msg$kind, "ll")) {
    return(list(ll = calc_ll_manager(msg$proposals, dadm = ctx$data[[msg$s]],
                                     model = ctx$model, component = msg$component,
                                     r_cores = 1L, varying = msg$varying,
                                     ll_cache_validated = isTRUE(ctx$ll_cache_validated))))
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

# Evaluate each subject's candidate rows (the same number for every subject)
# on the chain's existing workers: a vector with one row each, else a
# subjects x rows matrix. Replies keep subject order; missing replies are
# recomputed on the master before the decision, so transport failure cannot
# change it. Worker CPU accumulates in .emc_pool_state$subject_ll_cpu.
# `source` names the move whose likelihood failures are recorded.
.emc_wpool_subject_ll <- function(pool, rows, ctx, part = NULL, source) {
  n <- length(rows)
  if (is.null(pool) || !isTRUE(pool$alive)) {
    # Without workers, subjects served by the stored-argument path need none
    # of the dispatch below.
    if (n > 0L) {
      pooled <- if (is.data.frame(ctx$data[[1L]])) attr(ctx$data[[1L]], "pooled_values")
      fast <- .emc_ll_subjects_fast(rows, ctx$data, ctx$model, pooled, ctx)
      if (length(fast) == n && !any(vapply(fast, is.null, logical(1)))) {
        if (is.null(.emc_pool_state$subject_ll_cpu)) .emc_pool_state$subject_ll_cpu <- 0
        ll <- if (nrow(rows[[1L]]) == 1L) unlist(fast, use.names = FALSE) else
          do.call(rbind, fast)
        return(list(ll = ll, pool = pool))
      }
    }
    part <- list(seq_len(n))
  } else if (is.null(part)) {
    part <- .emc_lpt_partition(.emc_subject_cost(ctx$data), pool$n)
  }
  if (sum(lengths(part)) != n || !setequal(unlist(part), seq_len(n)) ||
      (isTRUE(pool$alive) && length(part) != pool$n))
    stop("Subject likelihood partition must cover each subject exactly once")
  # Workers hold the data of their template, so the current pooled values
  # travel with each request.
  pooled <- if (is.data.frame(ctx$data[[1L]])) attr(ctx$data[[1L]], "pooled_values")
  claim_map <- if (isTRUE(pool$alive) && pool$n > 1L && !is.null(pool$dir) &&
                   isTRUE(getOption("emc2.worker_queue", TRUE))) {
    .emc_wpool_claim_attach(file.path(pool$dir, "claim"))
  } else NULL
  # Over the shared subject queue every subject's rows go out once in a
  # broadcast file, and each worker claims the next subject (costliest first)
  # from the pool's counter as it frees up, replying once when the queue is
  # empty. A subject no valid reply holds is computed here.
  if (!is.null(claim_map)) {
    queued <- emc_wpool_subject_ll_queue_cpp(pool, rows, ctx, part, source, topenv())
    if (!is.null(queued)) return(queued)
  }
  msgs <- lapply(part, function(subs)
    list(kind = "subject_ll", subs = subs, rows = rows[subs], notify = FALSE,
         pooled = pooled))
  sent <- logical(length(part))
  if (isTRUE(pool$alive)) {
    for (w in seq_along(part)) {
      if (!length(part[[w]])) next
      sent[w] <- .emc_wpool_request(pool, w, msgs[[w]])
      if (!sent[w]) {
        .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
        pool$alive <- FALSE
        break
      }
    }
  }
  replies <- vector("list", length(part))
  deadline <- .emc_wpool_deadline()
  k <- nrow(rows[[1L]])
  valid <- function(res, subs) {
    is.list(res$ll) && length(res$ll) == length(subs) &&
      all(vapply(res$ll, function(x) is.numeric(x) && length(x) == k &&
                   !anyNA(x) && !any(x == Inf), logical(1)))
  }
  worker_cpu <- 0
  for (w in which(sent)) {
    res <- .emc_wpool_reply(pool, w, deadline)
    if (is.null(res) || !is.null(res$failed) || !valid(res, part[[w]])) {
      .emc_wpool_degraded("subject likelihood worker gave no valid reply")
      pool$alive <- FALSE
    } else {
      replies[[w]] <- res
      if (is.numeric(res$cpu)) worker_cpu <- worker_cpu + res$cpu
    }
  }
  ll <- matrix(0, n, k)
  for (w in seq_along(part)) {
    if (!length(part[[w]])) next
    res <- replies[[w]]
    if (is.null(res)) res <- .emc_wpool_compute(msgs[[w]], ctx)
    for (e in res$failures) {
      cls <- .emc_reject_record(e, source)
      if (identical(.emc_failure_action(cls), "abort"))
        .emc_failure_abort(cls, source, e)
    }
    ll[part[[w]], ] <- do.call(rbind, res$ll)
  }
  .emc_pool_state$subject_ll_cpu <-
    (if (is.null(.emc_pool_state$subject_ll_cpu)) 0 else .emc_pool_state$subject_ll_cpu) +
    worker_cpu
  list(ll = if (k == 1L) ll[, 1L] else ll, pool = pool)
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

.emc_wpool_shared_chol <- function(shared) {
  group_chol <- shared$group_chol
  if (is.null(group_chol) && !is.null(shared$group_var)) {
    group_chol <- build_group_chol_cache(shared$group_var, shared$idx_list,
                                         shared$marginal_idx,
                                         include_full = if (is.null(shared$include_full))
                                           TRUE else shared$include_full)
  }
  group_chol
}

# Map the shared subject-queue counter at `path` (src/claim_counter.cpp),
# once per process; NULL where it cannot be mapped.
.emc_wpool_claim_attach <- function(path) {
  cached <- .emc_pool_state$claim_map
  if (!is.null(cached) && identical(cached$path, path)) return(cached)
  ptr <- tryCatch(emc_claim_open(path), error = function(e) NULL)
  if (is.null(ptr)) return(NULL)
  .emc_pool_state$claim_map <- list(path = path, ptr = ptr)
  .emc_pool_state$claim_map
}

# One worker's share of the shared subject queue: claim subjects from the
# counter until it runs out, then answer for all of them at once.
.emc_wpool_compute_claim <- function(msg, ctx) {
  t_recv <- Sys.time()
  started <- proc.time()
  map <- .emc_wpool_claim_attach(msg$claim_path)
  if (is.null(map)) stop("cannot map the subject queue counter")
  b <- .emc_wpool_shared_read(msg$shared_file)
  shared <- list(group_var = b$group_var, group_chol = .emc_wpool_shared_chol(b))
  inp <- b$inputs
  one_subject <- function(s)
    list(subs = s, alpha = inp$alpha[, s, drop = FALSE],
         population_mu = inp$population_mu[, s, drop = FALSE],
         pm = inp$pm[s], prev_ll = inp$prev_ll[s], seeds = inp$seeds[s],
         extra = if (is.null(inp$extra)) NULL else inp$extra[s],
         pooled = inp$pooled, rows_only = inp$rows_only,
         skip = if (is.null(inp$skip)) NULL else inp$skip[s])
  whole <- .emc_wpool_claim_native(msg, ctx, map, inp, shared, one_subject)
  if (!is.null(whole)) {
    spent <- proc.time() - started
    return(c(whole, list(
      t_recv = t_recv, t_done = Sys.time(),
      cpu = unname(spent[["user.self"]] + spent[["sys.self"]]),
      elapsed = unname(spent[["elapsed"]]), rejects = whole$rejects,
      ensemble_rejects = whole$ensemble_rejects))[
        c("subs", "props", "pm", "seeds", "times", "extra_ll", "t_recv", "t_done", "cpu",
          "elapsed", "rejects", "ensemble_rejects")])
  }
  got <- integer(0)
  outs <- list()
  repeat {
    i <- emc_claim_next(map$ptr, msg$key, msg$n)
    if (i < 0L) break
    s <- inp$order[[i + 1L]]
    outs[[length(outs) + 1L]] <- .emc_wpool_compute_particle(one_subject(s), ctx, shared)
    got <- c(got, s)
  }
  spent <- proc.time() - started
  rejects <- if (length(outs)) Reduce(`+`, lapply(outs, `[[`, "rejects"))
             else stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                                  .EMC_REJECT_CLASSES)
  ensemble_rejects <- list(
    counts = stats::setNames(integer(length(.EMC_REJECT_CLASSES)), .EMC_REJECT_CLASSES),
    msgs = list())
  for (out in outs) if (!is.null(out$ensemble_rejects)) {
    ensemble_rejects$counts <- ensemble_rejects$counts + out$ensemble_rejects$counts
    for (cls in names(out$ensemble_rejects$msgs))
      if (is.null(ensemble_rejects$msgs[[cls]]))
        ensemble_rejects$msgs[[cls]] <- out$ensemble_rejects$msgs[[cls]]
  }
  list(subs = got,
       props = if (length(outs)) do.call(cbind, lapply(outs, `[[`, "props"))
               else matrix(0, ctx$n_pars + 1L, 0L),
       pm = do.call(c, lapply(outs, `[[`, "pm")),
       seeds = do.call(c, lapply(outs, `[[`, "seeds")),
       times = as.numeric(unlist(lapply(outs, `[[`, "times"))),
       extra_ll = if (is.null(inp$extra)) NULL else do.call(c, lapply(outs, `[[`, "extra_ll")),
       t_recv = t_recv, t_done = Sys.time(),
       cpu = unname(spent[["user.self"]] + spent[["sys.self"]]),
       elapsed = unname(spent[["elapsed"]]), rejects = rejects,
       ensemble_rejects = ensemble_rejects)
}

# The same share with the subjects claimed inside the native particle kernel:
# one call steps every subject this worker claims and evaluates its candidate
# rows. A claimed subject the kernel leaves to R runs through
# .emc_wpool_compute_particle() from its untouched stream, and candidates
# the kernel does not evaluate take the checked path.
# NULL, before anything is claimed, when the round cannot run this way or no
# subject is left.
.emc_wpool_claim_native <- function(msg, ctx, map, inp, shared, one_subject) {
  n <- NCOL(inp$alpha)
  if (!isTRUE(getOption("emc2.native_claim", TRUE)) || isTRUE(inp$rows_only) ||
      !is.matrix(inp$alpha) || n < 1L || length(inp$pm) != n || length(inp$seeds) != n ||
      length(inp$prev_ll) != n || length(ctx$data) != n ||
      (!is.null(inp$extra) && length(inp$extra) != n) ||
      (!is.null(inp$skip) && length(inp$skip) != n))
    return(NULL)
  all_subjects <- list(subs = seq_len(n), alpha = inp$alpha,
                       population_mu = inp$population_mu, pm = inp$pm, skip = inp$skip,
                       pooled = inp$pooled, rows_only = inp$rows_only)
  args <- .emc_native_particle_inputs(all_subjects, ctx, shared, shared$group_chol,
                                      statics = FALSE)
  if (is.null(args)) return(NULL)
  # A subject's block-fixed parts are built when this worker claims it.
  par_names <- rownames(inp$alpha)
  get_static <- function(s)
    tryCatch(.emc_native_static(s, ctx, par_names), error = function(e) NULL)
  extra_cached <- isTRUE(ctx$ll_cache_validated) && is.null(inp$pooled) &&
    is.null(ctx$marginalise)
  res <- .emc_with_preserved_rng(emc2_native_particle_claim(
    get_static, inp$alpha, as.numeric(inp$prev_ll), inp$population_mu, args$group_var,
    args$group, inp$pm, inp$seeds, args$skip, args$stage, args$tune_in, map$ptr,
    as.integer(msg$key), as.integer(unlist(inp$order))[seq_len(msg$n)], inp$extra,
    extra_cached, ctx$data, ctx$model))
  got <- res$claimed
  if (!length(got)) return(NULL)
  res <- .emc_native_particle_settings(res, inp$pm, args$stage)
  done <- res$done[got] %in% TRUE
  skipped <- !done & args$skip[got]
  complete <- done | skipped
  unscored <- if (is.null(inp$extra)) logical(length(got)) else
    complete & !vapply(inp$extra[got], is.null, logical(1)) &
      vapply(res$extra_ll[got], is.null, logical(1))
  rows <- seq_len(ctx$n_pars)
  props <- matrix(0, ctx$n_pars + 1L, length(got))
  pm <- vector("list", length(got))
  seeds <- vector("list", length(got))
  times <- numeric(length(got))
  extra_ll <- if (is.null(inp$extra)) NULL else vector("list", length(got))
  props[rows, done] <- res$alpha[, got[done]]
  props[ctx$n_pars + 1L, done] <- res$ll[got[done]]
  pm[done] <- res$pm_settings[got[done]]
  seeds[done] <- res$seeds[got[done]]
  times[done] <- res$cpu[got[done]]
  props[rows, skipped] <- inp$alpha[, got[skipped]]
  props[ctx$n_pars + 1L, skipped] <- inp$prev_ll[got[skipped]]
  pm[skipped] <- inp$pm[got[skipped]]
  seeds[skipped] <- inp$seeds[got[skipped]]
  if (!is.null(extra_ll)) extra_ll[complete] <- res$extra_ll[got[complete]]
  rejects <- stats::setNames(integer(length(.EMC_REJECT_CLASSES)), .EMC_REJECT_CLASSES)
  ensemble_rejects <- list(counts = rejects, msgs = list())
  note_ensemble <- function(counts, msgs) {
    ensemble_rejects$counts <<- ensemble_rejects$counts + counts
    for (cls in names(msgs))
      if (is.null(ensemble_rejects$msgs[[cls]]))
        ensemble_rejects$msgs[[cls]] <<- msgs[[cls]]
  }
  # Claim order, as the replies' rejections are reported.
  for (j in which(!complete | unscored)) {
    s <- got[j]
    if (complete[j]) {
      before <- .emc_reject_counts("ensemble")
      extra_ll[[j]] <- .emc_ensemble_ll_candidate(
        inp$extra[[s]], .emc_pooled_set_one(ctx$data[[s]], inp$pooled), ctx$model,
        ll_cache_validated = isTRUE(ctx$ll_cache_validated))
      delta <- .emc_reject_delta(before, "ensemble")
      if (any(delta > 0L)) {
        recorded <- .emc_reject_messages()
        note_ensemble(delta, stats::setNames(
          lapply(names(delta), function(cls) recorded[[paste("ensemble", cls, sep = ".")]]),
          names(delta))[delta > 0L])
      }
      next
    }
    out <- .emc_wpool_compute_particle(one_subject(s), ctx, shared)
    props[, j] <- out$props
    pm[j] <- out$pm
    seeds[j] <- out$seeds
    times[j] <- out$times
    if (!is.null(extra_ll)) extra_ll[j] <- out$extra_ll
    rejects <- rejects + out$rejects
    note_ensemble(out$ensemble_rejects$counts, out$ensemble_rejects$msgs)
  }
  list(subs = got, props = props, pm = pm, seeds = seeds, times = times,
       extra_ll = extra_ll, rejects = rejects, ensemble_rejects = ensemble_rejects)
}

.emc_wpool_compute_particle <- function(msg, ctx, shared = NULL) .emc_with_preserved_rng({
  if (is.null(shared)) {
    shared <- if (!is.null(msg$shared_file)) {
      .emc_wpool_shared_read(msg$shared_file)
    } else {
      unserialize(msg$shared)
    }
  }
  group_chol <- .emc_wpool_shared_chol(shared)
  subs <- msg$subs
  props <- matrix(0, ctx$n_pars + 1L, length(subs))
  pm <- vector("list", length(subs))
  seeds <- vector("list", length(subs))
  times <- numeric(length(subs))
  extra_ll <- if (is.null(msg$extra)) NULL else vector("list", length(subs))
  t_recv <- Sys.time()
  w_started <- proc.time()
  rejects_before <- .emc_reject_counts("particle")
  ensemble_rejects_before <- .emc_reject_counts("ensemble")
  native <- .emc_native_particle_batch(msg, ctx, shared, group_chol)
  extra_fast <- if (!is.null(extra_ll) && !isTRUE(msg$rows_only))
    .emc_ll_subjects_fast(msg$extra, ctx$data[subs], ctx$model, msg$pooled, ctx)
  extra_candidate <- function(k, data_s) {
    if (!is.null(extra_fast[[k]])) return(extra_fast[[k]])
    .emc_ensemble_ll_candidate(msg$extra[[k]], data_s, ctx$model,
                               ll_cache_validated = isTRUE(ctx$ll_cache_validated))
  }
  # Every subject either stepped natively or was skipped, with its candidate
  # likelihoods already in hand: the results are assembled whole.
  whole <- FALSE
  if (!is.null(native) && !isTRUE(msg$rows_only)) {
    done <- native$done %in% TRUE
    skipped <- !done & (if (is.null(msg$skip)) FALSE else msg$skip %in% TRUE)
    whole <- all(done | skipped) &&
      (is.null(extra_ll) || (length(extra_fast) == length(subs) &&
                               !any(vapply(extra_fast, is.null, logical(1))) &&
                               !any(vapply(msg$extra, is.null, logical(1)))))
  }
  if (whole) {
    rows <- seq_len(ctx$n_pars)
    props[rows, done] <- native$alpha[, done]
    props[ctx$n_pars + 1L, done] <- native$ll[done]
    pm[done] <- native$pm_settings[done]
    seeds[done] <- native$seeds[done]
    times[done] <- native$cpu[done]
    if (any(skipped)) {
      props[rows, skipped] <- msg$alpha[, skipped]
      props[ctx$n_pars + 1L, skipped] <- msg$prev_ll[skipped]
      pm[skipped] <- msg$pm[skipped]
      seeds[skipped] <- msg$seeds[skipped]
    }
    if (!is.null(extra_ll)) extra_ll <- extra_fast
  } else for (k in seq_along(subs)) {
    s <- subs[k]
    data_s <- .emc_pooled_set_one(ctx$data[[s]], msg$pooled)
    if (isTRUE(native$done[k])) {
      props[, k] <- c(native$alpha[, k], native$ll[k])
      pm[[k]] <- native$pm_settings[[k]]
      seeds[[k]] <- native$seeds[[k]]
      times[k] <- native$cpu[k]
      if (!is.null(extra_ll) && !is.null(msg$extra[[k]]))
        extra_ll[[k]] <- extra_candidate(k, data_s)
      next
    }
    if (isTRUE(msg$rows_only)) {
      # Pooled-parameter candidates (R/pooled.R): no particle step, no RNG.
      props[, k] <- c(msg$alpha[, k], msg$prev_ll[k])
      pm[[k]] <- msg$pm[[k]]
      seeds[[k]] <- msg$seeds[[k]]
      extra_ll[[k]] <- .emc_pooled_rows_ll(
        msg$extra[[k]], data_s, ctx$model,
        ll_cache_validated = isTRUE(ctx$ll_cache_validated))
      next
    }
    # A subject whose particle step the allocation skipped this iteration
    # (R/move_allocation.R) keeps its state, settings and stream.
    if (isTRUE(msg$skip[k])) {
      props[, k] <- c(msg$alpha[, k], msg$prev_ll[k])
      pm[[k]] <- msg$pm[[k]]
      seeds[[k]] <- msg$seeds[[k]]
      if (!is.null(extra_ll) && !is.null(msg$extra[[k]]))
        extra_ll[[k]] <- extra_candidate(k, data_s)
      next
    }
    assign(".Random.seed", msg$seeds[[k]], envir = globalenv())
    t0 <- sum(proc.time()[c("user.self", "sys.self")])
    out <- safe_new_particle(
      s = s, data = data_s, pm_settings = msg$pm[[k]],
      eff_mu = ctx$eff_mu[[s]], eff_var = ctx$eff_var[[s]],
      chains_mu = ctx$chains_mu[[s]], chains_var = ctx$chains_var[[s]],
      prev_ll = msg$prev_ll[k], parameters = NULL, model = ctx$model,
      stage = ctx$stage, type = ctx$type, tune = ctx$tune,
      marginalise = ctx$marginalise, r_cores = ctx$r_cores,
      chol_cache = ctx$chol_caches[[s]], group_chol = group_chol,
      current_alpha = .emc_subject_col(msg$alpha, k),
      population_mu = msg$population_mu[, k],
      population_var = shared$group_var,
      ll_cache_validated = isTRUE(ctx$ll_cache_validated)
    )
    times[k] <- sum(proc.time()[c("user.self", "sys.self")]) - t0
    props[, k] <- c(out$proposal, out$ll)
    pm[[k]] <- out$pm_settings
    seeds[[k]] <- get(".Random.seed", envir = globalenv())
    # Ensemble pool members for this subject (no RNG use).
    if (!is.null(extra_ll) && !is.null(msg$extra[[k]]))
      extra_ll[[k]] <- extra_candidate(k, data_s)
  }
  w_spent <- proc.time() - w_started
  ensemble_delta <- .emc_reject_delta(ensemble_rejects_before, "ensemble")
  ensemble_msgs <- stats::setNames(list(), character(0))
  if (any(ensemble_delta > 0L)) {
    reject_messages <- .emc_reject_messages()
    ensemble_msgs <- lapply(names(ensemble_delta), function(cls)
      reject_messages[[paste("ensemble", cls, sep = ".")]])
    names(ensemble_msgs) <- names(ensemble_delta)
    ensemble_msgs <- ensemble_msgs[ensemble_delta > 0L]
  }
  list(props = props, pm = pm, seeds = seeds, times = times, extra_ll = extra_ll,
       t_recv = t_recv, t_done = Sys.time(),
       cpu = unname(w_spent[["user.self"]] + w_spent[["sys.self"]]),
       elapsed = unname(w_spent[["elapsed"]]),
       rejects = .emc_reject_delta(rejects_before, "particle"),
       ensemble_rejects = list(counts = ensemble_delta, msgs = ensemble_msgs))
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
  # Decisions use the latest .EMC_LL_PROBE_N timings. Without a pool the
  # mode stays "probe" for the life of the process, so the record is capped.
  keep <- function(x) {
    x <- c(x, per)
    n <- length(x)
    if (n > .EMC_LL_PROBE_N) x[(n - .EMC_LL_PROBE_N + 1L):n] else x
  }
  if (pooled) st$t_pool <- keep(st$t_pool) else st$t_serial <- keep(st$t_serial)
  if (pooled && identical(st$mode, "pool")) {
    if (dynamic) st$t_dynamic <- keep(st$t_dynamic)
    else st$t_static <- keep(st$t_static)
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
                            prev_ll, seeds, extra = NULL, pooled = NULL,
                            rows_only = FALSE, skip = NULL) {
  profile <- isTRUE(getOption("emc2.sampler_profile", FALSE))
  iter_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  n_subjects <- length(pm_settings)
  props <- matrix(0, ctx$n_pars + 1L, n_subjects)
  times <- numeric(n_subjects)
  cpu_total <- 0
  extra_ll <- if (is.null(extra)) NULL else vector("list", n_subjects)

  population_mu <- if (identical(ctx$type, "standard") && is.matrix(pars$subj_mu) &&
                       is.double(pars$subj_mu) &&
                       identical(dim(pars$subj_mu), c(as.integer(ctx$n_pars), n_subjects))) {
    matrix(as.numeric(pars$subj_mu), ctx$n_pars, n_subjects)
  } else matrix(vapply(seq_len(n_subjects), function(s) {
    as.numeric(get_group_level(pars, s, ctx$type)$mu)
  }, numeric(ctx$n_pars)), ctx$n_pars, n_subjects)
  alpha <- matrix(pars$alpha, ctx$n_pars, n_subjects, dimnames = dimnames(pars$alpha))
  rownames(population_mu) <- rownames(alpha)
  shared_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  shared <- list(group_var = group_chol$ref, idx_list = group_chol$idx_list,
                 marginal_idx = group_chol$marginal_idx,
                 include_full = group_chol$include_full,
                 group_chol = group_chol)
  was_alive <- isTRUE(pool$alive)
  watch <- !is.null(pool$done)
  # Shared subject queue: per-subject cost varies between iterations with the
  # proposals, so a partition planned from last iteration's times leaves
  # workers idle while the unluckiest finishes.  Every subject's inputs go out
  # once in the broadcast, and each worker claims the next subject (costliest
  # first) from a shared counter as it frees up, replying once when the queue
  # is empty.  Each subject carries its own RNG stream, so the result does not
  # depend on which worker computes it.
  claim_map <- if (was_alive && watch && pool$n > 1L &&
                   isTRUE(getOption("emc2.worker_queue", TRUE))) {
    .emc_wpool_claim_attach(file.path(pool$dir, "claim"))
  } else NULL
  claim <- !is.null(claim_map)
  # The unprofiled round over the queue runs natively (src/pool_exchange.cpp).
  if (claim && !profile) {
    queued <- emc_wpool_particle_queue_cpp(pool, ctx, part, alpha, population_mu,
                                           group_chol, pm_settings, prev_ll, seeds, extra,
                                           pooled, skip, as.integer(ctx$n_pars),
                                           isTRUE(rows_only), topenv())
    if (!is.null(queued)) return(queued)
  }
  order <- if (claim) .emc_wpool_queue_order(part) else NULL
  shared_raw <- if (was_alive) {
    serialize(c(shared[c("group_var", "idx_list", "marginal_idx", "include_full")],
                if (claim) list(inputs = list(
                  alpha = alpha, population_mu = population_mu,
                  pm = pm_settings, prev_ll = prev_ll, seeds = seeds,
                  extra = extra, pooled = pooled, rows_only = rows_only,
                  skip = skip, order = order))),
              NULL)
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
  generation <- .emc_wpool_one_int(pool$generation)
  request <- .emc_wpool_next_request()
  subject_msg <- function(subs, w) {
    list(subs = subs, shared_file = shared_file,
         alpha = alpha[, subs, drop = FALSE],
         population_mu = population_mu[, subs, drop = FALSE],
         pm = pm_settings[subs], prev_ll = prev_ll[subs], seeds = seeds[subs],
         extra = if (is.null(extra)) NULL else extra[subs],
         pooled = pooled, rows_only = rows_only,
         skip = if (is.null(skip)) NULL else skip[subs],
         notify = watch, w = w, generation = generation, request = request)
  }
  msgs <- lapply(seq_len(pool$n), function(w) {
    if (claim) {
      list(kind = "claim", shared_file = shared_file, claim_path = claim_map$path,
           key = request, n = length(order), notify = TRUE, w = w,
           generation = generation, request = request)
    } else subject_msg(part[[w]], w)
  })
  has_work <- if (claim) rep(TRUE, pool$n) else lengths(part) > 0L
  sent <- rep(FALSE, pool$n)
  dispatch_at <- finish_at <- rep(NA_real_, pool$n)
  dispatch_wall <- if (profile) rep(Sys.time()[NA], pool$n) else NULL
  worker_cpu <- worker_elapsed <- queue_delay <- rep(NA_real_, pool$n)
  fallback_workers <- 0L
  fallback_subjects <- 0L
  rejects <- .EMC_REJECT_ZERO
  ensemble_rejects <- list(counts = .EMC_REJECT_ZERO, msgs = list())
  send_started <- if (profile) proc.time()[["elapsed"]] else NA_real_
  if (was_alive && !is.null(shared_file)) {
    if (claim) emc_claim_reset(claim_map$ptr, request)
    send_deadline <- .emc_wpool_deadline()
    for (w in seq_len(pool$n)) {
      if (!has_work[w]) next
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
  queue_delays <- numeric(0)
  assigned <- if (claim) vector("list", pool$n) else part
  stored <- logical(n_subjects)
  store <- function(w, res, source = "worker", subs = part[[w]]) {
    props[, subs] <<- res$props
    times[subs] <<- res$times
    pm_settings[subs] <<- res$pm
    seeds[subs] <<- res$seeds
    stored[subs] <<- TRUE
    if (!is.null(extra) && !is.null(res$extra_ll)) extra_ll[subs] <<- res$extra_ll
    if (!is.null(res$rejects)) rejects <<- rejects + res$rejects
    if (!is.null(res$ensemble_rejects)) {
      ensemble_rejects$counts <<- ensemble_rejects$counts + res$ensemble_rejects$counts
      for (cls in names(res$ensemble_rejects$msgs))
        if (is.null(ensemble_rejects$msgs[[cls]]))
          ensemble_rejects$msgs[[cls]] <<- res$ensemble_rejects$msgs[[cls]]
    }
    if (identical(source, "master") && was_alive) {
      fallback_workers <<- fallback_workers + 1L
      fallback_subjects <<- fallback_subjects + length(subs)
    }
    if (identical(source, "worker") && is.numeric(res$cpu)) cpu_total <<- cpu_total + res$cpu
    if (claim && !is.na(w)) assigned[[w]] <<- c(assigned[[w]], subs)
    if (!profile || is.na(w)) return(invisible(NULL))
    finish_at[w] <<- proc.time()[["elapsed"]] - iter_started
    if (!is.null(res$cpu)) worker_cpu[w] <<- res$cpu
    if (!is.null(res$elapsed)) worker_elapsed[w] <<- res$elapsed
    if (identical(source, "worker") && !is.null(res$t_recv) &&
        !is.na(dispatch_wall[w])) {
      queue_delay[w] <<- as.numeric(difftime(res$t_recv, dispatch_wall[w],
                                             units = "secs"))
      queue_delays <<- c(queue_delays, queue_delay[w])
    }
    invisible(NULL)
  }
  # A claim reply must name subjects of this iteration that no other reply
  # holds, with one result column each.
  claim_valid <- function(res) {
    subs <- res$subs
    is.numeric(subs) && !anyNA(subs) && all(subs >= 1L & subs <= n_subjects) &&
      !anyDuplicated(subs) && !any(stored[subs]) &&
      is.matrix(res$props) && ncol(res$props) == length(subs) &&
      length(res$pm) == length(subs) && length(res$seeds) == length(subs)
  }
  # Under the shared queue the subjects a failed worker had claimed are not
  # known, so they are recomputed with the unclaimed ones after the replies.
  recompute <- function(w) {
    if (!claim) store(w, .emc_wpool_compute(msgs[[w]], ctx, shared = shared),
                      source = "master")
  }
  stop_claims <- function() if (claim) emc_claim_reset(claim_map$ptr, 0L)

  pending <- sent & has_work
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
        stop_claims()
        if (identical(rec$code, "timeout")) {
          stuck <- which(pending)
          idle <- stuck[!picked_up[stuck]]
          stuck_subs <- if (claim) integer(0) else unlist(part[stuck])
          .emc_wpool_degraded(sprintf(
            "worker%s %s produced no reply within %.0f s%s%s",
            if (length(stuck) > 1L) "s" else "", paste(stuck, collapse = ", "),
            deadline,
            if (length(stuck_subs)) {
              sprintf(" (subject%s %s)", if (length(stuck_subs) > 1L) "s" else "",
                      paste(stuck_subs, collapse = ", "))
            } else "",
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
        for (k in which(pending)) recompute(k)
        pending[] <- FALSE
        break
      }
      res <- .emc_wpool_reply(pool, w, deadline)
      pending[w] <- FALSE
      if (is.null(res)) {
        stop_claims()
        .emc_wpool_degraded(.emc_wpool_transport_error("worker gave no reply"))
        pool$alive <- FALSE
        recompute(w)
      } else if (!is.null(res$failed) || (claim && !claim_valid(res))) {
        recompute(w)
      } else if (claim) {
        store(w, res, subs = as.integer(res$subs))
      } else {
        store(w, res)
      }
    }
    for (w in seq_len(pool$n)) if (has_work[w] && !sent[w]) recompute(w)
  } else {
    deadline <- .emc_wpool_deadline()
    for (w in seq_len(pool$n)) {
      if (!has_work[w]) next
      res <- if (sent[w]) .emc_wpool_reply(pool, w, deadline) else NULL
      if (is.null(res) || !is.null(res$failed)) {
        if (sent[w] && is.null(res)) {
          .emc_wpool_degraded("worker gave no reply")
          pool$alive <- FALSE
        }
        recompute(w)
      } else {
        store(w, res)
      }
    }
  }
  # Subjects no worker answered for (claimed by a worker that failed, or never
  # claimed because the pool degraded) run here.
  if (claim) {
    for (s in which(!stored)) {
      store(NA_integer_, .emc_wpool_compute(subject_msg(s, NA_integer_), ctx,
                                            shared = shared),
            source = "master", subs = s)
    }
  }
  receive_elapsed <- proc.time()[["elapsed"]] - receive_started
  if (was_alive) .emc_wpool_record_receive(receive_elapsed)
  measure_wire <- .emc_profile_expensive()
  private_bytes <- if (measure_wire) {
    sum(vapply(msgs[sent], function(msg) {
      msg$shared <- NULL
      msg$shared_file <- NULL
      length(serialize(msg, NULL))
    }, integer(1)))
  } else NA_real_
  shared_path_bytes <- if (measure_wire && !is.null(shared_file)) {
    sum(vapply(msgs[sent], function(msg) {
      length(serialize(list(shared_file = msg$shared_file), NULL))
    }, integer(1)))
  } else if (measure_wire) 0 else NA_real_
  worker_load <- vapply(assigned, function(subs) sum(times[subs]), numeric(1))
  list(props = props, pm_settings = pm_settings, seeds = seeds, times = times,
       extra_ll = extra_ll, alive = pool$alive, rejects = rejects,
       ensemble_rejects = ensemble_rejects, worker_cpu = cpu_total,
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
         workers_active = sum(vapply(assigned, length, integer(1)) > 0L),
         dispatch_first = .emc_range_or_na(dispatch_at, min),
         dispatch_last = .emc_range_or_na(dispatch_at, max),
         finish_first = .emc_range_or_na(finish_at, min),
         finish_last = .emc_range_or_na(finish_at, max),
         queue_delay_max = .emc_range_or_na(queue_delays, max),
         queue_delay_mean = .emc_range_or_na(queue_delays, mean),
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

# What the native exchange (src/pool_exchange.cpp) leaves to R when the pool
# fails part-way: the lost pool's warning with the workers that never
# answered, which are stopped; one subject's step in this process; and the
# failures the replies carried.
.emc_wpool_iter_timeout <- function(pool, stuck, idle, deadline) {
  .emc_wpool_degraded(sprintf(
    "worker%s %s produced no reply within %.0f s%s",
    if (length(stuck) > 1L) "s" else "", paste(stuck, collapse = ", "), deadline,
    if (length(idle)) {
      sprintf("; worker%s %s never collected the request",
              if (length(idle) > 1L) "s" else "", paste(idle, collapse = ", "))
    } else ""))
  .emc_wpool_terminate_jobs(pool$jobs[stuck], wait = TRUE, terminate = TRUE)
  invisible(NULL)
}

.emc_wpool_iter_missing <- function(s, ctx, alpha, population_mu, pm_settings, prev_ll,
                                    seeds, extra, pooled, skip, group_chol,
                                    rows_only = FALSE) {
  shared <- list(group_var = group_chol$ref, idx_list = group_chol$idx_list,
                 marginal_idx = group_chol$marginal_idx,
                 include_full = group_chol$include_full, group_chol = group_chol)
  .emc_wpool_compute(
    list(subs = s, alpha = alpha[, s, drop = FALSE],
         population_mu = population_mu[, s, drop = FALSE], pm = pm_settings[s],
         prev_ll = prev_ll[s], seeds = seeds[s],
         extra = if (is.null(extra)) NULL else extra[s], pooled = pooled,
         rows_only = rows_only, skip = if (is.null(skip)) NULL else skip[s]),
    ctx, shared = shared)
}

.emc_wpool_note_failures <- function(failures, source) {
  for (e in failures) {
    cls <- .emc_reject_record(e, source)
    if (identical(.emc_failure_action(cls), "abort"))
      .emc_failure_abort(cls, source, e)
  }
  invisible(NULL)
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
  out <- lapply(out, function(x) if (is.null(x)) integer(0) else x)
  # The shared subject queue in .emc_wpool_iter() hands subjects out in this
  # order (costliest first), so the tail of an iteration is the cheapest work.
  attr(out, "order") <- order(cost, decreasing = TRUE)
  out
}

# Costliest-first dispatch order for the subjects in `part`: the partitioner's
# recorded order when present, restricted to (and completed by) the subjects
# the partition actually holds.
.emc_wpool_queue_order <- function(part) {
  subs <- unlist(part, use.names = FALSE)
  o <- attr(part, "order")
  if (is.null(o)) return(subs)
  o <- o[o %in% subs]
  c(o, setdiff(subs, o))
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
