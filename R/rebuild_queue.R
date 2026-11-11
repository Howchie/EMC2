# Rebuild jobs carry their own RNG stream or are deterministic. A shared
# counter balances uneven jobs without changing their result order. Forks
# inherit the inputs and exit after each rebuild, so model closures do not
# accumulate in the persistent sampling workers.
.emc_rebuild_queue <- new.env(parent = emptyenv())
.emc_rebuild_queue$costs <- list()

.emc_rebuild_cost_id <- function(cost_key, cores) {
  paste(as.character(cost_key)[1L], as.integer(cores)[1L], sep = "\r")
}

.emc_rebuild_makespan <- function(job_seconds, cores) {
  job_seconds <- as.numeric(job_seconds)
  if (!length(job_seconds)) return(0)
  loads <- numeric(max(1L, min(as.integer(cores), length(job_seconds))))
  for (seconds in sort(job_seconds, decreasing = TRUE)) {
    worker <- which.min(loads)
    loads[worker] <- loads[worker] + seconds
  }
  max(loads)
}

.emc_rebuild_prediction <- function(cost_key, n, cores) {
  if (is.null(cost_key) || !length(cost_key) || !nzchar(cost_key[1L]))
    return(NULL)
  model <- .emc_rebuild_queue$costs[[.emc_rebuild_cost_id(cost_key, cores)]]
  if (!is.list(model) || !length(model$job_seconds)) return(NULL)
  jobs <- model$job_seconds
  if (length(jobs) != n) jobs <- rep(stats::median(jobs), n)
  list(serial = sum(jobs),
       parallel = .emc_rebuild_makespan(jobs, cores) + model$startup_seconds,
       jobs = jobs,
       startup_seconds = model$startup_seconds)
}

.emc_rebuild_update_cost <- function(cost_key, cores, job_seconds,
                                     startup_seconds = NULL) {
  if (is.null(cost_key) || !length(cost_key) || !nzchar(cost_key[1L]) ||
      !length(job_seconds) || any(!is.finite(job_seconds)) ||
      any(job_seconds < 0)) return(invisible(NULL))
  id <- .emc_rebuild_cost_id(cost_key, cores)
  old <- .emc_rebuild_queue$costs[[id]]
  if (is.list(old) && length(old$job_seconds) == length(job_seconds)) {
    job_seconds <- 0.5 * old$job_seconds + 0.5 * job_seconds
  }
  if (!is.null(startup_seconds) && is.finite(startup_seconds) &&
      startup_seconds >= 0) {
    if (is.list(old) && is.finite(old$startup_seconds))
      startup_seconds <- 0.5 * old$startup_seconds + 0.5 * startup_seconds
  } else {
    startup_seconds <- if (is.list(old) && is.finite(old$startup_seconds))
      old$startup_seconds else 0
  }
  .emc_rebuild_queue$costs[[id]] <- list(job_seconds = job_seconds,
                                        startup_seconds = startup_seconds,
                                        updated = proc.time()[["elapsed"]])
  invisible(NULL)
}

.emc_rebuild_serial <- function(X, FUN, order) {
  .emc_with_preserved_rng({
    n <- length(X)
    out <- vector("list", n)
    elapsed <- numeric(n)
    for (i in as.integer(order)) {
      started <- proc.time()[["elapsed"]]
      out[i] <- list(FUN(X[[i]]))
      elapsed[i] <- proc.time()[["elapsed"]] - started
    }
    names(out) <- names(X)
    list(values = out, job_seconds = elapsed)
  })
}

.emc_rebuild_timed_job <- function(i, X, FUN) {
  started <- proc.time()[["elapsed"]]
  value <- FUN(X[[i]])
  list(index = i, value = value,
       seconds = proc.time()[["elapsed"]] - started)
}

.emc_rebuild_lapply <- function(X, FUN, n_cores, order = seq_along(X),
                                cost_key = NULL) {
  n <- length(X)
  if (!n) return(vector("list", 0L))
  cores <- max(1L, min(n, as.integer(n_cores)[1L]))
  stopifnot(length(order) == n,
            identical(sort(as.integer(order)), seq_len(n)))

  prediction <- .emc_rebuild_prediction(cost_key, n, cores)
  margin <- getOption("emc2.rebuild_serial_margin", 0.90)
  use_serial <- !is.null(prediction) &&
    !isTRUE(.emc_rebuild_queue$busy) &&
    is.finite(margin) && margin > 0 && margin < 1 &&
    prediction$serial < margin * prediction$parallel

  if (use_serial) {
    result <- .emc_rebuild_serial(X, FUN, order)
    .emc_rebuild_update_cost(cost_key, cores, result$job_seconds)
    .emc_rebuild_queue$last_decision <- list(
      key = cost_key, route = "serial", serial_seconds = prediction$serial,
      parallel_seconds = prediction$parallel, n = n, cores = cores)
    return(result$values)
  }

  if (is.null(cost_key)) {
    if (.Platform$OS.type != "unix" || cores <= 1L || n <= 2L * cores ||
        isTRUE(.emc_rebuild_queue$busy))
      return(auto_mclapply(X, FUN, mc.cores = cores, mc.set.seed = FALSE))
  } else if (.Platform$OS.type != "unix" || cores <= 1L ||
             n <= 2L * cores || isTRUE(.emc_rebuild_queue$busy)) {
    started <- proc.time()[["elapsed"]]
    items <- auto_mclapply(as.integer(order), function(i)
      .emc_rebuild_timed_job(i, X, FUN), mc.cores = cores,
      mc.set.seed = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started
    out <- vector("list", n)
    job_seconds <- numeric(n)
    for (item in items) {
      out[item$index] <- list(item$value)
      job_seconds[item$index] <- item$seconds
    }
    names(out) <- names(X)
    startup <- if (cores > 1L && .Platform$OS.type == "unix")
      max(0, elapsed - .emc_rebuild_makespan(job_seconds, cores)) else NULL
    .emc_rebuild_update_cost(cost_key, cores, job_seconds, startup)
    .emc_rebuild_queue$last_decision <- list(
      key = cost_key, route = if (cores > 1L) "parallel" else "serial",
      serial_seconds = sum(job_seconds),
      parallel_seconds = .emc_rebuild_makespan(job_seconds, cores) +
        (startup %||% 0), n = n, cores = cores)
    return(out)
  }

  .emc_rebuild_queue$busy <- TRUE
  on.exit(.emc_rebuild_queue$busy <- FALSE, add = TRUE)
  map <- .emc_rebuild_queue$map
  pid <- Sys.getpid()
  if (is.null(map) || !identical(.emc_rebuild_queue$pid, pid)) {
    path <- tempfile("emc-rebuild-claim-")
    on.exit(unlink(path), add = TRUE)
    map <- emc_claim_open(path)
    unlink(path)
    .emc_rebuild_queue$map <- map
    .emc_rebuild_queue$pid <- pid
  }
  # A fresh key per round: a fork left over from an earlier round still maps
  # the counter, and must not claim this round's jobs.
  key <- .emc_rebuild_queue$key <- (.emc_rebuild_queue$key %||% 0L) %% 1000000L + 1L
  emc_claim_reset(map, key)
  started <- proc.time()[["elapsed"]]
  parts <- auto_mclapply(seq_len(cores), function(worker) {
    indices <- integer(n); values <- vector("list", n)
    job_seconds <- numeric(n); used <- 0L
    repeat {
      claim <- emc_claim_next(map, key, n)
      if (claim < 0L) break
      i <- as.integer(order[claim + 1L])
      item <- .emc_rebuild_timed_job(i, X, FUN)
      used <- used + 1L
      indices[used] <- i
      values[used] <- list(item$value)
      job_seconds[used] <- item$seconds
    }
    list(indices = indices[seq_len(used)], values = values[seq_len(used)],
         job_seconds = job_seconds[seq_len(used)],
         worker_seconds = sum(job_seconds[seq_len(used)]))
  }, mc.cores = cores, mc.set.seed = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  out <- vector("list", n)
  job_seconds <- numeric(n)
  worker_seconds <- numeric(length(parts))
  for (w in seq_along(parts)) {
    part <- parts[[w]]
    out[part$indices] <- part$values
    job_seconds[part$indices] <- part$job_seconds
    worker_seconds[w] <- part$worker_seconds
  }
  names(out) <- names(X)
  startup <- max(0, elapsed - max(worker_seconds))
  .emc_rebuild_update_cost(cost_key, cores, job_seconds, startup)
  .emc_rebuild_queue$last_decision <- list(
    key = cost_key, route = "parallel", serial_seconds = sum(job_seconds),
    parallel_seconds = max(worker_seconds) + startup, n = n, cores = cores)
  out
}
