# One profiling record schema for the sampler.
#
# Every low-level timing number the sampler produces goes through this file.
# Before it existed the worker pool built an ad-hoc list, `run_stage()` copied
# that list field-by-field into a data frame, and each benchmark script printed
# its own selection with its own labels.  Three places had to agree on a name
# and none of them declared it, which is how `worker_max` came to be documented
# in `WorkingTests/bench_worker_pool.R` as "slowest worker's own work" while it
# actually held the slowest *subject's* CPU time.  With several subjects per
# worker those differ by roughly the number of subjects per worker, and the
# difference was being read as transport overhead.
#
# The schema below is the single declaration.  `.emc_profile_row()` builds rows
# from it, `.emc_profile_report()` prints them, and both benchmark scripts call
# that one reporting path.  Adding a measurement means adding one entry here.
#
# Reserved names for later work, so that they are not invented twice with
# different spellings: `block_proposals`, `block_eff_proposals`,
# `block_check_progress`, `block_save`, `block_concat` (per-block costs, absent
# from per-iteration timing today), and `kernel_reuse_*` (kernel preparation
# reuse rates).  Add them here when the commit that populates them lands.

# --- schema -----------------------------------------------------------------

# `parent` places a field in the time breakdown: fields whose parent is "total"
# are components of the iteration, and fields naming another field are
# components of *it*.  The report indents accordingly and takes each percentage
# against its own parent, so nested numbers do not silently sum past 100%.
.emc_profile_field <- function(name, type, group, unit = "", parent = NA_character_,
                               desc = "") {
  list(name = name, type = type, group = group, unit = unit,
       parent = parent, desc = desc)
}

.emc_profile_schema <- local({
  reject_fields <- function(source) {
    lapply(c("numerical", "infrastructure", "programming", "unknown"),
           function(cls) .emc_profile_field(
             sprintf("reject_%s_%s", source, substr(cls, 1L, 4L)),
             "integer", "health", "count",
             desc = sprintf("%s rejections classified as %s", source, cls)))
  }
  fields <- c(
    list(
      # -- identity ---------------------------------------------------------
      .emc_profile_field("iteration", "integer", "id", "index",
                         desc = "sample store index written this iteration"),
      .emc_profile_field("stage", "character", "id", "",
                         desc = "sampling stage"),
      .emc_profile_field("subjects", "integer", "id", "count",
                         desc = "subjects updated"),
      .emc_profile_field("workers", "integer", "id", "count",
                         desc = "pool width"),
      .emc_profile_field("workers_active", "integer", "id", "count",
                         desc = "workers given at least one subject"),
      .emc_profile_field("particles_max", "integer", "id", "count",
                         desc = "largest effective particle count over subjects"),
      .emc_profile_field("particles_sum", "integer", "id", "count",
                         desc = "effective particle counts summed over subjects"),

      # -- time -------------------------------------------------------------
      .emc_profile_field("total", "numeric", "time", "s",
                         desc = "whole iteration, master elapsed"),
      .emc_profile_field("gibbs", "numeric", "time", "s", "total",
                         desc = "group-level Gibbs update"),
      .emc_profile_field("group_cache", "numeric", "time", "s", "total",
                         desc = "group covariance factorisation"),
      .emc_profile_field("recycle", "numeric", "time", "s", "total",
                         desc = "worker re-fork"),
      .emc_profile_field("grow", "numeric", "time", "s", "total",
                         desc = "worker pool growth"),
      .emc_profile_field("particle", "numeric", "time", "s", "total",
                         desc = "particle step, master elapsed"),
      .emc_profile_field("fill", "numeric", "time", "s", "total",
                         desc = "writing the iteration into the sample store"),
      # Not a component of any iteration: pool start-up happens once, before
      # the loop, and is reported on the first row so that dropping that row
      # (which every benchmark does) does not also hide what it cost.
      .emc_profile_field("startup", "numeric", "time", "s",
                         desc = "worker pool start-up, charged to the first iteration"),
      .emc_profile_field("shared_serialize", "numeric", "time", "s", "particle",
                         desc = "serialising and publishing the shared group draw"),
      .emc_profile_field("request_send", "numeric", "time", "s", "particle",
                         desc = "writing per-worker requests"),
      .emc_profile_field("response_wait", "numeric", "time", "s", "particle",
                         desc = "waiting for and reading replies"),

      # -- work distribution ------------------------------------------------
      # `subject_*` are per-subject CPU times; `worker_*` aggregate those over
      # each worker's assigned partition.  Keeping both is the point: the ratio
      # is the number that says whether a partition is balanced, and reporting
      # only the per-subject maximum understates worker cost.
      .emc_profile_field("subject_cpu_max", "numeric", "work", "s",
                         desc = "largest single subject's CPU time"),
      .emc_profile_field("subject_cpu_sum", "numeric", "work", "s",
                         desc = "subject CPU time summed over subjects"),
      .emc_profile_field("worker_cpu_max", "numeric", "work", "s",
                         desc = "largest CPU workload assigned to one worker"),
      .emc_profile_field("worker_cpu_sum", "numeric", "work", "s",
                         desc = "worker-reported CPU time summed over workers"),
      .emc_profile_field("worker_elapsed_max", "numeric", "work", "s",
                         desc = "longest worker-reported elapsed time"),
      .emc_profile_field("worker_elapsed_sum", "numeric", "work", "s",
                         desc = "worker-reported elapsed time summed over workers"),
      .emc_profile_field("dispatch_first", "numeric", "work", "s",
                         desc = "first request written, relative to iteration start"),
      .emc_profile_field("dispatch_last", "numeric", "work", "s",
                         desc = "last request written, relative to iteration start"),
      .emc_profile_field("finish_first", "numeric", "work", "s",
                         desc = "first reply stored, relative to iteration start"),
      .emc_profile_field("finish_last", "numeric", "work", "s",
                         desc = "last reply stored, relative to iteration start"),
      .emc_profile_field("queue_delay_max", "numeric", "work", "s",
                         desc = "longest wait between a request being written and picked up"),
      .emc_profile_field("queue_delay_mean", "numeric", "work", "s",
                         desc = "mean wait between a request being written and picked up"),

      # -- wire (opt-in: each of these walks and serialises a message) -------
      .emc_profile_field("shared_bytes", "numeric", "wire", "B",
                         desc = "serialised shared group draw"),
      .emc_profile_field("private_bytes", "numeric", "wire", "B",
                         desc = "serialised per-worker request payloads"),
      .emc_profile_field("wire_bytes", "numeric", "wire", "B",
                         desc = "total bytes crossing the boundary this iteration"),

      # -- memory -----------------------------------------------------------
      .emc_profile_field("private_memory", "numeric", "memory", "B",
                         desc = "process-tree private memory (PSS where available)"),
      .emc_profile_field("private_memory_method", "character", "memory", "",
                         desc = "how private_memory was obtained"),

      # -- per-block costs ---------------------------------------------------
      # C16.  `run_emc` does work between blocks that never appeared in
      # per-iteration timing: building chain and efficient proposals,
      # convergence diagnostics, the checkpoint write and the history join.  A
      # block is not an iteration, so these carry their own total rather than
      # nesting under one.
      .emc_profile_field("block", "integer", "block", "index",
                         desc = "block number within the stage"),
      .emc_profile_field("block_iterations", "integer", "block", "count",
                         desc = "iterations requested in this block"),
      .emc_profile_field("block_chains", "integer", "block", "count",
                         desc = "chains run in this block"),
      .emc_profile_field("block_total", "numeric", "block", "s",
                         desc = "whole block, master elapsed"),
      .emc_profile_field("block_prepare", "numeric", "block", "s", "block_total",
                         desc = "trimming, resetting and subsetting before sampling"),
      .emc_profile_field("block_proposals", "numeric", "block", "s", "block_total",
                         desc = "create_chain_proposals and sub-blocking"),
      .emc_profile_field("block_eff_proposals", "numeric", "block", "s", "block_total",
                         desc = "create_eff_proposals"),
      .emc_profile_field("block_sample", "numeric", "block", "s", "block_total",
                         desc = "the chains' sampling call, start-up included"),
      .emc_profile_field("block_concat", "numeric", "block", "s", "block_total",
                         desc = "thinning and joining the block into the fit"),
      .emc_profile_field("block_check_progress", "numeric", "block", "s", "block_total",
                         desc = "convergence diagnostics"),
      .emc_profile_field("block_save", "numeric", "block", "s", "block_total",
                         desc = "checkpoint write"),

      # -- particle adaptation ----------------------------------------------
      # What `update_pm_settings` adapts the particle count from: the effective
      # size of one proposal cloud, NOT the chain's effective sample size. The
      # two are different questions and the audit's experiment is whether the
      # first tracks the second.
      .emc_profile_field("particle_weight_ess", "numeric", "work", "count",
                         desc = "importance-weight ESS of the proposal cloud, mean over subjects"),

      # -- kernel reuse (opt-in: emc_kernel_stats(TRUE)) ---------------------
      # How much of a model kernel's arithmetic is recomputed per trial that a
      # cell-resolution version would compute once.  `rows / cells` is the
      # reuse available; near 1 means there is nothing to hoist.
      .emc_profile_field("kernel_rows", "numeric", "kernel", "count",
                         desc = "per-trial evaluations of the reusable subexpression"),
      .emc_profile_field("kernel_cells", "numeric", "kernel", "count",
                         desc = "evaluations the same subexpression needs at cell resolution"),
      .emc_profile_field("kernel_seconds", "numeric", "kernel", "s",
                         desc = "wall time inside the model kernel"),

      # -- health -----------------------------------------------------------
      .emc_profile_field("fallback_subjects", "integer", "health", "count",
                         desc = "subjects recomputed in the master after a failure"),
      .emc_profile_field("fallback_workers", "integer", "health", "count",
                         desc = "worker shares that failed or were never sent"),
      .emc_profile_field("degraded", "integer", "health", "count",
                         desc = "1 if the pool was marked dead during this iteration")
    ),
    reject_fields("particle"),
    reject_fields("gibbs")
  )
  names(fields) <- vapply(fields, `[[`, character(1), "name")
  fields
})

.emc_profile_default <- function(type) {
  switch(type,
         integer = NA_integer_,
         character = NA_character_,
         NA_real_)
}

# --- options ----------------------------------------------------------------

.emc_profile_enabled <- function() {
  isTRUE(getOption("emc2.sampler_profile", FALSE))
}

# Serialising every request to measure its size costs one extra object walk per
# worker per iteration, and the tree memory read touches /proc for every
# process.  Both stay behind their own switch so an ordinary profile run
# measures the sampler rather than the profiler.
.emc_profile_expensive <- function() {
  .emc_profile_enabled() &&
    isTRUE(getOption("emc2.sampler_profile_expensive", FALSE))
}

# --- row construction -------------------------------------------------------

# Build one profile row.  Unknown field names are an error rather than a silent
# extra column: the whole purpose of the schema is that a typo in a benchmark
# cannot quietly create `worker_max` again.
.emc_profile_row <- function(...) {
  vals <- list(...)
  if (length(vals) && is.null(names(vals))) stop("profile fields must be named")
  unknown <- setdiff(names(vals), names(.emc_profile_schema))
  if (length(unknown)) {
    stop("unknown profile field(s): ", paste(unknown, collapse = ", "),
         ". Declare them in .emc_profile_schema.")
  }
  row <- lapply(.emc_profile_schema, function(f) {
    v <- vals[[f$name]]
    if (is.null(v) || length(v) != 1L) return(.emc_profile_default(f$type))
    switch(f$type,
           integer = if (is.na(v)) NA_integer_ else as.integer(v),
           character = as.character(v),
           as.numeric(v))
  })
  structure(row, class = "data.frame", row.names = 1L)
}

.emc_profile_bind <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    empty <- .emc_profile_row()
    return(empty[0L, , drop = FALSE])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- tree-aware private memory ----------------------------------------------

# Summing RSS over a forked tree counts every shared page once per process, so
# a pool of eight workers reports roughly eight times the data it holds.  PSS
# divides each page by the number of processes mapping it, which is the number
# that answers "how much memory does this fit need".  Linux exposes it in
# smaps_rollup; elsewhere we report RSS and say so, rather than reporting a
# number whose meaning depends on the platform without labelling it.
# Every process a fit is actually using, root included.
#
# Two mechanisms have to be followed, because the pool uses both.  Chain
# processes are forked, so they are ordinary descendants and a walk over
# /proc's parent links finds them.  Worker pool *templates* are not: the
# default backend starts them with `system2(wait = FALSE)`, whose shell exits
# immediately, so init adopts the template and the fit's own subtree no longer
# contains it or any of the workers it forks.  On an eight-subject fit that is
# most of the processes and most of the memory.
#
# The pool publishes the ownership we need: `.emc_wpool_template_start` puts
# `EMC2-worker-pool[<dir>|master=<pid>]` in the template's command line, and
# the workers inherit it by forking.  `adopt` claims exactly the templates
# whose recorded master is already in the tree, so this is reading an
# ownership record rather than guessing from a process name.
#
# Linux only; elsewhere the caller's own pid list is all there is.
.emc_profile_tree_pids <- function(root = Sys.getpid(), adopt = TRUE) {
  root <- as.integer(root)
  root <- root[!is.na(root)]
  if (!length(root) || !dir.exists("/proc")) return(root)
  all <- suppressWarnings(as.integer(list.files("/proc")))
  all <- all[!is.na(all)]
  if (!length(all)) return(root)

  # A process can exit between listing /proc and reading its stat file.  That
  # is the normal case, not an exceptional one, so it must be silent: the
  # connection warning would otherwise be emitted once per departed worker.
  stat <- lapply(all, function(pid) {
    tryCatch(suppressWarnings(
      readLines(file.path("/proc", pid, "stat"), n = 1L, warn = FALSE)),
      error = function(e) character())
  })
  # The comm field is parenthesised and may itself contain spaces and
  # brackets, so split after the *last* ')' rather than on whitespace.
  parent <- vapply(stat, function(line) {
    if (!length(line)) return(NA_integer_)
    fields <- strsplit(trimws(substring(line, regexpr("\\)[^)]*$", line) + 1L)),
                       "[ ]+")[[1L]]
    if (length(fields) < 2L) return(NA_integer_)
    suppressWarnings(as.integer(fields[2L]))
  }, integer(1))
  comm <- vapply(stat, function(line) {
    if (!length(line)) return(NA_character_)
    open <- regexpr("(", line, fixed = TRUE)
    close <- regexpr("\\)[^)]*$", line)
    if (open < 1L || close <= open) return(NA_character_)
    substring(line, open + 1L, close - 1L)
  }, character(1))

  # Bounded: each pass adds at least one process or stops, and the process
  # table is finite, so this cannot spin on a cycle in a corrupted /proc.
  descend <- function(out) {
    frontier <- out
    for (depth in seq_along(all)) {
      kids <- setdiff(all[!is.na(parent) & parent %in% frontier], out)
      if (!length(kids)) break
      out <- c(out, kids)
      frontier <- kids
    }
    out
  }
  out <- descend(root)
  if (!isTRUE(adopt)) return(unique(out))

  # Only R processes can be pool templates, and reading a command line for
  # every process on a busy machine costs more than the walk itself.
  cand <- all[!is.na(comm) & comm %in% c("R", "Rscript", "exec") & !(all %in% out)]
  if (!length(cand)) return(unique(out))
  owner <- vapply(cand, function(pid) {
    raw <- tryCatch(suppressWarnings(
      readBin(file.path("/proc", pid, "cmdline"), "raw", 16384L)),
      error = function(e) raw())
    if (!length(raw)) return(NA_integer_)
    txt <- rawToChar(raw[raw != as.raw(0L)])
    hit <- regmatches(txt, regexpr("EMC2-worker-pool\\[[^]]*\\|master=[0-9]+\\]", txt))
    if (!length(hit)) return(NA_integer_)
    suppressWarnings(as.integer(sub(".*master=([0-9]+)\\]$", "\\1", hit[[1L]])))
  }, integer(1))
  adopted <- cand[!is.na(owner) & owner %in% out]
  if (!length(adopted)) return(unique(out))
  unique(descend(c(out, adopted)))
}

# `include_self = FALSE` measures somebody else's tree.  The sampler asking
# about its own workers is part of the answer and keeps the default; a
# benchmark that forks the whole fit and watches it from outside is not, and
# adding the observer's pages to the reading would inflate every number it
# prints by whatever the observer happens to be holding.
.emc_profile_memory <- function(pids = NULL, tree = FALSE, include_self = TRUE,
                                adopt = TRUE) {
  pids <- as.integer(pids[!is.na(pids)])
  if (isTRUE(include_self)) pids <- c(Sys.getpid(), pids)
  pids <- unique(pids)
  # One /proc scan for the whole root set, not one per root: the walk already
  # takes a vector, and scanning the process table repeatedly is most of the
  # cost of polling a live tree.
  if (isTRUE(tree)) pids <- .emc_profile_tree_pids(pids, adopt = adopt)
  if (!length(pids)) return(list(bytes = NA_real_, method = NA_character_))
  rollup <- file.path("/proc", pids, "smaps_rollup")
  have <- file.exists(rollup)
  if (any(have)) {
    kb <- vapply(rollup[have], function(f) {
      # `file.exists` above and this read are not one operation: a worker can
      # exit in between, and while the error is caught the connection warning
      # is not, so a tree poll would print one line per departed process.
      lines <- tryCatch(suppressWarnings(readLines(f, warn = FALSE)),
                        error = function(e) character())
      hit <- grep("^Pss:", lines, value = TRUE)
      if (!length(hit)) return(NA_real_)
      as.numeric(sub("^Pss:\\s*([0-9]+).*$", "\\1", hit[1L]))
    }, numeric(1))
    if (any(!is.na(kb))) {
      # A pid that has already exited is not a gap in the reading: it holds no
      # memory, and leaving it out is right.  A process that is still there and
      # whose rollup could not be read is a gap, and the sum is then short by
      # an unknown amount, which is worth saying.
      short <- sum(is.na(kb))
      return(list(bytes = sum(kb, na.rm = TRUE) * 1024,
                  method = if (short) "pss-partial" else "pss"))
    }
  }
  status <- file.path("/proc", pids, "status")
  have <- file.exists(status)
  if (any(have)) {
    kb <- vapply(status[have], function(f) {
      lines <- tryCatch(suppressWarnings(readLines(f, warn = FALSE)),
                        error = function(e) character())
      hit <- grep("^VmRSS:", lines, value = TRUE)
      if (!length(hit)) return(NA_real_)
      as.numeric(sub("^VmRSS:\\s*([0-9]+).*$", "\\1", hit[1L]))
    }, numeric(1))
    if (any(!is.na(kb))) {
      short <- sum(is.na(kb))
      return(list(bytes = sum(kb, na.rm = TRUE) * 1024,
                  method = if (short) "rss-sum-partial" else "rss-sum"))
    }
  }
  if (.Platform$OS.type == "unix") {
    out <- tryCatch(
      suppressWarnings(system2("ps", c("-o", "rss=", "-p",
                                       paste(pids, collapse = ",")),
                               stdout = TRUE, stderr = FALSE)),
      error = function(e) character())
    kb <- suppressWarnings(as.numeric(trimws(out)))
    if (length(kb) && any(!is.na(kb))) {
      return(list(bytes = sum(kb, na.rm = TRUE) * 1024, method = "rss-sum"))
    }
  }
  list(bytes = NA_real_, method = NA_character_)
}

# --- rejection counters -----------------------------------------------------

# Counting only.  Nothing here changes what the sampler does with a failure:
# `safe_new_particle()` still repeats the previous state, which is what
# preserves the PMwG invariant under numerical rejection.  Deciding which
# classified failures should instead propagate is a separate policy question.
.emc_profile_state <- new.env(parent = emptyenv())

.EMC_REJECT_CLASSES <- c("numerical", "infrastructure", "programming", "unknown")

.emc_reject_reset <- function(source = NULL) {
  if (is.null(source)) {
    .emc_profile_state$rejects <- NULL
    .emc_profile_state$reject_msgs <- NULL
    return(invisible(NULL))
  }
  counts <- .emc_profile_state$rejects
  if (!is.null(counts)) counts[[source]] <- NULL
  .emc_profile_state$rejects <- counts
  invisible(NULL)
}

# The classification is a diagnostic, not a contract.  It reads the condition's
# class and message, so a numerical failure that R reports as a plain
# `simpleError` with an unfamiliar message lands in "unknown" rather than being
# mislabelled.  "unknown" staying large is itself the useful signal: it means
# the classifier needs another pattern before any policy could act on it.
.emc_classify_failure <- function(cond) {
  if (is.null(cond)) return("unknown")
  msg <- tryCatch(conditionMessage(cond), error = function(e) "")
  if (!length(msg)) msg <- ""
  msg <- tolower(paste(msg, collapse = " "))
  cls <- class(cond)
  if (any(c("emc_transport_error", "emc_pool_error") %in% cls)) {
    return("infrastructure")
  }
  numerical <- paste(
    "not positive definite", "leading minor", "chol", "singular",
    "non-finite", "not finite", "infinite", "nan", "na/nan",
    "system is exactly singular", "numerically", "cannot be inverted",
    "missing value where", "0 \\(non-na\\)", sep = "|")
  infrastructure <- paste(
    "cannot open", "connection", "fifo", "pipe", "no such file",
    "unserialize", "serialize", "cannot allocate", "ignoring sigpipe",
    "worker", "socket", "timed out", "unable to fork", sep = "|")
  programming <- paste(
    "could not find function", "object '", "unused argument",
    "argument \".*\" is missing", "subscript out of bounds",
    "non-conformable", "incorrect number of dimensions",
    "invalid 'type'", "arguments imply differing number", sep = "|")
  # Infrastructure first: an allocation failure inside a numerical routine is
  # still an allocation failure, and treating it as ordinary numerical
  # rejection is exactly the mistake that hides a broken pool as slow mixing.
  if (grepl(infrastructure, msg)) return("infrastructure")
  if (grepl(programming, msg)) return("programming")
  if (grepl(numerical, msg)) return("numerical")
  "unknown"
}

.emc_reject_record <- function(cond, source = "particle") {
  cls <- .emc_classify_failure(cond)
  counts <- .emc_profile_state$rejects
  if (is.null(counts)) counts <- list()
  cur <- counts[[source]]
  if (is.null(cur)) {
    cur <- stats::setNames(integer(length(.EMC_REJECT_CLASSES)), .EMC_REJECT_CLASSES)
  }
  cur[[cls]] <- cur[[cls]] + 1L
  counts[[source]] <- cur
  .emc_profile_state$rejects <- counts
  # Keep the first message per (source, class) so a nonzero count can be chased
  # without re-running the fit with a debugger attached.
  key <- paste(source, cls, sep = ".")
  msgs <- .emc_profile_state$reject_msgs
  if (is.null(msgs)) msgs <- list()
  if (is.null(msgs[[key]])) {
    msgs[[key]] <- tryCatch(conditionMessage(cond), error = function(e) "")
    .emc_profile_state$reject_msgs <- msgs
  }
  # Counting is still all this does to the sampler.  Announcing is the policy's
  # job and lives in R/failure_policy.R; it is called from here so that no
  # caller can count a failure without the policy having seen it.
  .emc_failure_announce(cls, source, cond)
  invisible(cls)
}

.emc_reject_counts <- function(source = "particle") {
  cur <- .emc_profile_state$rejects[[source]]
  if (is.null(cur)) {
    stats::setNames(integer(length(.EMC_REJECT_CLASSES)), .EMC_REJECT_CLASSES)
  } else cur
}

# Counters live in the process that caught the error, so a worker's counts must
# travel back in its reply.  Snapshot before the work, difference after.
# Every source's totals at once, in the shape `.emc_failure_summary()` reads.
# Per-block records, collected in the master.  `run_emc` runs there, so unlike
# the per-iteration record -- which is built inside a forked chain and has to
# travel back -- these can simply be accumulated.
.emc_block_reset <- function() {
  .emc_profile_state$blocks <- NULL
  invisible(NULL)
}

.emc_block_record <- function(...) {
  if (!.emc_profile_enabled()) return(invisible(NULL))
  rows <- .emc_profile_state$blocks
  if (is.null(rows)) rows <- list()
  rows[[length(rows) + 1L]] <- .emc_profile_row(...)
  .emc_profile_state$blocks <- rows
  invisible(NULL)
}

.emc_block_profile <- function() {
  rows <- .emc_profile_state$blocks
  if (is.null(rows) || !length(rows)) return(NULL)
  .emc_profile_bind(rows)
}

# Kernel reuse totals, summed over whatever models ran.  NA when the C++
# counters are not compiled in or instrumentation was never switched on, so a
# profile from an ordinary run carries no kernel columns rather than zeros that
# would read as "measured, and there is no reuse".
.emc_kernel_totals <- function() {
  # Every worker request calls this twice, so the ordinary case -- measurement
  # off -- must not pay for the read: building the counter table as a data
  # frame costs ~200 us, the flag ~1 us.  A process that never switched it on
  # has recorded nothing, so this is the same answer, reached directly.
  if (!isTRUE(tryCatch(emc_kernel_stats(), error = function(e) FALSE))) {
    return(c(rows = NA_real_, cells = NA_real_, seconds = NA_real_))
  }
  st <- tryCatch(emc_kernel_stats_read(), error = function(e) NULL)
  if (is.null(st) || !nrow(st)) {
    return(c(rows = NA_real_, cells = NA_real_, seconds = NA_real_))
  }
  c(rows = sum(st$rows), cells = sum(st$cells), seconds = sum(st$seconds))
}

.emc_kernel_delta <- function(before) {
  now <- .emc_kernel_totals()
  if (anyNA(now) || anyNA(before)) {
    return(c(rows = NA_real_, cells = NA_real_, seconds = NA_real_))
  }
  now - before
}

.emc_reject_counts_all <- function() {
  counts <- .emc_profile_state$rejects
  if (is.null(counts)) list() else counts
}

.emc_reject_delta <- function(before, source = "particle") {
  .emc_reject_counts(source) - before
}

.emc_reject_messages <- function() {
  msgs <- .emc_profile_state$reject_msgs
  if (is.null(msgs)) list() else msgs
}

# Flatten a class-count vector into the schema's per-class columns.
.emc_reject_row_fields <- function(counts, source) {
  if (is.null(counts)) return(list())
  out <- as.list(as.integer(counts))
  names(out) <- sprintf("reject_%s_%s", source,
                        substr(names(counts), 1L, 4L))
  out
}

# --- reporting --------------------------------------------------------------

.emc_profile_fmt_bytes <- function(x) {
  if (!length(x) || all(is.na(x))) return("       n/a")
  x <- mean(x, na.rm = TRUE)
  units <- c("B", "kB", "MB", "GB", "TB")
  i <- 1L
  while (x >= 1024 && i < length(units)) { x <- x / 1024; i <- i + 1L }
  sprintf("%7.2f %s", x, units[i])
}

.emc_profile_mean <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# The per-block section, split out so a run with blocks but no per-iteration
# record -- which is every multi-block run -- still reports.
.emc_block_report <- function(blocks, elapsed = NULL, con = stdout()) {
  say <- function(...) cat(..., sep = "", file = con)
  say(sprintf("\nper-block costs (%d block%s)\n", nrow(blocks),
              if (nrow(blocks) > 1L) "s" else ""))
  btot <- .emc_profile_mean(blocks$block_total)
  say(sprintf("  %-22s %8.2f ms  (mean over blocks)\n", "block_total",
              1000 * btot))
  bf <- Filter(function(f) f$group == "block" && identical(f$parent, "block_total"),
               .emc_profile_schema)
  for (f in bf) {
    v <- .emc_profile_mean(blocks[[f$name]])
    if (is.na(v)) next
    pct <- if (is.na(btot) || btot <= 0) NA_real_ else 100 * v / btot
    say(sprintf("    %-20s %8.2f ms  %s\n", f$name, 1000 * v,
                if (is.na(pct)) "" else sprintf("%5.1f%%", pct)))
  }
  total_blocks <- sum(blocks$block_total, na.rm = TRUE)
  if (!is.null(elapsed) && elapsed > 0) {
    say(sprintf("  %-22s %8.2f s   %5.1f%% of the run\n", "all blocks",
                total_blocks, 100 * total_blocks / elapsed))
  }
  invisible(blocks)
}

#' Print a sampler profile
#'
#' The single reporting path for the profiling record described in
#' `.emc_profile_schema`.  Both `WorkingTests/bench_worker_pool.R` and
#' `WorkingTests/bench_large_fit_sampling.R` call this, so a field only has to
#' be formatted once and the two benchmarks cannot drift apart in what they
#' call a number.
#'
#' @param profile Data frame produced by the sampler with
#'   `options(emc2.sampler_profile = TRUE)`.
#' @param elapsed Optional wall clock for the whole run, in seconds.
#' @param drop_first Drop the first iteration, which pays pool start-up.
#' @param blocks Optional per-block record from `.emc_block_profile()`.
#' @param con Connection to write to.
#' @noRd
.emc_profile_report <- function(profile, elapsed = NULL, drop_first = TRUE,
                                con = stdout(), blocks = NULL) {
  say <- function(...) cat(..., sep = "", file = con)
  if (is.null(profile) || !nrow(profile)) {
    # Blocks are recorded in the master and survive; the per-iteration record is
    # built inside a forked chain and does not outlive `concat_emc`.  Having one
    # without the other is the ordinary case for a multi-block run, not an
    # error, so report what there is.
    if (!is.null(blocks) && nrow(blocks)) {
      .emc_block_report(blocks, elapsed = elapsed, con = con)
    } else {
      say("no profile recorded\n")
    }
    return(invisible(NULL))
  }
  # `full` keeps every iteration; `profile` drops the first for the per-iteration
  # averages, because iteration 1 pays pool start-up and would distort them.
  # Anything that is *not* a per-iteration average -- start-up itself, the
  # sampled memory readings, and the health totals -- is taken from `full`, or
  # dropping the first row would also hide the cost it was dropped for.
  full <- profile
  if (drop_first && nrow(profile) > 1L) {
    profile <- profile[-1L, , drop = FALSE]
  }
  m <- function(nm) .emc_profile_mean(profile[[nm]])
  tot <- m("total")

  say("\n")
  say(sprintf("iterations %d", nrow(profile)))
  if (!is.na(m("subjects"))) say(sprintf("   subjects %d", as.integer(m("subjects"))))
  if (!is.na(m("workers"))) {
    say(sprintf("   workers %d", as.integer(max(profile$workers, na.rm = TRUE))))
    act <- m("workers_active")
    if (!is.na(act)) say(sprintf(" (%.1f active)", act))
  }
  say("\n")
  if (!is.null(elapsed)) say(sprintf("wall %.2f s", elapsed), "   ")
  say(sprintf("iteration %.2f ms\n\n", 1000 * tot))

  # Time breakdown, driven by the schema rather than a hand-kept list, so a new
  # timing field appears here as soon as it is declared and populated.
  time_fields <- Filter(function(f) f$group == "time" && f$name != "total",
                        .emc_profile_schema)
  emit <- function(fields, parent, indent) {
    for (f in fields) {
      if (!identical(f$parent, parent)) next
      v <- m(f$name)
      if (is.na(v)) next
      base <- if (identical(parent, "total")) tot else m(parent)
      pct <- if (is.na(base) || base <= 0) NA_real_ else 100 * v / base
      say(sprintf("%s%-*s %8.2f ms  %s\n", strrep(" ", indent),
                  20 - indent, f$name, 1000 * v,
                  if (is.na(pct)) "" else sprintf("%5.1f%%", pct)))
      emit(fields, f$name, indent + 2L)
    }
  }
  emit(time_fields, "total", 2L)
  # Time fields that are not components of an iteration (pool start-up) are
  # reported separately rather than folded into a percentage of one iteration.
  oneoff <- Filter(function(f) is.na(f$parent), time_fields)
  for (f in oneoff) {
    v <- sum(full[[f$name]], na.rm = TRUE)
    if (!is.finite(v) || v <= 0) next
    say(sprintf("  %-20s %8.2f ms  (%s)\n", f$name, 1000 * v, f$desc))
  }

  say("\nwork distribution\n")
  # These two lines are the correction this schema exists for.  The first is
  # the slowest single subject; the second is the slowest worker's whole
  # assigned partition.  Comparing `response_wait` against the first one makes
  # every subject after the first in a partition look like transport cost.
  say(sprintf("  %-20s %8.2f ms  (slowest single subject)\n",
              "subject_cpu_max", 1000 * m("subject_cpu_max")))
  say(sprintf("  %-20s %8.2f ms  (slowest worker's assigned partition)\n",
              "worker_cpu_max", 1000 * m("worker_cpu_max")))
  say(sprintf("  %-20s %8.2f ms  (summed over subjects)\n",
              "subject_cpu_sum", 1000 * m("subject_cpu_sum")))
  we <- m("worker_elapsed_max")
  if (!is.na(we)) {
    say(sprintf("  %-20s %8.2f ms  (worker-reported)\n",
                "worker_elapsed_max", 1000 * we))
  }
  wait <- m("response_wait")
  wmax <- m("worker_cpu_max")
  if (!is.na(wait) && !is.na(wmax)) {
    say(sprintf("  %-20s %8.2f ms  (response_wait - worker_cpu_max)\n",
                "unexplained", 1000 * (wait - wmax)))
  }
  if (!is.na(m("queue_delay_max"))) {
    say(sprintf("  %-20s %8.2f ms  max, %.2f ms mean\n", "queue_delay",
                1000 * m("queue_delay_max"), 1000 * m("queue_delay_mean")))
  }
  we <- m("particle_weight_ess")
  if (!is.na(we)) {
    say(sprintf("  %-20s %8.1f     (one proposal cloud, NOT chain ESS)\n",
                "weight_ess", we))
  }
  if (!is.na(m("particles_sum"))) {
    say(sprintf("  %-20s %8.0f     (max %.0f per subject)\n", "particles_sum",
                m("particles_sum"), m("particles_max")))
  }

  if (!all(is.na(profile$wire_bytes))) {
    say("\nwire per iteration\n")
    for (nm in c("shared_bytes", "private_bytes", "wire_bytes")) {
      say(sprintf("  %-20s %s\n", nm, .emc_profile_fmt_bytes(profile[[nm]])))
    }
  }
  if (!all(is.na(full$private_memory))) {
    meth <- full$private_memory_method
    meth <- meth[!is.na(meth)]
    say(sprintf("\nprivate memory       %s  (%s, peak %s, %d of %d iterations sampled)\n",
                .emc_profile_fmt_bytes(full$private_memory),
                if (length(meth)) meth[[1L]] else "unknown",
                .emc_profile_fmt_bytes(max(full$private_memory, na.rm = TRUE)),
                sum(!is.na(full$private_memory)), nrow(full)))
  }

  kr <- m("kernel_rows")
  kc <- m("kernel_cells")
  if (!is.na(kr) && !is.na(kc) && kc > 0) {
    say("\nkernel reuse per iteration\n")
    say(sprintf("  %-20s %10.0f   per-trial evaluations\n", "kernel_rows", kr))
    say(sprintf("  %-20s %10.0f   at cell resolution\n", "kernel_cells", kc))
    say(sprintf("  %-20s %10.1fx  (1.0 means nothing to hoist)\n", "reuse", kr / kc))
    ks <- m("kernel_seconds")
    if (!is.na(ks)) {
      say(sprintf("  %-20s %10.2f ms  %s\n", "kernel_seconds", 1000 * ks,
                  if (!is.na(tot) && tot > 0) sprintf("%5.1f%% of the iteration",
                                                      100 * ks / tot) else ""))
    }
  }

  if (!is.null(blocks) && nrow(blocks)) {
    .emc_block_report(blocks, elapsed = elapsed, con = con)
  }

  health <- Filter(function(f) f$group == "health", .emc_profile_schema)
  totals <- vapply(health, function(f) sum(full[[f$name]], na.rm = TRUE),
                   numeric(1))
  if (any(totals > 0)) {
    say("\nhealth (totals over all iterations)\n")
    for (i in which(totals > 0)) {
      say(sprintf("  %-20s %8.0f   %s\n", health[[i]]$name, totals[[i]],
                  health[[i]]$desc))
    }
    for (key in names(.emc_reject_messages())) {
      say(sprintf("    %-18s %s\n", key, .emc_reject_messages()[[key]]))
    }
  }
  invisible(profile)
}
