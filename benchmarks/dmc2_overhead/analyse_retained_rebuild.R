args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: analyse_retained_rebuild.R OUT_DIR SEED")
out_dir <- args[[1L]]
seed <- as.integer(args[[2L]])
baseline <- readRDS(file.path(out_dir,
  sprintf("retained_rebuild_%d_baseline.rds", seed)))
candidate <- readRDS(file.path(out_dir,
  sprintf("retained_rebuild_%d_candidate.rds", seed)))

same_samples <- function(a, b) {
  if (length(a) != length(b)) return(FALSE)
  all(vapply(seq_along(a), function(j) {
    fields <- union(names(a[[j]]$samples), names(b[[j]]$samples))
    all(vapply(fields, function(nm)
      identical(a[[j]]$samples[[nm]], b[[j]]$samples[[nm]]), logical(1)))
  }, logical(1)))
}
strip_timing <- function(x) {
  if (!is.list(x)) return(x)
  nms <- names(x)
  if (!is.null(nms)) {
    timing <- grepl("cpu|cost|wall|elapsed|timing|seconds", nms,
                    ignore.case = TRUE)
    x <- x[!timing]
  }
  x[] <- lapply(x, strip_timing)
  x
}
same_field <- function(a, b, field) {
  length(a) == length(b) && all(vapply(seq_along(a), function(j)
    identical(strip_timing(a[[j]][[field]]),
              strip_timing(b[[j]][[field]])), logical(1)))
}
fields <- setdiff(union(names(baseline$fit[[1L]]),
                        names(candidate$fit[[1L]])),
                  c("data", "prior", "model"))
state_fields <- setdiff(fields, c("alloc", "subject_hmc"))
same_continuation <- all(vapply(state_fields, function(nm)
  same_field(baseline$fit, candidate$fit, nm), logical(1)))
same_rng <- same_field(baseline$fit, candidate$fit, "rng")
same_caller_rng <- identical(baseline$caller_seed, candidate$caller_seed)

allocation_schedule <- function(chain) {
  record <- function(x) x[intersect(c("stage", "B", "choice", "rounds",
    "collapsed_rounds", "rate", "particle_rate", "particle_fixed",
    "choice_fixed"), names(x))]
  alloc <- chain$alloc
  hmc <- chain$subject_hmc
  hmc$alloc <- NULL
  list(particle = alloc$particle[setdiff(names(alloc$particle), "cost")],
       pool = alloc$pool,
       current = record(alloc$last),
       history = lapply(alloc$history, record),
       shear = chain$shear$plan[intersect(c("dirs", "lam", "w", "log_step",
         "rounds", "stage"), names(chain$shear$plan))],
       collapsed = chain$collapsed$plan[intersect(c("C", "eps", "rounds",
         "dim", "stage"), names(chain$collapsed$plan))],
       hmc = strip_timing(hmc))
}
same_schedule <- all(vapply(seq_along(baseline$fit), function(j)
  identical(allocation_schedule(baseline$fit[[j]]),
            allocation_schedule(candidate$fit[[j]])), logical(1)))
exact <- same_samples(baseline$fit, candidate$fit) && same_rng &&
  same_continuation && same_caller_rng && same_schedule

trace <- candidate$rebuild_worker_trace
trace_parts <- strsplit(trimws(trace), "[[:space:]]+")
trace_ok <- length(trace_parts) >= 3L && all(vapply(trace_parts, function(x)
  length(x) >= 4L && identical(x[[3L]], "subject_pool") &&
    length(strsplit(x[[4L]], ",", fixed = TRUE)[[1L]]) >= 2L,
  logical(1)))
worker_pids_stable <- if (trace_ok) {
  grouped <- split(seq_along(trace_parts),
                   vapply(trace_parts, `[[`, character(1), 1L))
  worker_sets <- lapply(grouped, function(ii)
    unique(vapply(trace_parts[ii], `[[`, character(1), 4L)))
  length(grouped) == 3L && all(lengths(grouped) >= 2L) &&
    all(lengths(worker_sets) == 1L)
} else FALSE
memory <- lapply(trace_parts, function(x) {
  if (length(x) < 5L) return(NULL)
  values <- suppressWarnings(as.numeric(strsplit(x[[5L]], ",", fixed = TRUE)[[1L]]))
  if (length(values) != 6L) return(NULL)
  c(chain_pid = as.numeric(x[[1L]]),
    total_pss_mib = sum(values[c(2L, 5L)], na.rm = TRUE) / 1024,
    total_private_mib = sum(values[c(3L, 6L)], na.rm = TRUE) / 1024)
})
memory <- Filter(Negate(is.null), memory)
memory_summary <- if (length(memory)) {
  memory <- as.data.frame(do.call(rbind, memory))
  groups <- split(memory, memory$chain_pid)
  paste(vapply(names(groups), function(pid) {
    x <- groups[[pid]]
    sprintf("chain %s PSS %.1f->%.1f MiB, private %.1f->%.1f MiB",
      pid, x$total_pss_mib[[1L]], x$total_pss_mib[[nrow(x)]],
      x$total_private_mib[[1L]], x$total_private_mib[[nrow(x)]])
  }, character(1)), collapse = "; ")
} else "not recorded"
mode_failures <- sum(vapply(candidate$fit, function(chain)
  sum(chain$burn_mode_stats$laplace_status == "mode failed", na.rm = TRUE),
  numeric(1)))

cat(sprintf(paste0("seed %d: baseline %.3f s, candidate %.3f s, speedup %.3fx; ",
                   "samples=%s rng=%s continuation_without_timing=%s ",
                   "schedule=%s caller_rng=%s; ",
                   "worker_calls=%d retained_worker_pids_stable=%s mode_failures=%d\n"),
            seed, baseline$elapsed[["elapsed"]], candidate$elapsed[["elapsed"]],
            baseline$elapsed[["elapsed"]] / candidate$elapsed[["elapsed"]],
            same_samples(baseline$fit, candidate$fit), same_rng,
            same_continuation, same_schedule, same_caller_rng, length(trace),
            worker_pids_stable, mode_failures))
cat("rebuild-entry retained memory (chain + subject workers): ",
    memory_summary, "\n", sep = "")
if (!exact) {
  different <- state_fields[!vapply(state_fields, function(nm)
    same_field(baseline$fit, candidate$fit, nm), logical(1))]
  cat("different continuation fields:", paste(different, collapse = ", "), "\n")
  if (!same_schedule) cat("allocation/HMC schedule differs\n")
  stop("baseline and retained-rebuild candidate did not preserve exact state")
}
if (!trace_ok || !worker_pids_stable)
  stop("candidate rebuilds did not stay on stable retained subject workers")
if (mode_failures > 0L) stop("candidate Laplace refresh reported mode failures")
if (!isTRUE(baseline$fixed_allocation) || !isTRUE(candidate$fixed_allocation))
  stop("exact-state comparison requires fixed-allocation arms")
