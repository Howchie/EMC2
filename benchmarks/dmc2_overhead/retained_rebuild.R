# Replay two fixed 100-transition sample blocks from the same saved state.
# Persistent chains are enabled in both arms; the candidate additionally runs
# Laplace rebuild/check jobs on the subject workers retained by each chain.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
input <- Sys.getenv("FIXTURE")
output <- Sys.getenv("OUT")
trace_file <- Sys.getenv("REBUILD_TRACE")
seed <- as.integer(Sys.getenv("SEED", "101"))
cores_per_chain <- as.integer(Sys.getenv("CORES_PER_CHAIN", "2"))
active_chains <- as.integer(Sys.getenv("CHAINS", "2"))
if (!nzchar(input) || !file.exists(input) || !nzchar(output))
  stop("set FIXTURE and OUT to existing/desired RDS paths")

options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE,
        emc2.persistent_chains = TRUE,
        emc2.persistent_chain_restart_blocks = 25L,
        emc2.rebuild_worker_trace_file = if (nzchar(trace_file)) trace_file else NULL)
if (nzchar(trace_file)) unlink(trace_file)

# Hold the transition schedule fixed so timing noise in the adaptive allocator
# cannot change draws between arms. Cost measurements are still collected.
assignInNamespace(".emc_alloc_solve_units", function(u) {
  hmc <- vapply(u$units, function(unit)
    identical(unit$kind, "subject") && !is.null(unit$h), logical(1))
  x <- vapply(u$units, function(unit)
    if (identical(unit$kind, "shear")) 1 else 0, numeric(1))
  x[hmc] <- 0.5
  list(x = x, q = ifelse(hmc, 0.5, 1),
       joint = if (!is.null(u$joint)) 1,
       rho = NA_real_, tau = NA_real_, cost = u$C_rest,
       objective = NA_real_, objective0 = NA_real_, tau0 = NA_real_)
}, asNamespace("EMC2"))

# Log from each retained chain process so the benchmark confirms that rebuild
# tasks reached its retained subject workers rather than falling back locally.
if (exists(".emc_chain_process_rebuild", envir = asNamespace("EMC2"),
           inherits = FALSE) && nzchar(trace_file)) {
  tracer <- quote({
    trace_path <- getOption("emc2.rebuild_worker_trace_file", "")
    subject_pool <- .emc_pool_state$subject_pool
    pool_name <- "subject_pool"
    if (is.null(subject_pool) || !isTRUE(subject_pool$alive)) {
      subject_pool <- .emc_pool_state$subject_ll_pool
      pool_name <- "subject_ll_pool"
    }
    worker_pids <- if (is.null(subject_pool) || !length(subject_pool$jobs)) {
      integer()
    } else {
      vapply(subject_pool$jobs, function(job) as.integer(job$pid), integer(1))
    }
    read_smaps <- function(pid) {
      lines <- tryCatch(readLines(sprintf("/proc/%d/smaps_rollup", pid),
                                 warn = FALSE), error = function(e) character())
      value <- function(key) {
        hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
        if (!length(hit)) return(NA_real_)
        suppressWarnings(as.numeric(sub(
          "^[^:]+:\\s*([0-9]+).*$", "\\1", hit[[1L]])))
      }
      c(rss = value("Rss"), pss = value("Pss"),
        private = value("Private_Clean") + value("Private_Dirty"))
    }
    self_mem <- read_smaps(Sys.getpid())
    worker_mem <- if (length(worker_pids)) {
      do.call(rbind, lapply(worker_pids, read_smaps))
    } else matrix(numeric(), nrow = 0L, ncol = 3L,
                  dimnames = list(NULL, c("rss", "pss", "private")))
    worker_total <- vapply(colnames(worker_mem), function(nm) {
      x <- worker_mem[, nm]
      if (!length(x) || all(!is.finite(x))) NA_real_ else sum(x, na.rm = TRUE)
    }, numeric(1))
    if (length(trace_path) && nzchar(trace_path)) {
      cat(Sys.getpid(), length(tasks), pool_name,
          paste(worker_pids, collapse = ","),
          paste(c(round(self_mem), round(worker_total)), collapse = ","),
          "\n",
          file = trace_path, append = TRUE, sep = " ")
    }
  })
  trace(EMC2:::.emc_chain_process_rebuild, tracer = tracer, print = FALSE)
}

saved <- readRDS(input)
emc <- if (!is.null(saved$fit)) saved$fit else saved
emc <- EMC2:::restore_custom_kernel_pointers(emc, quiet = TRUE)
sample_n <- sum(emc[[1L]]$samples$stage[
  seq_len(emc[[1L]]$samples$idx)] == "sample")
if (sample_n < 100L) stop("fixture needs at least 100 sample draws")
# Retain a 100-draw conditioning window and position the rebuild clock so the
# second sample block refreshes while the persistent chain pool is alive.
emc <- subset(emc, stage = "sample", filter = sample_n - 100L,
              keep_stages = TRUE)
for (j in seq_along(emc)) {
  emc[[j]]$sample_transitions <- 100L
  emc[[j]]$sample_proposals_next <- 100L
}
input_md5 <- unname(tools::md5sum(input))
set.seed(seed)
started <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(
  emc, stage = "sample", stop_criteria = list(iter = 200L),
  step_size = 100L, trim = FALSE, verbose = FALSE,
  cores_per_chain = cores_per_chain, cores_for_chains = active_chains))
elapsed <- proc.time() - started
trace_lines <- if (nzchar(trace_file) && file.exists(trace_file))
  readLines(trace_file, warn = FALSE) else character()
saveRDS(list(fit = out, caller_seed = .Random.seed, elapsed = elapsed,
             seed = seed, fixture_md5 = input_md5,
             fixed_allocation = TRUE,
             rebuild_worker_trace = trace_lines,
             blocks = EMC2:::.emc_block_profile(),
             iterations = EMC2:::.emc_iteration_profile(),
             session = sessionInfo()), output)
cat("elapsed", elapsed[["elapsed"]], "seed", seed,
    "fixture_md5", input_md5, "rebuild_worker_calls", length(trace_lines), "\n")
