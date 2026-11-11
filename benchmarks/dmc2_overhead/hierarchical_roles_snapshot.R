# Replay a saved hierarchical chain through two sample blocks.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
role_dir <- Sys.getenv('ROLE_DIR')
stopifnot(nzchar(role_dir), dir.exists(role_dir))
ns <- asNamespace('EMC2')
record_role <- function(role, pid) {
  if (!length(pid) || is.na(pid)) return(invisible(NULL))
  lines <- tryCatch(readLines(sprintf('/proc/%s/smaps_rollup', pid), warn = FALSE),
                    error = function(e) character())
  value <- function(key) {
    line <- grep(paste0('^', key, ':'), lines, value = TRUE)
    if (!length(line)) return(NA_real_)
    as.numeric(strsplit(trimws(line), '[[:space:]]+')[[1]][2])
  }
  cat(role, as.integer(pid), value('Rss'), value('Pss'),
      value('Private_Clean') + value('Private_Dirty'), "\n",
      file = file.path(role_dir, paste0(Sys.getpid(), '.log')), append = TRUE)
}
record_role('master', Sys.getpid())
orig_stage <- get('run_stage', ns)
wrap_stage <- function(...) { record_role('chain', Sys.getpid()); on.exit(record_role('chain', Sys.getpid())); orig_stage(...) }
assignInNamespace('run_stage', wrap_stage, ns)
orig_session <- get('.emc_wpool_session_start', ns)
wrap_session <- function(...) {
  result <- orig_session(...)
  if (!is.null(result)) for (t in result$templates) if (!is.null(t))
    record_role('template', ns$.emc_wpool_job_pid(t$job))
  result
}
assignInNamespace('.emc_wpool_session_start', wrap_session, ns)
orig_spawn <- get('.emc_wpool_spawn', ns)
wrap_spawn <- function(...) {
  result <- orig_spawn(...)
  if (!is.null(result)) for (j in result$jobs) record_role('worker', ns$.emc_wpool_job_pid(j))
  result
}
assignInNamespace('.emc_wpool_spawn', wrap_spawn, ns)
input <- Sys.getenv('FIXTURE')
stopifnot(nzchar(input), file.exists(input))
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE, emc2.sampler_profile_expensive = TRUE)
emc <- EMC2:::restore_custom_kernel_pointers(readRDS(input), quiet = TRUE)
set.seed(as.integer(Sys.getenv('SEED', '101')))
start <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(emc, stage = 'sample',
  stop_criteria = list(iter = 20), step_size = 10, verbose = FALSE,
  cores_per_chain = as.integer(Sys.getenv('CORES_PER_CHAIN', '2')),
  cores_for_chains = as.integer(Sys.getenv('CHAINS', '2'))))
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, timing = timing,
             blocks = EMC2:::.emc_block_profile(), iterations = EMC2:::.emc_iteration_profile(),
             fixture_md5 = unname(tools::md5sum(input)), session = sessionInfo()), Sys.getenv('OUT'))
cat('timing', timing, '\n')
