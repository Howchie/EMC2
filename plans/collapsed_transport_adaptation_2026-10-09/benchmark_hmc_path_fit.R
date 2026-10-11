# Paired whole-sampler continuation from the same LBA checkpoint. This
# includes the forced sample proposal rebuild and its HMC allocation pilots.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
label <- args[[1L]]
trajectory_length <- as.numeric(args[[2L]])
out_file <- args[[3L]]
continuation <- as.integer(Sys.getenv("HMC_FIT_DRAWS", "100"))
cores_per_chain <- as.integer(Sys.getenv("HMC_FIT_CORES", "4"))
block_text <- Sys.getenv("HMC_FIT_BLOCKS", "")
hmc_blocks <- if (nzchar(block_text))
  list(strsplit(block_text, ",", fixed = TRUE)[[1L]]) else NULL
seed_text <- Sys.getenv("HMC_FIT_SEED", "")
benchmark_seed <- if (nzchar(seed_text)) as.integer(seed_text) else NA_integer_
pkgload::load_all(".", quiet = TRUE, compile = FALSE)

snapshot <- Sys.getenv("HMC_FIT_SNAPSHOT",
  "/tmp/lba_sampling_audit/calibrated_snapshot.rds")
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)
if (is.finite(benchmark_seed)) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(benchmark_seed)
  chain_seed <- get(".Random.seed", envir = globalenv())
  for (j in seq_along(emc)) {
    chain_seed <- parallel::nextRNGStream(chain_seed)
    emc[[j]]$rng <- EMC2:::.emc_chain_rng(chain_seed, emc[[j]]$n_subjects)
  }
}
baseline_hmc <- lapply(emc, function(ch) ch$subject_hmc$totals)
options(emc2.subject_hmc = TRUE,
        emc2.subject_hmc_length = trajectory_length,
        emc2.subject_hmc_blocks = hmc_blocks,
        emc2.collapsed_move = FALSE)
for (j in seq_along(emc)) {
  emc[[j]]$production_options$emc2.subject_hmc <- TRUE
  emc[[j]]$production_options$emc2.subject_hmc_length <- trajectory_length
  emc[[j]]$production_options$emc2.subject_hmc_blocks <- hmc_blocks
  emc[[j]]$production_options$emc2.collapsed_move <- FALSE
  # Make both arms rebuild their sample proposals at the same starting state.
  emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions
}

cg_file <- file.path("/sys/fs/cgroup",
  sub("^0::", "", readLines("/proc/self/cgroup")[1L]), "cpu.stat")
cg_cpu <- function() {
  x <- tryCatch(readLines(cg_file), error = function(e) character())
  u <- grep("^usage_usec", x, value = TRUE)
  if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_
}
cpu0 <- cg_cpu()
wall0 <- proc.time()[["elapsed"]]
started <- proc.time()[["elapsed"]]
emc <- EMC2::fit(emc, stage = "sample",
  stop_criteria = list(sample = list(iter = continuation)),
  step_size = continuation, cores_for_chains = length(emc),
  cores_per_chain = cores_per_chain, particle_factor = 10,
  max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
wall <- proc.time()[["elapsed"]] - started
cpu <- cg_cpu() - cpu0

n <- min(continuation, vapply(emc, function(ch)
  sum(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), integer(1)))
subject <- "52091"
params <- c("A", "B", "B_lR_target", "t0_LoadHigh", "t0_LoadLow")
pn <- emc[[1L]]$par_names
group <- person <- log_sd <- array(NA_real_, c(n, length(emc), length(pn)),
  dimnames = list(NULL, NULL, pn))
for (j in seq_along(emc)) {
  ch <- emc[[j]]
  ix <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n)
  s <- match(subject, dimnames(ch$samples$alpha)[[2L]])
  group[, j, ] <- t(ch$samples$theta_mu[, ix, drop = FALSE])
  person[, j, ] <- t(matrix(ch$samples$alpha[, s, ix], nrow = length(pn), ncol = n))
  log_sd[, j, ] <- t(vapply(ix, function(k)
    0.5 * log(diag(as.matrix(ch$samples$theta_var[, , k]))), numeric(length(pn))))
}
ess_summary <- function(x) {
  out <- posterior::summarise_draws(posterior::as_draws_array(x),
    posterior::ess_bulk, posterior::ess_tail, posterior::rhat)
  names(out)[match("posterior::ess_bulk", names(out))] <- "ess_bulk"
  names(out)[match("posterior::ess_tail", names(out))] <- "ess_tail"
  names(out)[match("posterior::rhat", names(out))] <- "rhat"
  out
}
gd <- ess_summary(group)
pd <- ess_summary(person)
sd <- ess_summary(log_sd)
key_names <- c(paste0("mu:", params), paste0(subject, ":", params),
  paste0("log_sd:", params))
key <- array(NA_real_, c(n, length(emc), length(key_names)),
  dimnames = list(NULL, NULL, key_names))
for (j in seq_along(emc)) {
  ch <- emc[[j]]
  ix <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n)
  s <- match(subject, dimnames(ch$samples$alpha)[[2L]])
  key[, j, seq_along(params)] <- t(ch$samples$theta_mu[params, ix, drop = FALSE])
  key[, j, length(params) + seq_along(params)] <-
    t(matrix(ch$samples$alpha[params, s, ix], nrow = length(params), ncol = n))
  key[, j, 2L * length(params) + seq_along(params)] <-
    matrix(log_sd[, j, params, drop = FALSE], nrow = n, ncol = length(params))
}
kd <- ess_summary(key)
mov <- Map(function(ch, before) {
  after <- ch$subject_hmc$totals
  Map(function(a, b) a - b, after, before)
}, emc, baseline_hmc)
result <- list(label = label, trajectory_length = trajectory_length,
  hmc_blocks = hmc_blocks, rng_seed = benchmark_seed,
  continuation = n, elapsed_seconds = wall, cgroup_cpu_seconds = cpu,
  group_mu_ess = gd, subject_52091_ess = pd, group_log_sd_ess = sd,
  key_ess = kd, hmc_totals = mov,
  allocation_history = lapply(emc, function(ch) ch$alloc$history),
  hmc_block_names = lapply(emc, function(ch) lapply(
    ch$subject_hmc$subjects, `[[`, "block_names")),
  hmc_block_evidence = lapply(emc, function(ch) ch$subject_hmc$block_evidence),
  hmc_length_by_subject = lapply(emc, function(ch)
    vapply(ch$subject_hmc$subjects, `[[`, numeric(1), "trajectory_length")))
saveRDS(result, out_file)
cat(sprintf("%s length %.4f: %d draws/chain; wall %.1fs, cgroup CPU %.1fs\n",
  label, trajectory_length, n, wall, cpu))
for (label2 in c("group_mu_ess", "subject_52091_ess", "group_log_sd_ess", "key_ess")) {
  d <- result[[label2]]
  cat(sprintf("  %-20s min bulk %.1f, median %.1f; max Rhat %.3f\n", label2,
    min(d$ess_bulk), median(d$ess_bulk), max(d$rhat)))
}
for (j in seq_along(mov)) cat(sprintf(
  "  chain %d: HMC acceptance %.3f, abort %.3f, attempts %d, CPU %.1fs calibration %.1fs\n",
  j, mov[[j]]$accepted / max(1, mov[[j]]$attempts),
  mov[[j]]$aborts / max(1, mov[[j]]$attempts), mov[[j]]$attempts,
  mov[[j]]$cpu, mov[[j]]$calib_cpu))
