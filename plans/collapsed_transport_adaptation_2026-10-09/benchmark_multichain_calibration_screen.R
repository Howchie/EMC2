# Compare the exact full pilot search with the warm-start screen on the same
# frozen three-chain LBA checkpoint. Pilot cost is included in fit wall/CPU.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[[1L]] %in% c("full", "screen"))
arm <- args[[1L]]
out_file <- args[[2L]]
continuation <- as.integer(Sys.getenv("COLL_SCREEN_DRAWS", "200"))
cores_per_chain <- as.integer(Sys.getenv("COLL_SCREEN_CORES", "4"))
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
snapshot <- Sys.getenv("COLL_SCREEN_SNAPSHOT",
  "/tmp/lba_sampling_audit/calibrated_snapshot.rds")
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)
base_gibbs <- lapply(emc, function(ch) ch$rng$gibbs)
base_collapsed <- lapply(emc, function(ch)
  ch$collapsed[c("accepted", "rounds", "iterations")])
for (j in seq_along(emc)) {
  emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions
  emc[[j]]$production_options$emc2.subject_hmc <- FALSE
  emc[[j]]$production_options$emc2.collapsed_move <- TRUE
  # Missing calibration triggers the full search; a calibrated plan is
  # warm-started by the neighborhood screen at this same-stage rebuild.
  if (identical(arm, "full"))
    emc[[j]]$collapsed$plan$calibration <- NULL
}
options(emc2.subject_hmc = FALSE, emc2.collapsed_move = TRUE)

cg_file <- file.path("/sys/fs/cgroup",
  sub("^0::", "", readLines("/proc/self/cgroup")[1L]), "cpu.stat")
cg_cpu <- function() {
  x <- tryCatch(readLines(cg_file), error = function(e) character())
  u <- grep("^usage_usec", x, value = TRUE)
  if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_
}
cpu0 <- cg_cpu()
started <- proc.time()[["elapsed"]]
emc <- EMC2:::add_proposals(emc, "sample",
  n_cores = length(emc) * cores_per_chain, n_blocks = 1L,
  block_iter = continuation)
rebuild_wall <- proc.time()[["elapsed"]] - started
calibration_cpu <- vapply(emc, function(ch)
  as.numeric(ch$collapsed$plan$calibration$elapsed_core_seconds %||%
    ch$collapsed$plan$calibration$elapsed %||% NA_real_), numeric(1))
# Pilots are discarded. Start both arms' retained transitions from the same
# Gibbs stream so calibration randomness does not confound sampling draws.
for (j in seq_along(emc)) {
  emc[[j]]$rng$gibbs <- base_gibbs[[j]]
  emc[[j]]$sample_proposals_next <- Inf
}
sample_started <- proc.time()[["elapsed"]]
fit <- EMC2::fit(emc, stage = "sample",
  stop_criteria = list(sample = list(iter = continuation)),
  step_size = continuation, cores_for_chains = length(emc),
  cores_per_chain = cores_per_chain, particle_factor = 10,
  max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
sample_wall <- proc.time()[["elapsed"]] - sample_started
wall <- proc.time()[["elapsed"]] - started
cpu <- cg_cpu() - cpu0

pn <- fit[[1L]]$par_names
p <- length(pn)
n <- min(continuation, vapply(fit, function(ch)
  sum(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), integer(1)))
stopifnot(n == continuation)
ids <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
make_array <- function(getter, variables) {
  out <- array(NA_real_, c(n, length(fit), length(variables)),
    dimnames = list(NULL, NULL, variables))
  for (j in seq_along(fit)) {
    ch <- fit[[j]]
    ix <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n)
    out[, j, ] <- getter(ch, ix)
  }
  posterior::as_draws_array(out)
}
summarise <- function(a) {
  d <- posterior::summarise_draws(a, posterior::ess_bulk,
    posterior::ess_tail, posterior::rhat)
  names(d)[match("posterior::ess_bulk", names(d))] <- "ess_bulk"
  names(d)[match("posterior::ess_tail", names(d))] <- "ess_tail"
  names(d)[match("posterior::rhat", names(d))] <- "rhat"
  d
}
group <- summarise(make_array(function(ch, ix)
  t(ch$samples$theta_mu[, ix, drop = FALSE]), pn))
logsd <- summarise(make_array(function(ch, ix)
  t(vapply(ix, function(k) log(sqrt(diag(ch$samples$theta_var[, , k]))),
    numeric(p))), paste0("logsd:", pn)))
subject_id <- "52091"
s <- match(subject_id, ids)
subject <- summarise(make_array(function(ch, ix)
  t(matrix(ch$samples$alpha[, s, ix], nrow = p, ncol = n)), pn))
all_person_names <- as.vector(outer(pn, ids, paste, sep = "@"))
all_person <- summarise(make_array(function(ch, ix) {
  z <- matrix(NA_real_, nrow = n, ncol = p * length(ids))
  for (ss in seq_along(ids))
    z[, ((ss - 1L) * p + 1L):(ss * p)] <-
      t(matrix(ch$samples$alpha[, ss, ix], nrow = p, ncol = n))
  z
}, all_person_names))
families <- list(group_mu = group, group_logsd = logsd,
  subject_52091 = subject, all_subject = all_person)
family_summary <- do.call(rbind, lapply(names(families), function(nm) {
  d <- families[[nm]]
  data.frame(family = nm, min_bulk_ess = min(d$ess_bulk),
    median_bulk_ess = median(d$ess_bulk), mean_bulk_ess = mean(d$ess_bulk),
    max_rhat = max(d$rhat, na.rm = TRUE),
    frac_rhat_over_1_05 = mean(d$rhat > 1.05, na.rm = TRUE))
}))
per_chain <- lapply(seq_along(fit), function(j) {
  ch <- fit[[j]]
  b <- base_collapsed[[j]]
  data.frame(chain = j,
    scale = ch$collapsed$plan$scale,
    t0_mean_scale = ch$collapsed$plan$scale * ch$collapsed$plan$t0_ratio,
    t0_hyper_scale = ch$collapsed$plan$scale * ch$collapsed$plan$t0_hyper_ratio,
    mean_scale = ch$collapsed$plan$scale * ch$collapsed$plan$mean_ratio,
    path_length = ch$collapsed$plan$length,
    rounds_per_iteration = ch$collapsed$plan$rounds,
    calibration_mode = ch$collapsed$plan$calibration$search_mode,
    pilots = ch$collapsed$plan$calibration$pilots,
    exact_acceptance = if (ch$collapsed$plan$block$n > 0)
      ch$collapsed$plan$block$acc / ch$collapsed$plan$block$n else NA_real_,
    accepted_delta = ch$collapsed$accepted - b$accepted,
    rounds_delta = ch$collapsed$rounds - b$rounds,
    calibration_seconds = calibration_cpu[j])
})
result <- list(arm = arm, draws_per_chain = n,
  total_wall_seconds = wall, rebuild_wall_seconds = rebuild_wall,
  sample_wall_seconds = sample_wall, cgroup_cpu_seconds = cpu,
  family_summary = family_summary,
  per_variable = lapply(families, function(d)
    d[, c("variable", "ess_bulk", "ess_tail", "rhat")]),
  chain = do.call(rbind, per_chain))
saveRDS(result, out_file)
cat(sprintf("%s: %d draws/chain, wall %.1fs (rebuild %.1fs, sample %.1fs), CPU %.1fs\n",
  arm, n, wall, rebuild_wall, sample_wall, cpu))
print(family_summary, row.names = FALSE)
print(result$chain, row.names = FALSE)
