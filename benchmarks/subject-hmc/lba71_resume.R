# Production replays from the saved lba71sim checkpoint (subject 11, chain 1
# in a tail region; see README).
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/lba71_resume.R <mode> <arm> <out.rds>
# mode resume: an explicit, labelled proposal rebuild at the checkpoint's
#   transition count (both arms), then N_ITER production transitions.
# mode late: as resume for WARM transitions; then subject 11 of chain 2 is
#   moved to chain 1's tail state (likelihood recomputed exactly), the rebuild
#   schedule is set to rebuild REBUILD_AFTER transitions later, and N_ITER more
#   transitions run. This is an imposed test perturbation, not a posterior draw.
# arm: base (HMC off) | hmc.
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(EMC2, lib.loc = Sys.getenv("EMC_LIB")); library(posterior) })
cat("EMC2 from", find.package("EMC2"), "\n")
mode <- args[1]; arm <- args[2]; out <- args[3]
n_iter <- as.integer(Sys.getenv("N_ITER", "1000"))
warm <- as.integer(Sys.getenv("WARM", "100"))
rebuild_after <- as.integer(Sys.getenv("REBUILD_AFTER", "100"))
ckpt <- "benchmarks/production-adaptation/results_shear/shear_v4/emc_lba71sim_doubling_shearK4_2.rds"
stopifnot(unname(tools::md5sum(ckpt)) == "38c06da04e0c21fdda3539305f6c2fd9")
options(emc2.subject_hmc = arm == "hmc", emc2.group_move = "ensemble")
e <- EMC2:::restore_duplicates(readRDS(ckpt))
s <- 11L
for (j in seq_along(e)) e[[j]]$sample_proposals_next <- e[[j]]$sample_transitions
go <- function(e, n) suppressWarnings(suppressMessages(run_emc(e, stage = "sample",
  stop_criteria = list(iter = n), cores_for_chains = 3, cores_per_chain = 5, step_size = 100,
  max_tries = 1000, verbose = FALSE)))
start_idx <- e[[1]]$samples$idx
t0 <- proc.time()[["elapsed"]]
marks <- list()
if (mode == "late") {
  e <- go(e, warm)
  e <- EMC2:::restore_duplicates(e)
  st1 <- e[[1]]$samples
  tail_state <- readRDS(ckpt)[[1]]$samples$alpha[, s, start_idx]
  j <- 2L; i <- e[[j]]$samples$idx
  e[[j]]$samples$alpha[, s, i] <- tail_state
  ll <- EMC2:::calc_ll_manager(matrix(tail_state, 1, dimnames = list(NULL, e[[j]]$par_names)),
                               dadm = e[[j]]$data[[s]], model = e[[j]]$model)
  e[[j]]$samples$subj_ll[s, i] <- as.numeric(ll)
  n_now <- e[[1]]$sample_transitions
  sched <- c(n_now + rebuild_after, n_now + 2 * rebuild_after, n_now + 4 * rebuild_after,
             n_now + 8 * rebuild_after)
  for (jj in seq_along(e)) {
    e[[jj]]$production_options$emc2.sample_rebuilds <- sched
    e[[jj]]$sample_proposals_next <- n_now + rebuild_after
  }
  marks$imposed_at <- n_now
  marks$rate_before <- if (!is.null(e[[j]]$subject_hmc$alloc)) e[[j]]$subject_hmc$alloc$rate[s] else NA
}
e <- go(e, n_iter)
elapsed <- proc.time()[["elapsed"]] - t0
e <- EMC2:::restore_duplicates(e)
pn <- e[[1]]$par_names
new_idx <- function(ch) which(ch$samples$stage == "sample")
n <- min(vapply(e, function(ch) length(new_idx(ch)), 1L))
keep <- function(ch) utils::tail(new_idx(ch), n)
traj <- sapply(e, function(ch) ch$samples$alpha["A", s, keep(ch)])
blocks <- split(seq_len(n), ceiling(seq_len(n) / 100))
cat(sprintf("%s %s: %d transitions, elapsed %.0f s\n", mode, arm, n, elapsed))
for (b in blocks) cat(sprintf("  its %4d-%4d  subj 11 A mean (sd) by chain: %s\n", min(b), max(b),
  paste(sprintf("%.3f (%.3f)", colMeans(traj[b, , drop = FALSE]), apply(traj[b, , drop = FALSE], 2, sd)), collapse = "  ")))
arr <- function(f) simplify2array(lapply(e, function(ch) t(vapply(keep(ch), function(t) f(ch, t), numeric(length(f(ch, keep(ch)[1])))))))
al <- arr(function(ch, t) as.vector(ch$samples$alpha[, , t]))
mu <- arr(function(ch, t) ch$samples$theta_mu[, t])
dg <- function(a) { a <- aperm(a, c(1, 3, 2)); list(rhat = apply(a, 3, rhat), ess = apply(a, 3, ess_bulk)) }
da <- dg(al); dm <- dg(mu)
sub11 <- (s - 1) * length(pn) + seq_along(pn)
cat(sprintf("  subject 11: max R-hat %.3f min ESS %.0f | all alpha: max R-hat %.3f, subjects with any R-hat > 1.1: %d | mu: max R-hat %.3f min ESS %.0f\n",
  max(da$rhat[sub11]), min(da$ess[sub11]), max(da$rhat),
  sum(tapply(da$rhat > 1.1, rep(seq_len(e[[1]]$n_subjects), each = length(pn)), any)), max(dm$rhat), min(dm$ess)))
for (j in seq_along(e)) {
  h <- e[[j]]$subject_hmc
  if (is.null(h)) next
  cat(sprintf("  chain %d HMC: attempts %d accepted %d aborts %d rows %.3g calib %.3g\n", j, h$totals$attempts,
    h$totals$accepted, h$totals$aborts, h$totals$rows, h$totals$calib_rows))
  for (r in h$history) cat(sprintf("    v%d at %s: subj 11 rate %.3f eps %.3f rhoP %.3f rhoH %s | rates>0 %d, mean %.3f | tau_G %.1f tau0 max %.1f -> %.1f\n",
    r$version, format(r$transitions), r$rate[s], r$eps[s], r$rhoP_max[s], format(round(r$rhoH_max[s], 3)),
    sum(r$rate > 0), mean(r$rate), r$tau_G, max(r$tau0), r$tau))
}
saveRDS(list(mode = mode, arm = arm, traj = traj, rhat = da$rhat, ess = da$ess, mu_rhat = dm$rhat, mu_ess = dm$ess,
             elapsed = elapsed, marks = marks, hmc = lapply(e, function(ch) ch$subject_hmc[c("totals", "history")])), out)
