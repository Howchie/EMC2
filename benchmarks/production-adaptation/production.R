# One production arm from a saved warm-up state (warmup.R).
#   CASE=wide15 SEED=1 ARM=doubling EMC_LIB=<lib> Rscript benchmarks/production-adaptation/production.R
# Arms (internal option emc2.sample_rebuilds):
#   doubling  rebuild after 100, 200, 400, 800 production transitions (default)
#   frozen    proposals built at production start, never rebuilt
#   early     rebuild after 100 and 300 transitions, then frozen
#   extended  400 more adapt iterations (discarded), then frozen
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
suppressPackageStartupMessages(library(posterior))
out <- Sys.getenv("PADAPT_OUT", "benchmarks/production-adaptation/results")
case <- Sys.getenv("CASE")
seed <- as.integer(Sys.getenv("SEED", "1"))
arm <- Sys.getenv("ARM", "doubling")
# fit()'s sample-stage rule: at least SAMPLE_ITER iterations AND max R-hat < 1.1
# on alpha and on mu, giving up after max_tries blocks.
n_iter <- as.integer(Sys.getenv("SAMPLE_ITER", "1000"))
max_tries <- as.integer(Sys.getenv("MAX_TRIES", "20"))
cores <- as.integer(Sys.getenv("CORES_PER_CHAIN", "2"))
emc <- readRDS(file.path(out, sprintf("warmup_%s_%d.rds", case, seed)))
warm <- attr(emc, "warmup_info")
# Same production stream for every arm of a (case, seed).
set.seed(5000L + seed)
# INDEP=1 runs the chains fully independently (no pooled windows, no shared
# modes); TAG names the variant in the output file.
if (Sys.getenv("INDEP", "0") == "1") options(emc2.independent_chains = TRUE)
# GROUP_MOVE=shear follows the ensemble group update with the shear move.
options(emc2.group_move = Sys.getenv("GROUP_MOVE", "ensemble"))
# SAVE_EMC=1 also saves the fitted object (all draws) next to the summary.
tag <- Sys.getenv("TAG", "")
options(emc2.sample_rebuilds = switch(arm,
  doubling = NULL, frozen = numeric(0), early = c(100, 300), extended = numeric(0),
  stop("unknown ARM")))

cpu <- function() { p <- proc.time(); sum(p[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) }
run <- function(emc, stage, sc) suppressWarnings(run_emc(emc, stage = stage,
  stop_criteria = sc, cores_for_chains = 3, cores_per_chain = cores,
  step_size = 100, max_tries = max_tries, verbose = FALSE, trim = FALSE))
ext_cpu <- 0; ext_elapsed <- 0
if (arm == "extended") {
  c0 <- cpu(); t0 <- proc.time()[["elapsed"]]
  emc <- run(emc, "adapt", list(iter = 400, min_unique = 150))
  ext_cpu <- cpu() - c0; ext_elapsed <- proc.time()[["elapsed"]] - t0
}
c0 <- cpu(); t0 <- proc.time()[["elapsed"]]
emc <- run(emc, "sample", list(iter = n_iter, max_gd = 1.1, omit_mpsrf = TRUE,
                               selection = c("alpha", "mu")))
prod_cpu <- cpu() - c0; prod_elapsed <- proc.time()[["elapsed"]] - t0

# How the sample stage ended. The last convergence-log entry holds the R-hat
# rules at the final check; the retained draws can be fewer than the
# transitions run, because a failed check may drop the first third.
log <- Filter(function(e) identical(e$stage, "sample"), emc[[1]]$convergence_log)
last <- log[[length(log)]]
gd_final <- last$gd$stats
stop_info <- data.frame(
  converged = !is.null(gd_final) && all(gd_final$value < gd_final$threshold),
  transitions = emc[[1]]$sample_transitions,
  retained = chain_n(emc)[1, "sample"],
  checks = length(log),
  rhat_alpha_max = if (is.null(gd_final)) NA else gd_final$value[gd_final$selection == "alpha"],
  rhat_mu_max = if (is.null(gd_final)) NA else gd_final$value[gd_final$selection == "mu"])

# Production draws, iterations x chains x variables.
emc <- EMC2:::restore_duplicates(emc)
sidx <- lapply(emc, function(ch) which(ch$samples$stage == "sample" &
                                         seq_along(ch$samples$stage) <= ch$samples$idx))
n <- min(lengths(sidx))
arr <- function(get) {
  x <- lapply(seq_along(emc), function(j) get(emc[[j]], tail(sidx[[j]], n)))
  out <- array(NA_real_, c(n, length(emc), ncol(x[[1]])),
               dimnames = list(NULL, NULL, colnames(x[[1]])))
  for (j in seq_along(x)) out[, j, ] <- x[[j]]
  out
}
mu <- arr(function(ch, i) t(ch$samples$theta_mu[, i, drop = FALSE]))
lsd <- arr(function(ch, i) {
  v <- ch$samples$theta_var
  m <- t(vapply(i, function(k) 0.5 * log(diag(as.matrix(v[, , k]))), numeric(dim(v)[1])))
  colnames(m) <- paste0("logsd_", rownames(ch$samples$theta_mu)); m
})
alpha <- arr(function(ch, i) {
  a <- ch$samples$alpha[, , i, drop = FALSE]
  m <- t(matrix(a, prod(dim(a)[1:2])))
  colnames(m) <- paste(rep(dimnames(a)[[1]], dim(a)[2]), rep(dimnames(a)[[2]], each = dim(a)[1]), sep = "|")
  m
})
diag_of <- function(x, what) {
  d <- as_draws_array(x)
  s <- summarise_draws(d, mean, sd, ess_bulk, ess_tail, rhat, mcse_mean)
  s$what <- what
  as.data.frame(s)
}
per_var <- rbind(diag_of(mu, "mu"), diag_of(lsd, "logsd"), diag_of(alpha, "alpha"))

# Movement per production segment (between rebuild times of the doubling arm):
# share of subject-iterations whose state changed, and the mean squared jump
# per iteration relative to the posterior variance.
seg_breaks <- c(0, 100, 200, 400, 800, 1600, n)
seg_breaks <- unique(pmin(seg_breaks, n))
pvar <- function(x) apply(x, 3, function(v) stats::var(as.vector(v)))
pv_alpha <- pvar(alpha)
pv_mu <- pvar(mu)
segs <- do.call(rbind, lapply(seq_len(length(seg_breaks) - 1L), function(k) {
  rows <- (seg_breaks[k] + 1L):seg_breaks[k + 1L]
  da <- apply(alpha[rows, , , drop = FALSE], c(2, 3), diff)
  dm <- apply(mu[rows, , , drop = FALSE], c(2, 3), diff)
  # Subject moved: any of its parameters changed.
  subj <- sub("^.*\\|", "", dimnames(alpha)[[3]])
  moved <- sapply(unique(subj), function(s) mean(apply(da[, , subj == s, drop = FALSE] != 0, c(1, 2), any)))
  data.frame(segment = sprintf("%d-%d", seg_breaks[k] + 1L, seg_breaks[k + 1L]),
             move_rate = mean(moved), move_rate_min = min(moved),
             esjd_alpha = mean(apply(da^2, 3, mean) / pv_alpha),
             esjd_mu = mean(apply(dm^2, 3, mean) / pv_mu))
}))

q <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
summ <- do.call(rbind, lapply(split(per_var, per_var$what), function(d) data.frame(
  what = d$what[1], ess_bulk_min = min(d$ess_bulk), ess_bulk_q05 = q(d$ess_bulk, 0.05),
  ess_bulk_med = median(d$ess_bulk), ess_tail_min = min(d$ess_tail),
  rhat_max = max(d$rhat))))
res <- list(case = case, seed = seed, arm = arm, tag = tag, n_iter = n, warmup = warm,
  ext_cpu = ext_cpu, ext_elapsed = ext_elapsed, prod_cpu = prod_cpu,
  prod_elapsed = prod_elapsed, stop = stop_info, summary = summ, segments = segs,
  per_var = per_var[, c("variable", "what", "mean", "sd", "ess_bulk", "ess_tail", "rhat", "mcse_mean")],
  rebuild_next = emc[[1]]$sample_proposals_next,
  transitions = emc[[1]]$sample_transitions)
print(stop_info); print(summ); print(segs)
cat(sprintf("%s %s seed %d: production %.0f s CPU (%.0f s elapsed)\n", case, arm, seed, prod_cpu, prod_elapsed))
stem <- sprintf("%s_%s%s_%d", case, arm, if (nzchar(tag)) paste0("_", tag) else "", seed)
saveRDS(res, file.path(out, paste0("prod_", stem, ".rds")))
if (Sys.getenv("SAVE_EMC", "0") == "1") saveRDS(emc, file.path(out, paste0("emc_", stem, ".rds")))
