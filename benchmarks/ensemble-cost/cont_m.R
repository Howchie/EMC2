# Continuation (3 chains x 3 cores) from a saved sample-stage checkpoint with
# a fixed ensemble pool size M and rounds R, exactly `iter` iterations
# (max_tries = 1), rebuilding every 100 transitions. Reports worst-coordinate
# ESS of the group means and log sds per wall-second and per process-tree
# CPU-second, the share of subjects the move refreshes, and lag-k
# autocorrelation of the worst coordinate against the checkpoint variance.
#   EMC_LIB=<lib> Rscript cont_m.R <fit.rds> <M> <R> <iter> <out.rds> [seed]
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(EMC2, lib.loc = Sys.getenv("EMC_LIB")); library(posterior) })
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
fx <- args[1]; M <- as.integer(args[2]); R <- as.integer(args[3]); iter <- as.integer(args[4]); out <- args[5]
seed <- if (length(args) >= 6) as.integer(args[6]) else 1L
arm <- paste0(M, "_", R)
emc <- readRDS(fx)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
p <- emc[[1]]$n_pars
n0 <- emc[[1]]$sample_transitions
if (is.null(n0)) n0 <- sum(emc[[1]]$samples$stage == "sample")
for (j in seq_along(emc)) {
  emc[[j]]$production_options$emc2.sample_rebuilds <- as.integer(seq(n0 + 100L, n0 + iter + 100L, by = 100L))
  emc[[j]]$production_options$emc2.ensemble_pool <- M
  emc[[j]]$production_options$emc2.ensemble_rounds <- R
  emc[[j]]$sample_proposals_next <- n0 + 100L
  emc[[j]]$ensemble_stats <- NULL
  # Fresh chain streams from the global seed, so seeds replicate.
  emc[[j]]$rng <- NULL
}
# Checkpoint variance of the group coordinates over all chains' sample draws.
coord_draws <- function(ch, w) rbind(ch$samples$theta_mu[, w, drop = FALSE],
  vapply(w, function(t) 0.5 * log(diag(as.matrix(ch$samples$theta_var[, , t]))), numeric(p)))
w0 <- lapply(emc, function(ch) which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx))
v0 <- apply(do.call(cbind, Map(coord_draws, emc, w0)), 1, var)
cg_file <- file.path("/sys/fs/cgroup", sub("^0::", "", readLines("/proc/self/cgroup")[1]), "cpu.stat")
cg_cpu <- function() { x <- tryCatch(readLines(cg_file), error = function(e) character(0))
  u <- grep("^usage_usec", x, value = TRUE); if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_ }
set.seed(seed)
c0 <- cg_cpu(); t0 <- proc.time()[["elapsed"]]; pt0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample", iter = iter, max_tries = 1,
                                             cores_per_chain = 3, cores_for_chains = 3, verbose = FALSE)))
wall <- proc.time()[["elapsed"]] - t0; cpu <- cg_cpu() - c0
# This process and its reaped descendants only (the cgroup may hold other work).
pt <- proc.time() - pt0
cpu_tree <- sum(pt[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
emc <- ns$restore_duplicates(emc)
n_new <- vapply(emc, function(ch) sum(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx), 1L) - vapply(w0, length, 1L)
n_new <- min(n_new)
draws <- lapply(emc, function(ch) {
  w <- which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx)
  w <- tail(w, n_new); coord_draws(ch, w)
})
arr <- aperm(simplify2array(draws), c(2, 3, 1))   # iterations x chains x coords
ess <- apply(arr, 3, ess_bulk)
alpha_draws <- lapply(emc, function(ch) {
  w <- which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx)
  w <- tail(w, n_new); matrix(ch$samples$alpha[, , w], ncol = n_new)
})
arr_a <- aperm(simplify2array(alpha_draws), c(2, 3, 1))
ess_alpha <- apply(arr_a, 3, ess_bulk)
names(ess) <- c(paste0("mu_", emc[[1]]$par_names), paste0("lsd_", emc[[1]]$par_names))
worst <- which.min(ess)
acf_k <- function(k) mean(vapply(seq_len(dim(arr)[2]), function(c) {
  x <- arr[, c, worst]; mean((x[-(1:k)] - mean(x)) * (x[seq_len(length(x) - k)] - mean(x))) / v0[worst] }, 0))
lag <- vapply(1:10, acf_k, 0)
dg <- ns$ensemble_diagnostics(emc)
res <- list(fixture = basename(fx), arm = arm, seed = seed, iter = iter, n_new = n_new, wall = wall, cpu = cpu, cpu_tree = cpu_tree,
            ess = ess, worst = names(ess)[worst], ess_min = min(ess), ess_min_per_cpu = min(ess) / cpu,
            ess_min_per_wall = min(ess) / wall, ess_min_per_cpu_tree = min(ess) / cpu_tree, lag = lag,
            ess_alpha_min = min(ess_alpha), ess_alpha_med = median(ess_alpha),
            M = M, R = R, moved = mean(dg$moved), wess = mean(dg$weight_ess),
            diag = dg, lib = find.package("EMC2"))
saveRDS(res, out)
cat(sprintf("%s %s seed %d: %d draws/chain, wall %.0f s, cpu %.0f s | worst %s ESS %.0f (per cpu-s %.3f, per wall-s %.2f) | lag1-5 of worst: %s\n",
  res$fixture, arm, seed, n_new, wall, cpu, res$worst, min(ess), min(ess) / cpu, min(ess) / wall, paste(sprintf("%.2f", lag[1:5]), collapse = " ")))
cat(sprintf("  process-tree cpu %.0f s (worst ESS per tree cpu-s %.3f)\n", cpu_tree, min(ess) / cpu_tree))
cat(sprintf("  min ESS mu %.0f, lsd %.0f; median %.0f | alpha min %.0f median %.0f\n", min(ess[1:p]), min(ess[-(1:p)]), median(ess),
            min(ess_alpha), median(ess_alpha)))
cat(sprintf("  M %d R %d: subjects refreshed per move %.3f, weight ESS/M %.3f\n", M, R, mean(dg$moved), mean(dg$weight_ess)))
