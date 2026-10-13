# Short continuation (3 chains x 3 cores) from a saved sample-stage checkpoint under one ensemble
# setting: fixed pool and rounds (arm "M_R", e.g. 32_10) or the between-block
# plan (arm "rule"). Rebuilds every 100 transitions in every arm. Reports
# worst-coordinate ESS of the group means and log sds per CPU-second and per
# wall-second, lag-k autocorrelation of the worst coordinate against the
# checkpoint's variance, and the plans chosen.
#   EMC_LIB=<lib> Rscript cont_mr.R <fit.rds> <arm> <iter> <out.rds> [seed]
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(EMC2, lib.loc = Sys.getenv("EMC_LIB")); library(posterior) })
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
fx <- args[1]; arm <- args[2]; iter <- as.integer(args[3]); out <- args[4]
seed <- if (length(args) >= 5) as.integer(args[5]) else 1L
emc <- readRDS(fx)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
p <- emc[[1]]$n_pars
n0 <- emc[[1]]$sample_transitions
if (is.null(n0)) n0 <- sum(emc[[1]]$samples$stage == "sample")
for (j in seq_along(emc)) {
  emc[[j]]$production_options$emc2.sample_rebuilds <- as.integer(seq(n0 + 100L, n0 + iter + 100L, by = 100L))
  emc[[j]]$sample_proposals_next <- n0 + 100L
  emc[[j]]$ensemble_stats$block <- NULL
  emc[[j]]$ensemble_plan_history <- NULL
  if (arm == "rule") emc[[j]]$ensemble_plan <- NULL else {
    mr <- as.integer(strsplit(arm, "_")[[1]])
    emc[[j]]$ensemble_plan <- list(pool = mr[1], rounds = mr[2])
  }
}
if (arm != "rule") {
  f <- function(sampler) sampler
  environment(f) <- ns
  assignInNamespace(".emc_ensemble_plan_update", f, ns)
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
c0 <- cg_cpu(); t0 <- proc.time()[["elapsed"]]
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample", iter = iter, cores_per_chain = 3,
                                             cores_for_chains = 3, verbose = FALSE)))
wall <- proc.time()[["elapsed"]] - t0; cpu <- cg_cpu() - c0
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
plans <- lapply(emc, function(ch) ch$ensemble_plan_history)
plan_tab <- do.call(rbind, lapply(seq_along(plans), function(j) do.call(rbind, lapply(plans[[j]], function(r)
  data.frame(chain = j, n = r$n, M = r$M, R = r$R, G = r$G, J_tot = r$J_tot, C = r$C, c_member = r$c_member, T_min = if (is.null(r$T_min)) NA else r$T_min, G_min = if (is.null(r$G_min)) NA else r$G_min, worst = if (is.null(r$worst)) NA else r$worst, eff = r$eff, moved = r$moved, pool = r$pool, rounds = r$rounds)))))
res <- list(fixture = basename(fx), arm = arm, seed = seed, iter = iter, n_new = n_new, wall = wall, cpu = cpu,
            ess = ess, worst = names(ess)[worst], ess_min = min(ess), ess_min_per_cpu = min(ess) / cpu,
            ess_min_per_wall = min(ess) / wall, lag = lag, plans = plan_tab,
            ess_alpha_min = min(ess_alpha), ess_alpha_med = median(ess_alpha),
            diag = ns$ensemble_diagnostics(emc), lib = find.package("EMC2"))
saveRDS(res, out)
cat(sprintf("%s %s seed %d: %d draws/chain, wall %.0f s, cpu %.0f s | worst %s ESS %.0f (per cpu-s %.3f, per wall-s %.2f) | lag1-5 of worst: %s\n",
  res$fixture, arm, seed, n_new, wall, cpu, res$worst, min(ess), min(ess) / cpu, min(ess) / wall, paste(sprintf("%.2f", lag[1:5]), collapse = " ")))
cat(sprintf("  min ESS mu %.0f, lsd %.0f; median %.0f | alpha min %.0f median %.0f\n", min(ess[1:p]), min(ess[-(1:p)]), median(ess),
            min(ess_alpha), median(ess_alpha)))
if (!is.null(plan_tab)) print(plan_tab, digits = 3, row.names = FALSE)
