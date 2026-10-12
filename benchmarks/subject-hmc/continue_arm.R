# Stationary-efficiency arm: continue a converged checkpoint for a fixed number
# of production iterations and measure ESS per CPU second on quantities fixed
# by the checkpoint, so every arm is scored on the same posterior and the same
# directions.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/continue_arm.R <ckpt.rds> <label> <seed> <iters> <out.rds>
# Environment: CORES_PER_CHAIN (default 9); DISCARD (new iterations dropped
# before scoring, default 100: the first block, whose plan the checkpoint's
# build allocated); N_PARTICLES (particle count per subject replacing the
# stored one, with fit()'s particle_factor set to match so the adaptation's
# bounds follow it); ARM_R (an R file defining arm_options(), a list passed to
# options(), and/or arm_hook(emc) returning the emc to continue).
# CALIBRATE (e.g. "100,300"): production rebuilds at these transition offsets
# after the start, so the allocator reprices from its own measured costs before
# the scored window, then the doubling schedule's next rebuild (twice the last
# calibration point), beyond the run. DISCARD defaults to the last offset.
# The seed reseeds every chain's RNG streams; the state is the checkpoint's.
# A production rebuild is forced at the start so each arm allocates under its
# own kernel from its first block.
args <- commandArgs(TRUE)
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages({ library(EMC2, lib.loc = lib); library(posterior) })
cat("EMC2 from", find.package("EMC2"), "\n")
ckpt <- args[1]; label <- args[2]; seed <- as.integer(args[3]); iters <- as.integer(args[4]); out <- args[5]
ns <- asNamespace("EMC2")
`%||%` <- function(a, b) if (is.null(a)) b else a
cores_per_chain <- as.integer(Sys.getenv("CORES_PER_CHAIN", "9"))
calibrate <- Sys.getenv("CALIBRATE", "")
calibrate <- if (nzchar(calibrate)) as.integer(strsplit(calibrate, ",")[[1]]) else integer(0)
stopifnot(all(is.finite(calibrate)), all(calibrate > 0), !is.unsorted(calibrate, strictly = TRUE))
discard <- Sys.getenv("DISCARD", "")
discard <- if (nzchar(discard)) as.integer(discard) else if (length(calibrate)) max(calibrate) else 100L
n_particles <- Sys.getenv("N_PARTICLES", "")
arm_file <- Sys.getenv("ARM_R", "")
# PROFILE=1 records fit()'s per-iteration and per-block timings (see fit_arm.R).
if (Sys.getenv("PROFILE") == "1") options(emc2.sampler_profile = TRUE)
arm_options <- function() list()
arm_hook <- function(emc) emc
if (nzchar(arm_file)) source(arm_file, local = TRUE)
opts <- arm_options()
if (length(opts)) do.call(options, opts)

emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(ckpt)), quiet = TRUE)
class(emc) <- "emc"
n_subjects <- emc[[1]]$n_subjects
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
sample_idx <- which(emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)] == "sample")
stopifnot(length(sample_idx) >= 100)

# Scoring directions, fixed by the checkpoint: the group means and log sds
# standardised by their checkpoint sds, rotated to the directions of largest
# lag-one autocorrelation (the generalised eigenvectors of the lag-one
# autocovariance against the covariance, chains pooled). Cached beside the
# checkpoint so every arm scores on the same directions.
group_draws <- function(e, idx) lapply(e, function(ch) {
  mu <- t(ch$samples$theta_mu[, idx, drop = FALSE])
  lsd <- t(0.5 * log(apply(ch$samples$theta_var[, , idx, drop = FALSE], 3, diag)))
  cbind(mu, lsd)
})
dirs_file <- sub("\\.rds$", "_dirs.rds", ckpt)
if (file.exists(dirs_file)) dirs <- readRDS(dirs_file) else {
  G <- group_draws(emc, sample_idx)
  centre <- Reduce(`+`, lapply(G, colMeans)) / length(G)
  scale <- sqrt(Reduce(`+`, lapply(G, function(g) apply(g, 2, var))) / length(G))
  Z <- lapply(G, function(g) sweep(sweep(g, 2, centre), 2, scale, "/"))
  C <- Reduce(`+`, lapply(Z, function(z) crossprod(z) / nrow(z))) / length(Z)
  C1 <- Reduce(`+`, lapply(Z, function(z) { n <- nrow(z); crossprod(z[-n, ], z[-1, ]) / (n - 1) })) / length(Z)
  C1 <- (C1 + t(C1)) / 2
  R <- chol(C + diag(1e-10 * mean(diag(C)), ncol(C)))
  M <- backsolve(R, t(backsolve(R, C1, transpose = TRUE)), transpose = TRUE)
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
  U <- backsolve(R, ev$vectors)
  dirs <- list(centre = centre, scale = scale, U = U, rho = ev$values,
               names = c(paste0("mu.", pn), paste0("lsd.", pn)))
  saveRDS(dirs, dirs_file)
}

base <- local({ set.seed(seed, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) {
  base <- parallel::nextRNGStream(base)
  emc[[j]]$rng <- ns$.emc_chain_rng(base, n_subjects)
  emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions
  if (nzchar(n_particles)) {
    pm <- attr(emc[[j]]$samples, "pm_settings")
    for (s in seq_along(pm)) for (k in seq_along(pm[[s]])) if (!is.null(pm[[s]][[k]]$n_particles)) {
      pm[[s]][[k]]$n_particles <- as.integer(n_particles)
    }
    attr(emc[[j]]$samples, "pm_settings") <- pm
  }
}
if (length(calibrate)) {
  n_tr <- emc[[1]]$sample_transitions
  options(emc2.sample_rebuilds = c(n_tr + calibrate, 2L * (n_tr + max(calibrate))))
  cat("Calibration rebuilds at transitions", n_tr + calibrate, "then", 2L * (n_tr + max(calibrate)), "\n")
}
particle_factor <- if (nzchar(n_particles)) as.numeric(n_particles) / sqrt(p) else 50
emc <- arm_hook(emc)
n0 <- emc[[1]]$samples$idx
totals0 <- lapply(emc, function(ch) list(hmc = ch$subject_hmc$totals$attempts %||% 0,
  shear = ch$shear$rounds %||% 0, collapsed = ch$collapsed$rounds %||% 0))

cg_file <- file.path("/sys/fs/cgroup", sub("^0::", "", readLines("/proc/self/cgroup")[1]), "cpu.stat")
cg_cpu <- function() {
  x <- tryCatch(readLines(cg_file), error = function(e) character(0))
  u <- grep("^usage_usec", x, value = TRUE)
  if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_
}
c0 <- cg_cpu(); t0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample",
  stop_criteria = list(sample = list(iter = iters)), step_size = 100L, max_tries = 1L,
  cores_for_chains = length(emc), cores_per_chain = cores_per_chain, trim = FALSE,
  particle_factor = particle_factor, verbose = FALSE, verboseProgress = FALSE)))
el <- proc.time() - t0
cpu <- cg_cpu() - c0
if (!is.finite(cpu)) cpu <- el[["user.self"]] + el[["user.child"]] + el[["sys.self"]] + el[["sys.child"]]
emc <- ns$restore_duplicates(emc)
n1 <- emc[[1]]$samples$idx
new <- (n0 + 1):n1
stopifnot(length(new) == iters)
score <- new[-seq_len(discard)]

arr3 <- function(mats) { a <- simplify2array(mats); aperm(a, c(1, 3, 2)) }  # iterations x chains x quantities
diag_of <- function(a) list(ess_bulk = apply(a, 3, ess_bulk), ess_tail = apply(a, 3, ess_tail), rhat = apply(a, 3, rhat))
G <- group_draws(emc, score)
Zd <- lapply(G, function(g) sweep(sweep(g, 2, dirs$centre), 2, dirs$scale, "/") %*% dirs$U)
d_dirs <- diag_of(arr3(Zd))
d_mu <- diag_of(arr3(lapply(G, function(g) g[, seq_len(p), drop = FALSE])))
d_lsd <- diag_of(arr3(lapply(G, function(g) g[, p + seq_len(p), drop = FALSE])))
d_alpha <- diag_of(arr3(lapply(emc, function(ch) t(matrix(ch$samples$alpha[, , score], p * n_subjects)))))
halves <- split(score, rep(1:2, each = ceiling(length(score) / 2))[seq_along(score)])
half_dirs <- lapply(halves, function(h) {
  Gh <- group_draws(emc, h)
  apply(arr3(lapply(Gh, function(g) sweep(sweep(g, 2, dirs$centre), 2, dirs$scale, "/") %*% dirs$U)), 3, ess_bulk)
})
totals1 <- lapply(emc, function(ch) list(hmc = ch$subject_hmc$totals$attempts %||% 0,
  shear = ch$shear$rounds %||% 0, collapsed = ch$collapsed$rounds %||% 0))
moves <- mapply(function(a, b) unlist(b) - unlist(a), totals0, totals1, SIMPLIFY = TRUE)
pm <- attr(emc[[1]]$samples, "pm_settings")
np <- unlist(lapply(pm, function(s) lapply(s, function(k) k$n_particles)))
last_plan <- lapply(emc, function(ch) { h <- ch$alloc$history; if (length(h)) h[[length(h)]] else NULL })
res <- list(ckpt = ckpt, label = label, seed = seed, iters = iters, discard = discard, scored = length(score),
            calibrate = calibrate, alloc_history = lapply(emc, function(ch) ch$alloc$history),
            lib = find.package("EMC2"), cpu = cpu, elapsed = el[["elapsed"]],
            cpu_per_iter_chain = cpu / iters / length(emc), n_particles = range(np),
            dirs = d_dirs, mu = d_mu, lsd = d_lsd, alpha = d_alpha, dir_rho = dirs$rho,
            half_dirs = half_dirs, moves = moves,
            plan = lapply(last_plan, function(l) if (is.null(l)) NULL else
              list(rounds = l$rounds, rate = l$rate, particle_rate = l$particle_rate, tau = l$tau, tau0 = l$tau0)),
            n_subjects = n_subjects, p = p,
            blocks = if (Sys.getenv("PROFILE") == "1") EMC2:::.emc_block_profile() else NULL,
            iterations = if (Sys.getenv("PROFILE") == "1") EMC2:::.emc_iteration_profile() else NULL)
saveRDS(res, out)
k_dir <- which.min(d_dirs$ess_bulk)
cat(sprintf("%s seed %d: %d iterations (%d scored) | cpu %.0f s (%.2f s/iter/chain) elapsed %.0f s | particles %s\n",
  label, seed, iters, length(score), cpu, res$cpu_per_iter_chain, res$elapsed, paste(unique(range(np)), collapse = "-")))
cat(sprintf("  slowest direction %d: ESS %.0f (halves %.0f / %.0f), %.2f per kCPU-s | direction ESS min/median %.0f / %.0f\n",
  k_dir, d_dirs$ess_bulk[k_dir], half_dirs[[1]][k_dir], half_dirs[[2]][k_dir], 1000 * d_dirs$ess_bulk[k_dir] / cpu,
  min(d_dirs$ess_bulk), median(d_dirs$ess_bulk)))
for (nm in c("mu", "lsd", "alpha")) { d <- res[[nm]]
  cat(sprintf("  %-5s ESS bulk min %.0f median %.0f  tail min %.0f | min per kCPU-s %.2f | max R-hat %.3f\n",
    nm, min(d$ess_bulk), median(d$ess_bulk), min(d$ess_tail), 1000 * min(d$ess_bulk) / cpu, max(d$rhat))) }
cat("  moves per chain (hmc attempts, shear rounds, collapsed rounds):\n"); print(moves)
