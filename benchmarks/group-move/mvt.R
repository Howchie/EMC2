# Minimal viable test: does the AM group move beat the ORIGINAL sampler
# (commit 940d8b40, before any group proposal existed)?
#
#   Rscript mvt.R data                      # simulate the shared datasets once
#   Rscript mvt.R fit <lib> <arm> <scenario> <seed> <cores_per_chain>
#   Rscript mvt.R summarize
#
# <lib> is an R library holding the EMC2 build to test: 940d8b40 for
# arm=original, the working tree for arm=am. Each fit is its own process so the
# two builds never share a session. Run original and am of one seed
# concurrently with equal cores so they see the same machine load.
#
# Scenarios (RDM, 6 group parameters, 80 subjects, standard hierarchy unless
# noted):
#   hard     40 trials/subject, group correlation 0.8: subject data are weak
#            relative to between-subject spread, the regime a group move targets
#   easy     300 trials/subject, group correlation 0.3: data dominate, group
#            move should be neutral; tests for regression
#   easy_dg  as easy, diagonal-gamma hierarchy

args <- commandArgs(TRUE)
mode <- args[1]
root <- "benchmarks/group-move/results/mvt"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

scenarios <- list(
  hard    = list(trials = 40L,  rho = 0.8, type = "standard"),
  easy    = list(trials = 300L, rho = 0.3, type = "standard"),
  easy_dg = list(trials = 300L, rho = 0.3, type = "diagonal-gamma")
)
truth <- c(v = log(2.5), v_lMd = 1, B = log(1.5), B_lRd = 0,
           A = log(0.3), t0 = log(0.2))
true_sd <- sqrt(log1p(0.2^2))   # make_random_effects' default spread for log pars
n_subjects <- 80L

make_des <- function(data = NULL) {
  AD <- matrix(c(-0.5, 0.5), ncol = 1L, dimnames = list(NULL, "d"))
  common <- list(model = EMC2::RDM, matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, B ~ lR, A ~ 1, t0 ~ 1, s ~ 1),
    contrasts = list(v = list(lM = AD), B = list(lR = AD)),
    constants = c(s = log(1)))
  suppressMessages(if (is.null(data))
    do.call(EMC2::design, c(common, list(factors = list(subjects = seq_len(n_subjects),
      S = c("left", "right")), Rlevels = c("left", "right"), report_p_vector = FALSE)))
  else do.call(EMC2::design, c(common, list(data = data, report_p_vector = FALSE))))
}

if (mode == "data") {
  suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
  for (sc in names(scenarios)) {
    cfg <- scenarios[[sc]]
    key <- sprintf("trials%d-rho%.1f", cfg$trials, cfg$rho)
    f <- file.path(root, paste0("data-", key, ".rds"))
    if (file.exists(f)) next
    des <- make_des()
    R <- matrix(cfg$rho, 6, 6); diag(R) <- 1
    Sigma <- true_sd^2 * R; dimnames(Sigma) <- list(names(truth), names(truth))
    set.seed(4242L + cfg$trials)
    subj <- make_random_effects(des, truth, n_subj = n_subjects, covariances = Sigma)
    dat <- make_data(subj, des, n_trials = cfg$trials %/% 2L)  # 2 S cells
    saveRDS(list(data = dat, truth = truth, Sigma = Sigma, subj = subj), f)
    cat("wrote", f, nrow(dat), "rows\n")
  }
}

if (mode == "fit") {
  lib <- args[2]; arm <- args[3]; sc <- args[4]; seed <- as.integer(args[5])
  cpc <- as.integer(args[6])
  suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
  cat("EMC2 loaded from", find.package("EMC2"), "\n")
  cfg <- scenarios[[sc]]
  key <- sprintf("trials%d-rho%.1f", cfg$trials, cfg$rho)
  sim <- readRDS(file.path(root, paste0("data-", key, ".rds")))
  if (arm == "am") options(emc2.group_move = TRUE, emc2.group_move_method = "am",
                           emc2.group_move_proposals = 1L, emc2.group_move_scale = TRUE)
  des <- make_des(sim$data)
  set.seed(1000L * seed + 7L)
  emc <- suppressMessages(make_emc(sim$data, des, type = cfg$type, n_chains = 3L,
                                   rt_resolution = NULL))
  timing <- c(); iters <- c()
  crit <- list(preburn = list(iter = 100L), burn = NULL, adapt = NULL,
               sample = list(iter = 2000L))
  for (stage in names(crit)) {
    sc_crit <- EMC2:::get_stop_criteria(stage, crit[[stage]], cfg$type)
    t0 <- proc.time()[["elapsed"]]
    emc <- suppressWarnings(run_emc(emc, stage = stage, stop_criteria = sc_crit,
      cores_for_chains = 3L, cores_per_chain = cpc, max_tries = 20L,
      step_size = 100L, verbose = FALSE, trim = FALSE))
    timing[stage] <- proc.time()[["elapsed"]] - t0
    cat(sprintf("%s %s %s seed%d: %s done, %.0fs\n", arm, sc, cfg$type, seed, stage, timing[stage]))
  }
  iters <- EMC2::chain_n(emc)[1L, ]
  idx <- which(emc[[1]]$samples$stage == "sample")
  pn <- rownames(emc[[1]]$samples$theta_mu)
  arr <- function(f, nm) {
    out <- array(NA_real_, c(length(idx), length(emc), length(nm)),
                 dimnames = list(NULL, NULL, nm))
    for (ch in seq_along(emc)) out[, ch, ] <- f(emc[[ch]]$samples)
    out
  }
  draws <- list(
    mu = arr(function(s) t(s$theta_mu[, idx]), paste0("mu_", pn)),
    logsd = arr(function(s) t(vapply(idx, function(i) 0.5 * log(diag(s$theta_var[, , i])),
                                     numeric(length(pn)))), paste0("logsd_", pn)))
  if (cfg$type == "standard") {
    pr <- which(upper.tri(diag(length(pn))), arr.ind = TRUE)
    draws$cor <- arr(function(s) t(vapply(idx, function(i)
      stats::cov2cor(s$theta_var[, , i])[pr], numeric(nrow(pr)))),
      paste0("cor_", pn[pr[, 1]], "__", pn[pr[, 2]]))
  }
  alpha_ess <- sapply(seq_len(dim(emc[[1]]$samples$alpha)[2]), function(sj)
    sapply(seq_along(pn), function(p) posterior::ess_bulk(
      sapply(emc, function(ch) ch$samples$alpha[p, sj, idx]))))
  gm <- tryCatch(EMC2::group_move_diagnostics(emc), error = function(e) NULL)
  saveRDS(list(arm = arm, scenario = sc, seed = seed, lib = find.package("EMC2"),
               timing = timing, iters = iters, draws = draws,
               alpha_ess = alpha_ess, gm = gm),
          file.path(root, sprintf("%s-%s-seed%d.rds", sc, arm, seed)))
}

if (mode == "summarize") {
  fs <- list.files(root, pattern = "^(hard|easy|easy_dg)-.*seed.*\\.rds$", full.names = TRUE)
  res <- lapply(fs, readRDS)
  fit_rows <- do.call(rbind, lapply(res, function(r) {
    ess <- unlist(lapply(r$draws, function(z) apply(z, 3, posterior::ess_bulk)))
    rh <- unlist(lapply(r$draws, function(z) apply(z, 3, posterior::rhat)))
    blk <- rep(names(r$draws), vapply(r$draws, function(z) dim(z)[3], numeric(1)))
    data.frame(scenario = r$scenario, arm = r$arm, seed = r$seed,
      burn_it = r$iters[["burn"]], adapt_it = r$iters[["adapt"]],
      t_pre_sample = sum(r$timing[c("preburn", "burn", "adapt")]),
      t_sample = r$timing[["sample"]],
      s_per_iter = r$timing[["sample"]] / 2000,
      max_rhat = max(rh), minESS_mu = min(ess[blk == "mu"]),
      minESS_logsd = min(ess[blk == "logsd"]),
      minESS_cor = if (any(blk == "cor")) min(ess[blk == "cor"]) else NA,
      medESS_alpha = median(r$alpha_ess),
      gm_acc = if (!is.null(r$gm)) mean(r$gm$sample_acceptance) else NA)
  }))
  fit_rows <- fit_rows[order(fit_rows$scenario, fit_rows$seed, fit_rows$arm), ]
  cat("== per fit ==\n"); print(fit_rows, digits = 3, row.names = FALSE)

  # Run-to-run reproducibility: for each variable, the spread of posterior means
  # across seeds in units of the posterior SD. A sampler that has really mixed
  # gives ~ sqrt(1/ESS) (a few %); a sampler that "converges" to a different
  # place each run gives a large value.
  cat("\n== cross-seed spread of posterior means / posterior sd (median, max over variables) ==\n")
  for (sc in unique(fit_rows$scenario)) for (arm in c("original", "am")) {
    rr <- Filter(function(r) r$scenario == sc && r$arm == arm, res)
    if (length(rr) < 2) next
    for (b in names(rr[[1]]$draws)) {
      m <- sapply(rr, function(r) apply(r$draws[[b]], 3, mean))
      s <- sapply(rr, function(r) apply(r$draws[[b]], 3, sd))
      z <- apply(m, 1, sd) / rowMeans(s)
      cat(sprintf("  %-8s %-9s %-6s median %.3f  max %.3f  (%s)\n", sc, arm, b,
                  median(z), max(z), names(z)[which.max(z)]))
    }
  }
  # Pooled agreement between arms (all seeds' chains pooled)
  cat("\n== pooled posterior agreement, am vs original: |mean diff| / pooled sd ==\n")
  for (sc in unique(fit_rows$scenario)) {
    get <- function(arm, b) do.call(abind3, lapply(Filter(function(r)
      r$scenario == sc && r$arm == arm, res), function(r) r$draws[[b]]))
    abind3 <- function(...) { x <- list(...); out <- array(NA, c(dim(x[[1]])[1],
      sum(sapply(x, function(a) dim(a)[2])), dim(x[[1]])[3]),
      dimnames = list(NULL, NULL, dimnames(x[[1]])[[3]])); k <- 0
      for (a in x) { out[, k + seq_len(dim(a)[2]), ] <- a; k <- k + dim(a)[2] }; out }
    for (b in names(res[[1]]$draws)) {
      if (!any(sapply(res, function(r) r$scenario == sc && !is.null(r$draws[[b]])))) next
      o <- get("original", b); a <- get("am", b)
      if (is.null(dim(o)) || is.null(dim(a))) next
      d <- abs(apply(a, 3, mean) - apply(o, 3, mean)) / apply(o, 3, sd)
      rs <- apply(a, 3, sd) / apply(o, 3, sd)
      cat(sprintf("  %-8s %-6s max|dz| %.3f   sd ratio %.2f-%.2f\n", sc, b, max(d),
                  min(rs), max(rs)))
    }
  }
}
