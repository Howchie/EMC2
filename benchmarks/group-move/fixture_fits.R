# First-pass fit test for the ensemble group update on data from
# WorkingTests/generate_correlated_group_rdm_fixture.R (RDM with A = 0, 13 coefficients, match drift ratio 1.5, 100 trials per cell, group SDs x3 so every SD is >= 1.4x the per-subject SE).
#
#   Rscript fixture_fits.R fit <lib> <case> <arm> <seed>
#   Rscript fixture_fits.R summarize
#
# <arm> is "ensemble" (emc2.ensemble = TRUE) or "baseline" (FALSE); both use the
# same build <lib>, so the only difference is the ensemble update. Package
# defaults for preburn/burn/adapt, then 1,500 fixed production iterations,
# 3 chains x 9 cores. Fixtures are generated once per case by the fixture
# script with the case's settings and cached.
args <- commandArgs(TRUE)
root <- Sys.getenv("FIXTURE_ROOT", "benchmarks/group-move/results/fixture_A0v15")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
# Each case: generator settings (gen), optionally a fixture shared with another
# case (fixture), and fit settings (type, par_groups, sample iterations).
gen_base <- c(CORR_TRIALS_PER_CELL = "100", CORR_SD_SCALE = "3")
cases <- list(
  diag   = list(gen = c(gen_base, CORR_STRENGTH = "0")),
  corr   = list(gen = c(gen_base, CORR_STRENGTH = "1")),
  strong = list(gen = c(gen_base, CORR_STRENGTH = "1.2")),
  # 80 subjects: pool and reselection cost grow with n.
  n80    = list(gen = c(gen_base, CORR_STRENGTH = "1", CORR_N_SUBJECTS = "80")),
  # Free A: per-subject A/B/t0 ridge. Longer production so both arms converge.
  freeA  = list(gen = c(gen_base, CORR_STRENGTH = "1", CORR_FREE_A = "1"), sample = 4000L),
  # Other hierarchies on existing fixtures.
  corr_blocked = list(fixture = "corr", blocks = TRUE),
  diag_dg      = list(fixture = "diag", type = "diagonal-gamma"),
  # Latent-factor hierarchies on the correlated fixture.
  corr_factor = list(fixture = "corr", type = "factor", n_factors = 2L),
  corr_infnt  = list(fixture = "corr", type = "infnt_factor", n_factors = 4L),
  # Long baseline-only run on the strongest case: settles which arm is off
  # when the two arms' posterior means differ by more than Monte Carlo error.
  strong_long = list(fixture = "strong", sample = 12000L),
  # Sanity check, the case where the ensemble should be least useful: no
  # group correlation and 4x the data per subject, so subject posteriors are
  # tight, the Gibbs step is already close to independent, and every pool
  # member costs a 4x-longer likelihood.
  tight  = list(gen = c(gen_base[-1], CORR_TRIALS_PER_CELL = "400", CORR_STRENGTH = "0")),
  # Second sanity check: wide subjects (24 parameters), no correlation. The
  # pool proposal's weights degenerate with dimension, so the ensemble should
  # rarely move anyone while still paying for its pools.
  wide   = list(wide = TRUE),
  # Real data, same property: 19 subjects with ~830 trials each.
  forstmann = list(data = "forstmann")
)
fixture_file <- function(case) file.path(root, paste0("fixture-",
  if (is.null(cases[[case]]$fixture)) case else cases[[case]]$fixture, ".rds"))

if (args[1] == "fit") {
  lib <- args[2]; case <- args[3]; arm <- args[4]; seed <- as.integer(args[5])
  cs <- cases[[case]]
  ff <- fixture_file(case)
  if (isTRUE(cs$wide) && !file.exists(ff)) {
    suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
    fac <- list(subjects = seq_len(30), E = c("speed", "neutral", "accuracy"),
                S = c("left", "right"))
    # Every group SD >= 1.4x its per-subject SE (a 27-parameter version with
    # v ~ lM * E * S at SD 0.3 put the v x S contrasts at 0.55-0.78x and no
    # arm converged).
    fo <- list(v ~ lM * E, B ~ E * lR * S, t0 ~ E * S)
    co <- c(s = log(1), A = log(0))
    dw <- suppressMessages(design(factors = fac, Rlevels = c("left", "right"),
      matchfun = function(d) d$S == d$lR, model = RDM, formula = fo, constants = co))
    gm <- sampled_pars(dw); gm[] <- 0
    gm[c("v", "v_lMTRUE", "B", "t0")] <- c(log(2), log(1.5), log(1.1), log(0.25))
    cv <- diag(0.5^2, length(gm)); dimnames(cv) <- list(names(gm), names(gm))
    set.seed(20260929)
    se <- make_random_effects(dw, gm, n_subj = 30, covariances = cv)
    dat <- make_data(se, dw, n_trials = 200, verbose = FALSE)
    saveRDS(list(data = dat, group_means = gm, intended_covariance = cv, subject_effects = se,
                 settings = list(formula = fo, constants = co)), ff)
  }
  if (is.null(cs$data) && !file.exists(ff)) {
    gen <- if (is.null(cs$fixture)) cs$gen else cases[[cs$fixture]]$gen
    env <- c(gen, CORR_FIXTURE_OUT = ff, EMC_LIB = lib, CORR_CHECK_N = "5000")
    st <- system2("Rscript", "WorkingTests/generate_correlated_group_rdm_fixture.R",
                  env = paste0(names(env), "=", env))
    if (st != 0 || !file.exists(ff)) stop("fixture generation failed")
  }
  suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
  cat("EMC2 from", find.package("EMC2"), "\n")
  fx <- if (is.null(cs$data)) readRDS(ff) else
    list(data = get(cs$data, asNamespace("EMC2")),
         settings = list(formula = list(v ~ lM * E, B ~ E * lR, t0 ~ 1),
                         constants = c(s = log(1), A = log(0))))
  options(emc2.ensemble = arm %in% c("ensemble", "ensemble_marg"),
          emc2.ensemble_latent_weight = if (arm == "ensemble_marg") "marginal" else "conditional")
  # Extra internal options for design experiments, e.g.
  # EMC2_OPTS="emc2.ensemble_rounds=5,emc2.ensemble_pool=16".
  for (kv in strsplit(strsplit(Sys.getenv("EMC2_OPTS"), ",")[[1]], "=")) {
    if (length(kv) == 2L) options(stats::setNames(list(type.convert(kv[2], as.is = TRUE)), kv[1]))
  }
  des <- suppressMessages(design(data = fx$data, model = RDM,
    matchfun = function(d) d$S == d$lR, formula = fx$settings$formula,
    constants = fx$settings$constants))
  type <- if (is.null(cs$type)) "standard" else cs$type
  pn0 <- names(sampled_pars(des))
  par_groups <- if (isTRUE(cs$blocks)) match(sub("_.*", "", pn0), unique(sub("_.*", "", pn0))) else NULL
  set.seed(1000L * seed + 3L)
  emc <- suppressMessages(do.call(make_emc, c(list(fx$data, des, type = type, n_chains = 3L,
                                   rt_resolution = NULL, par_groups = par_groups),
                                   if (!is.null(cs$n_factors)) list(n_factors = cs$n_factors))))
  n_sample <- if (is.null(cs$sample)) 1500L else cs$sample
  timing <- c()
  for (stage in c("preburn", "burn", "adapt", "sample")) {
    crit <- EMC2:::get_stop_criteria(stage, if (stage == "sample") list(iter = n_sample) else NULL,
                                     type)
    t0 <- proc.time()[["elapsed"]]
    emc <- run_emc(emc, stage = stage, stop_criteria = crit, cores_for_chains = 3L,
                   cores_per_chain = as.integer(Sys.getenv("CORES_PER_CHAIN", "9")),
                   particle_factor = as.numeric(Sys.getenv("PARTICLE_FACTOR", "50")),
                   verbose = FALSE)
    timing[stage] <- proc.time()[["elapsed"]] - t0
    cat(format(Sys.time(), "%T"), case, arm, seed, stage, round(timing[stage]), "s\n")
  }
  emc <- EMC2:::restore_duplicates(emc)
  iters <- sapply(c("preburn", "burn", "adapt", "sample"), function(st)
    sum(emc[[1]]$samples$stage == st))
  idx <- which(emc[[1]]$samples$stage == "sample")
  pn <- rownames(emc[[1]]$samples$theta_mu)
  grab <- function(f, nm) { out <- array(NA_real_, c(length(idx), length(emc), length(nm)),
    dimnames = list(NULL, NULL, nm)); for (ch in seq_along(emc)) out[, ch, ] <- t(f(emc[[ch]]$samples)); out }
  pr <- which(upper.tri(diag(length(pn))), arr.ind = TRUE)
  draws <- list(
    mu = grab(function(s) s$theta_mu[, idx], pn),
    logsd = grab(function(s) 0.5 * log(apply(s$theta_var[, , idx], 3, diag)), pn),
    cor = if (type == "diagonal-gamma") NULL else grab(function(s) apply(s$theta_var[, , idx], 3, function(v) cov2cor(v)[pr]),
               paste(pn[pr[, 1]], pn[pr[, 2]], sep = "~")))
  alpha_ess <- sapply(seq_len(dim(emc[[1]]$samples$alpha)[2]), function(sj)
    sapply(seq_along(pn), function(p) posterior::ess_bulk(
      sapply(emc, function(ch) ch$samples$alpha[p, sj, idx]))))
  ens <- tryCatch(EMC2:::ensemble_diagnostics(emc), error = function(e) NULL)
  ens_subj <- lapply(emc, function(ch) { st <- ch$ensemble_stats
    if (is.null(st$subj_moved)) NULL else cbind(moved = st$subj_moved, wess = st$subj_wess) / st$sweeps })
  alpha_draws <- if (Sys.getenv("SAVE_ALPHA") == "1" || case == "tight")
    lapply(emc, function(ch) ch$samples$alpha[, , idx, drop = FALSE])
  saveRDS(list(case = case, arm = arm, seed = seed, timing = timing,
               particle_factor = as.numeric(Sys.getenv("PARTICLE_FACTOR", "50")),
               emc2_opts = Sys.getenv("EMC2_OPTS"),
               iters = iters,
               draws = draws, alpha_ess = alpha_ess, ens = ens, ens_subj = ens_subj,
               alpha_draws = alpha_draws,
               type = type, par_groups = par_groups,
               truth = list(mu = fx$group_means, cov = fx$intended_covariance)),
          file.path(root, sprintf("%s-%s-seed%d.rds", case, arm, seed)))
  if (nzchar(Sys.getenv("SAVE_EMC")))
    save(emc, file = file.path(root, sprintf("emc-%s-%s-seed%d.RData", case, arm, seed)))
}

if (args[1] == "summarize") {
  fs <- list.files(root, pattern = "-seed[0-9]+\\.rds$", full.names = TRUE)
  res <- lapply(fs, readRDS)
  # Drop structurally constant columns (zero correlations across covariance
  # blocks or in the diagonal hierarchy).
  live <- function(draws) lapply(Filter(Negate(is.null), draws), function(z)
    z[, , apply(z, 3, stats::sd) > 0, drop = FALSE])
  rows <- do.call(rbind, lapply(res, function(r) {
    dr <- live(r$draws)
    e <- lapply(dr, function(z) apply(z, 3, posterior::ess_bulk))
    h <- unlist(lapply(dr, function(z) apply(z, 3, posterior::rhat)))
    mn <- function(x) if (length(x)) round(min(x)) else NA
    data.frame(case = r$case, arm = r$arm, seed = r$seed,
      t_pre = round(sum(r$timing[c("preburn", "burn", "adapt")])), t_sample = round(r$timing[["sample"]]),
      max_rhat = round(max(h), 3), n_rhat_gt_1.1 = sum(h > 1.1),
      min_mu = mn(e$mu), med_mu = round(median(e$mu)),
      min_sd = mn(e$logsd), min_cor = mn(e$cor),
      med_cor = if (length(e$cor)) round(median(e$cor)) else NA,
      worst = mn(unlist(e)),
      min_alpha = round(min(r$alpha_ess)), med_alpha = round(median(r$alpha_ess)),
      moved = if (!is.null(r$ens)) round(mean(r$ens$moved), 2) else NA,
      wess = if (!is.null(r$ens)) round(mean(r$ens$weight_ess), 3) else NA)
  }))
  rows <- rows[order(rows$case, rows$seed, rows$arm), ]
  print(rows, row.names = FALSE)
  cat("\nEnsemble / baseline, paired by seed. ESS per second of production on the worst",
      "\ngroup quantity and on the median subject parameter; total wall time; and",
      "\nagreement of posterior means (z = difference / combined MCSE over every live",
      "\ngroup quantity: sd(z) ~ 1 and ~5% beyond 2 if both arms target the same posterior).\n")
  key <- function(r) paste(r$case, r$seed)
  by_key <- split(res, vapply(res, key, ""))
  for (k in sort(names(by_key))) for (arm in c("ensemble", "ensemble_marg")) {
    g <- by_key[[k]]
    a <- Filter(function(r) r$arm == arm, g); b <- Filter(function(r) r$arm == "baseline", g)
    if (!length(a) || !length(b)) next
    a <- a[[1]]; b <- b[[1]]
    ra <- rows[rows$case == a$case & rows$seed == a$seed & rows$arm == arm, ]
    rb <- rows[rows$case == b$case & rows$seed == b$seed & rows$arm == "baseline", ]
    da <- live(a$draws); db <- live(b$draws)
    z <- unlist(lapply(names(da), function(bl) {
      nm <- intersect(dimnames(da[[bl]])[[3]], dimnames(db[[bl]])[[3]])
      vapply(nm, function(p) { x <- da[[bl]][, , p]; y <- db[[bl]][, , p]
        (mean(x) - mean(y)) / sqrt(posterior::mcse_mean(x)^2 + posterior::mcse_mean(y)^2) }, 0)
    }))
    cat(sprintf("  %-16s %-13s worst %5.2f  subj med %5.2f min %5.2f  time %4.2f  | z: sd %.2f, >2 %4.1f%%, max %.1f (n=%d)\n",
      k, arm, (ra$worst / ra$t_sample) / (rb$worst / rb$t_sample),
      (ra$med_alpha / ra$t_sample) / (rb$med_alpha / rb$t_sample),
      (ra$min_alpha / ra$t_sample) / (rb$min_alpha / rb$t_sample),
      (ra$t_pre + ra$t_sample) / (rb$t_pre + rb$t_sample),
      stats::sd(z), 100 * mean(abs(z) > 2), max(abs(z)), length(z)))
  }
}
