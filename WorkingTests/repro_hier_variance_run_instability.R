## Minimal reproduction: between-run instability of the hierarchical variance
## posterior (and hence pWAIC) for an unchanged model and data.
##
## Observed in an applied project (2026-09-25): two fits with identical data,
## priors, bounds and likelihood (RDMSWTN vs RDMSWTN_UT with u fixed at 0; the
## summed log-likelihood matched to machine precision at the same parameters)
## gave the same group means but systematically different between-subject SDs
## (every parameter larger in one run, e.g. 0.74 vs 0.89) and pWAIC 468 vs 556.
## Each run looked converged: max R-hat on sigma2 <= 1.02, min ESS ~500.
## Between-run disagreement far beyond within-run Monte Carlo error suggests a
## seeding, chain-initialisation or variance-mixing problem, not model
## differences.
##
## This script uses no project data or custom kernels. It simulates a plain
## hierarchical RDMSWTN and fits it under controlled conditions:
##   A  seed 1                      (reference)
##   B  seed 1 again                (must be bit-identical to A if seeding is
##                                   reproducible; with parallel chains under
##                                   Mersenne-Twister it may legitimately not be)
##   C  seed 2, D  seed 3           (between-run agreement within MCSE expected)
##   E  seed 1, RDMSWTN_UT(u = 0)   (identical model to A: should agree like C/D)
##   S1, S2  seed 1 twice with serial chains (cores_for_chains = 1): isolates
##          parallel RNG handling from sampler behaviour
##
## Environment knobs (defaults are small but not trivial):
##   REPRO_N_SUBJ (30), REPRO_N_TRIALS (100 per design cell), REPRO_CPC
##   (cores per chain, 2), REPRO_RUNS ("A,B,C,D,E,S1,S2"), REPRO_OUT (rds path)
## Run: Rscript WorkingTests/repro_hier_variance_run_instability.R
suppressMessages({library(EMC2)})

n_subj   <- as.integer(Sys.getenv("REPRO_N_SUBJ", "30"))
n_trials <- as.integer(Sys.getenv("REPRO_N_TRIALS", "100"))
cpc      <- as.integer(Sys.getenv("REPRO_CPC", "2"))
runs     <- strsplit(Sys.getenv("REPRO_RUNS", "A,B,C,D,E,S1,S2"), ",")[[1]]
out_file <- Sys.getenv("REPRO_OUT", "repro_hier_variance_run_instability.rds")
cat("RNGkind:", paste(RNGkind(), collapse = " / "), "\n")

## ---- Simulated data: 2-choice RDMSWTN, match effect on v, A = 0, s = 1 ----
matchfun <- function(d) d$S == d$lR
mk_design <- function(model) design(
  factors = list(subjects = seq_len(n_subj), S = c("left", "right")),
  Rlevels = c("left", "right"), matchfun = matchfun,
  formula = list(v ~ lM, B ~ 1, t0 ~ 1),
  constants = c(s = log(1)),
  model = model)
des_rdm <- mk_design(RDMSWTN)
des_ut  <- mk_design(function() RDMSWTN_UT())  # u defaults to log(0) = 0

group_means <- sampled_pars(des_rdm)
group_means[] <- c(log(1), log(2), log(1.5), log(0.25))   # v, v_lMTRUE, B, t0
print(group_means)
set.seed(20260925)
subj_pars <- make_random_effects(des_rdm, group_means, n_subj = n_subj,
                                 variance_proportion = 0.2)
dat <- make_data(subj_pars, des_rdm, n_trials = n_trials)
cat("data:", nrow(dat), "trials,", n_subj, "subjects\n")
true_sd <- apply(subj_pars, 2, sd)

## ---- One fit ------------------------------------------------------------
fit_one <- function(label, seed, des, serial = FALSE) {
  set.seed(seed)
  emc <- make_emc(dat, des, n_chains = 3)
  t_start <- Sys.time()
  emc <- fit(emc, cores_per_chain = cpc,
             cores_for_chains = if (serial) 1 else 3, verbose = FALSE)
  secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  mu  <- get_pars(emc, selection = "mu", stage = "sample",
                  merge_chains = TRUE, return_mcmc = FALSE)
  s2  <- get_pars(emc, selection = "sigma2", stage = "sample",
                  merge_chains = TRUE, return_mcmc = FALSE)
  sdd <- if (length(dim(s2)) == 3) sapply(dimnames(s2)[[1]], function(p) sqrt(s2[p, p, ]))
         else t(sqrt(s2))                     # iterations x parameters
  ## ESS of each SD computed directly from the merged draws (the names from
  ## ess_summary() do not line up with these columns).
  ess_sd <- coda::effectiveSize(coda::as.mcmc(sdd))
  ll  <- EMC2:::.ll_matrix_pooled(emc, stage = "sample", filter = 0, cores = cpc)
  lpd <- sum(matrixStats::colLogSumExps(ll) - log(nrow(ll)))
  pw  <- sum(matrixStats::colVars(ll))
  ## First retained sigma2 draw per chain: are chains starting in the same place?
  s2c <- get_pars(emc, selection = "sigma2", stage = "sample",
                  merge_chains = FALSE, return_mcmc = FALSE)
  first <- tryCatch(sapply(s2c, function(ch) {
    if (length(dim(ch)) == 3) sqrt(diag(ch[, , 1])) else sqrt(ch[, 1]) }),
    error = function(e) NULL)
  list(label = label, seed = seed, serial = serial, model = des$model()$c_name,
       secs = secs, mu_draws = mu, sd_draws = sdd, ess_sd = ess_sd,
       rhat_mu = max(unlist(gd_summary(emc, selection = "mu")), na.rm = TRUE),
       rhat_sd = max(unlist(gd_summary(emc, selection = "sigma2")), na.rm = TRUE),
       lpd = lpd, pwaic = pw, first_sd_by_chain = first)
}

spec <- list(A = list(1, des_rdm, FALSE), B = list(1, des_rdm, FALSE),
             C = list(2, des_rdm, FALSE), D = list(3, des_rdm, FALSE),
             E = list(1, des_ut, FALSE),
             S1 = list(1, des_rdm, TRUE), S2 = list(1, des_rdm, TRUE))
res <- list()
for (r in runs) {
  cat("\n== run", r, "\n")
  res[[r]] <- do.call(fit_one, c(list(label = r), setNames(spec[[r]], c("seed", "des", "serial"))))
  with(res[[r]], cat(sprintf("  %s seed %d serial %s: %.0fs  Rhat mu %.3f sd %.3f  lpd %.2f  pWAIC %.2f\n",
                             model, seed, serial, secs, rhat_mu, rhat_sd, lpd, pwaic)))
}
saveRDS(list(res = res, true_sd = true_sd, group_means = group_means), out_file)

## ---- Diagnostics -----------------------------------------------------------
sd_med <- sapply(res, function(x) apply(x$sd_draws, 2, median))
mcse   <- sapply(res, function(x) apply(x$sd_draws, 2, sd) / sqrt(pmax(x$ess_sd, 1)))
cat("\nPosterior median between-subject SD by run (true sample SD in first column):\n")
print(round(cbind(true = true_sd[rownames(sd_med)], sd_med), 3))
cat("\nlpd / pWAIC by run:\n")
print(round(sapply(res, function(x) c(lpd = x$lpd, pWAIC = x$pwaic)), 2))

pair_z <- function(a, b) (sd_med[, a] - sd_med[, b]) / sqrt(mcse[, a]^2 + mcse[, b]^2)
cat("\nBetween-run SD differences in MCSE units (|z| > 4 is not Monte Carlo noise):\n")
pairs <- combn(names(res), 2)
zs <- apply(pairs, 2, function(p) max(abs(pair_z(p[1], p[2]))))
print(data.frame(run1 = pairs[1, ], run2 = pairs[2, ], max_abs_z = round(zs, 1)), row.names = FALSE)

same <- function(a, b) if (all(c(a, b) %in% names(res)))
  isTRUE(all.equal(res[[a]]$mu_draws, res[[b]]$mu_draws, tolerance = 0)) else NA
cat("\nReproducibility (identical mu draws): A==B", same("A", "B"),
    "  S1==S2", same("S1", "S2"), "\n")
cat("First sigma draw per chain (identical columns => shared initialisation):\n")
for (r in names(res)) if (!is.null(res[[r]]$first_sd_by_chain)) {
  cat(" ", r, "\n"); print(round(res[[r]]$first_sd_by_chain, 3)) }

cat("\nFlags:\n")
if (isTRUE(any(zs > 4))) cat("  * Between-run SD posteriors disagree beyond MCSE (max |z| =",
                     round(max(zs), 1), ") -> variance posterior not reproducible across runs.\n")
if (isFALSE(same("S1", "S2"))) cat("  * Serial chains with the same seed are NOT reproducible -> seeding bug.\n")
if (isFALSE(same("A", "B"))) cat("  * Parallel chains with the same seed differ (expected under Mersenne-Twister + mclapply; check whether EMC2 intends otherwise).\n")
