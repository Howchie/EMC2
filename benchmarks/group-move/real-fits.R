# Decision benchmark for the group-move default: real EMC2 fits, run the way a
# user runs them, comparing
#   original  - no group move (emc2.group_move = FALSE): the pre-group-proposal
#               sampler, i.e. the baseline a default change must beat
#   legacy    - the committed pooled move, K = 4 (current default)
#   am        - regularized adaptive Metropolis, K = 1 (candidate default)
# on 5 models (3 simulated with 80 subjects, 2 on the real forstmann data)
# x 2 hierarchies x 3 seeds = 90 fits.
#
# Every arm uses package defaults for preburn/burn/adapt (burn stops on
# Rhat, adapt stops on min_unique and is extended by fit until AM is frozen),
# then a FIXED-length production stage. So two things are measured:
#   (1) cost to reach production (burn + adapt time and iterations), and
#   (2) production efficiency (ESS per second on group AND subject parameters),
# plus (3) correctness: all arms target the same posterior, so pooled means and
#     SDs must agree across arms within Monte Carlo error.
#
# Usage (repository root):
#   Rscript benchmarks/group-move/real-fits.R smoke          # mechanics only
#   Rscript benchmarks/group-move/real-fits.R run [filters]  # fit, cached per case
#   Rscript benchmarks/group-move/real-fits.R summarize      # tables + verdict
# Filters (optional, for sharding across processes): model=RDM6,LBA9
#   type=standard seed=1,2 arm=am; cores=8 sets cores_per_chain (default 8,
#   so one process uses 3 x 8 = 24 cores -- run shards one after another,
#   not concurrently, unless cores is lowered so the total fits the machine).
# Arms of one (model, type, seed) cell always run back-to-back in the same
# process so they see the same machine load. Do not run shards while other
# heavy jobs are on the machine; timing is the outcome.

suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
if (!requireNamespace("posterior", quietly = TRUE))
  stop("Install posterior for benchmark diagnostics")

args <- commandArgs(TRUE)
mode <- if (length(args)) args[[1L]] else "run"
stopifnot(mode %in% c("smoke", "run", "summarize"))
filters <- list()
for (a in args[-1L]) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  filters[[kv[1L]]] <- strsplit(kv[2L], ",", fixed = TRUE)[[1L]]
}
smoke <- mode == "smoke"
cores_per_chain <- as.integer(if (is.null(filters$cores)) 8L else filters$cores)
filters$cores <- NULL

outdir <- file.path("benchmarks/group-move/results",
                    if (smoke) "real-fits-smoke" else "real-fits")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- settings ---------------------------------------------------------------
arms <- c("original", "legacy", "am")
types <- c("standard", "diagonal-gamma")
seeds <- if (smoke) 1L else 1:3
n_subjects <- if (smoke) 6L else 80L
n_trials <- if (smoke) 10L else 50L         # per S x E cell: 300 per subject
n_chains <- 3L
production <- if (smoke) 100L else 2000L
preburn <- if (smoke) 25L else 100L
particle_factor <- 50L                      # package default
max_tries <- if (smoke) 2L else 20L
smoke_am <- list(emc2.group_move_warmup = 60L, emc2.group_move_settle = 20L)

# Any change to package code, this script, or run settings invalidates the
# cached fits. In particular, worker count changes the timing being measured.
run_settings <- c(smoke = as.integer(smoke), cores_per_chain = cores_per_chain,
  n_subjects = n_subjects, n_trials = n_trials, n_chains = n_chains,
  production = production, preburn = preburn, particle_factor = particle_factor,
  max_tries = max_tries)
fingerprint <- paste(
  paste(tools::md5sum(c(list.files("R", full.names = TRUE),
    list.files("src", pattern = "\\.(cpp|h)$", full.names = TRUE),
    "benchmarks/group-move/real-fits.R")), collapse = "|"),
  paste(names(run_settings), run_settings, sep = "=", collapse = "|"),
  sep = "::")

# ---- models and simulated data ---------------------------------------------
# Simulated from a known hierarchy, so there are enough subjects (forstmann has
# only 19) and the truth is available. Group means are on EMC2's sampled scale;
# the true group covariance has the package's default per-parameter spread
# (make_random_effects, variance_proportion = 0.2) and a 0.3 equicorrelation,
# so the standard hierarchy has real correlations to find. One dataset per
# model, shared by every seed and arm.
matchfun <- function(d) d$S == d$lR
ADmat <- matrix(c(-0.5, 0.5), ncol = 1L, dimnames = list(NULL, "d"))
factors <- list(subjects = seq_len(n_subjects), S = c("left", "right"),
                E = c("speed", "neutral", "accuracy"))

# Group dimensions (mu only; the move adds one scale coordinate per parameter):
#   RDM6 = 6, LBA9 = 9, RDM14 = 14. RDM14 is the realistic wide case.
models <- list(
  RDM6 = list(model = RDM,
    formula = list(v ~ lM, B ~ lR, A ~ 1, t0 ~ 1, s ~ 1),
    contrasts = list(v = list(lM = ADmat), B = list(lR = ADmat)),
    constants = c(s = log(1)),
    truth = c(v = log(2.5), v_lMd = 1, B = log(1.5), B_lRd = 0,
              A = log(0.3), t0 = log(0.2))),
  LBA9 = list(model = LBA,
    formula = list(v ~ lM, sv ~ lM, B ~ E + lR, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat), sv = list(lM = ADmat),
                     B = list(lR = ADmat)),
    constants = c(sv = log(1)),
    truth = c(v = 2.5, v_lMd = 2, sv_lMd = -0.2, B = log(1),
              B_Eneutral = 0.2, B_Eaccuracy = 0.4, B_lRd = 0,
              A = log(0.5), t0 = log(0.2))),
  RDM14 = list(model = RDM,
    formula = list(v ~ E * lM, B ~ E * lR, A ~ 1, t0 ~ 1, s ~ 1),
    contrasts = list(v = list(lM = ADmat), B = list(lR = ADmat)),
    constants = c(s = log(1)),
    truth = c(v = log(2.5), v_Eneutral = 0, v_Eaccuracy = -0.1, v_lMd = 1,
              "v_Eneutral:lMd" = 0.1, "v_Eaccuracy:lMd" = 0.2,
              B = log(1.2), B_Eneutral = 0.2, B_Eaccuracy = 0.4, B_lRd = 0,
              "B_Eneutral:lRd" = 0, "B_Eaccuracy:lRd" = 0,
              A = log(0.3), t0 = log(0.2)))
)
# Real data (forstmann: 19 subjects, all ~830 trials each). No known truth, but
# the correctness gate is cross-arm agreement, which does not need one.
models$RDM6_forstmann <- modifyList(models$RDM6, list(truth = NULL, source = "forstmann"))
models$RDM14_forstmann <- modifyList(models$RDM14, list(truth = NULL, source = "forstmann"))
if (smoke) models <- models[c("RDM6", "RDM6_forstmann")]

make_design <- function(cfg, data = NULL) suppressMessages(
  if (is.null(data)) design(factors = factors, Rlevels = c("left", "right"),
    model = cfg$model, matchfun = matchfun, formula = cfg$formula,
    contrasts = cfg$contrasts, constants = cfg$constants, report_p_vector = FALSE)
  else design(data = data, model = cfg$model, matchfun = matchfun,
    formula = cfg$formula, contrasts = cfg$contrasts, constants = cfg$constants,
    report_p_vector = FALSE))

simulate_model <- function(model_name) {
  cache <- file.path(outdir, paste0("data-", model_name, ".rds"))
  if (file.exists(cache)) {
    x <- readRDS(cache)
    if (identical(x$fingerprint, fingerprint)) return(x)
  }
  cfg <- models[[model_name]]
  if (identical(cfg$source, "forstmann")) {
    data("forstmann", package = "EMC2", envir = environment())
    keep_s <- unique(forstmann$subjects)[seq_len(if (smoke) 4L else 19L)]
    dat <- droplevels(forstmann[forstmann$subjects %in% keep_s, ])
    if (smoke) dat <- do.call(rbind, lapply(split(dat, dat$subjects), head, 60L))
    return(list(data = dat, fingerprint = fingerprint))
  }
  des <- make_design(cfg)
  mu <- cfg$truth[names(sampled_pars(des))]
  stopifnot(!anyNA(mu))
  # make_random_effects' default spread rule (variance_proportion = 0.2)
  sds <- vapply(names(mu), function(p) {
    f <- des$model()$transform$func[[get_p_types(p)]]
    if (identical(f, "exp")) sqrt(log1p(0.2^2)) else 0.2 * max(abs(mu[[p]]), 1)
  }, numeric(1))
  R <- matrix(0.3, length(mu), length(mu)); diag(R) <- 1
  Sigma <- diag(sds) %*% R %*% diag(sds)
  dimnames(Sigma) <- list(names(mu), names(mu))
  set.seed(20260926L + match(model_name, names(models)))
  subj <- make_random_effects(des, mu, n_subj = n_subjects, covariances = Sigma)
  dat <- make_data(subj, des, n_trials = n_trials)
  list(data = dat, truth_mu = mu, truth_Sigma = Sigma, subject_pars = subj,
       fingerprint = fingerprint)
}

keep <- function(name, values) {
  if (is.null(filters[[name]])) values else values[as.character(values) %in% filters[[name]]]
}
expected_arms <- keep("arm", arms)
expected_seeds <- keep("seed", seeds)

set_arm <- function(arm) {
  options(emc2.group_move = arm != "original",
          emc2.group_move_method = if (arm == "legacy") "legacy" else "am",
          emc2.group_move_proposals = if (arm == "legacy") 4L else 1L,
          emc2.group_move_scale = TRUE,
          emc2.group_move_warmup = NULL, emc2.group_move_settle = NULL)
  if (smoke && arm == "am") options(smoke_am)
}

# ---- one fit ----------------------------------------------------------------
run_case <- function(model_name, type, seed, arm) {
  cfg <- models[[model_name]]
  dat <- simulate_model(model_name)$data
  des <- make_design(cfg, dat)
  set_arm(arm)
  # Same seed across arms: identical start points and data, only the kernel differs.
  set.seed(1000L * seed + 7L)
  emc <- suppressMessages(make_emc(dat, des, type = type, n_chains = n_chains,
                                   compress = TRUE, rt_resolution = NULL))
  stage_crit <- list(
    preburn = list(iter = preburn),
    burn = NULL,                       # package default: mean_gd 1.1
    adapt = if (smoke) list(min_unique = 50L) else NULL,  # default min_unique
    sample = list(iter = production)   # fixed length, no Rhat stopping
  )
  timing <- c(preburn = NA, burn = NA, adapt = NA, sample = NA)
  cpu <- timing
  for (stage in names(stage_crit)) {
    crit <- EMC2:::get_stop_criteria(stage, stage_crit[[stage]], type)
    if (stage == "sample") { crit$max_gd <- NULL; crit$mean_gd <- NULL }
    t0 <- proc.time()
    emc <- suppressWarnings(run_emc(emc, stage = stage, stop_criteria = crit,
      cores_for_chains = n_chains, cores_per_chain = cores_per_chain,
      particle_factor = particle_factor, max_tries = max_tries,
      step_size = 100L, verbose = FALSE, trim = FALSE))
    dt <- proc.time() - t0
    timing[stage] <- dt[["elapsed"]]
    cpu[stage] <- sum(dt[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
  }
  iters <- chain_n(emc)[1L, ]
  list(draws = extract_draws(emc, type), timing = timing, cpu = cpu,
       iters = iters, gm = group_move_diagnostics(emc),
       dimension = nrow(emc[[1L]]$samples$theta_mu))
}

# iterations x chains x variables arrays for the production stage
extract_draws <- function(emc, type) {
  idx <- which(emc[[1L]]$samples$stage == "sample")
  pn <- rownames(emc[[1L]]$samples$theta_mu)
  grab <- function(f, names) {
    out <- array(NA_real_, c(length(idx), length(emc), length(names)),
                 dimnames = list(NULL, NULL, names))
    for (ch in seq_along(emc)) out[, ch, ] <- f(emc[[ch]]$samples)
    out
  }
  mu <- grab(function(s) t(s$theta_mu[, idx, drop = FALSE]), paste0("mu_", pn))
  lsd <- grab(function(s) t(vapply(idx, function(i)
    0.5 * log(diag(s$theta_var[, , i])), numeric(length(pn)))), paste0("logsd_", pn))
  out <- list(mu = mu, logsd = lsd)
  if (type == "standard" && length(pn) > 1L) {
    pairs <- which(upper.tri(diag(length(pn))), arr.ind = TRUE)
    cn <- paste0("cor_", pn[pairs[, 1L]], "__", pn[pairs[, 2L]])
    out$cor <- grab(function(s) t(vapply(idx, function(i)
      stats::cov2cor(s$theta_var[, , i])[pairs], numeric(nrow(pairs)))), cn)
  }
  an <- as.vector(outer(pn, dimnames(emc[[1L]]$samples$alpha)[[2L]], paste, sep = "|"))
  out$alpha <- grab(function(s) t(matrix(s$alpha[, , idx, drop = FALSE],
                                         ncol = length(idx))), paste0("alpha_", an))
  out
}

block_metrics <- function(z) {
  data.frame(variable = dimnames(z)[[3L]],
    rhat = apply(z, 3L, posterior::rhat),
    ess_bulk = apply(z, 3L, posterior::ess_bulk),
    ess_tail = apply(z, 3L, posterior::ess_tail),
    mean = apply(z, 3L, mean), sd = apply(z, 3L, stats::sd),
    mcse_mean = apply(z, 3L, posterior::mcse_mean),
    stringsAsFactors = FALSE)
}

case_rows <- function(model_name, type, seed, arm, res) {
  do.call(rbind, lapply(names(res$draws), function(b) {
    m <- block_metrics(res$draws[[b]])
    cbind(data.frame(model = model_name, type = type, seed = seed, arm = arm,
                     block = b, stringsAsFactors = FALSE), m,
          dimension = res$dimension,
          burn_iter = res$iters[["burn"]], adapt_iter = res$iters[["adapt"]],
          t_preburn = res$timing[["preburn"]], t_burn = res$timing[["burn"]],
          t_adapt = res$timing[["adapt"]], t_sample = res$timing[["sample"]],
          cpu_sample = res$cpu[["sample"]],
          gm_sample_acceptance = mean(res$gm$sample_acceptance, na.rm = TRUE))
  }))
}

# ---- run --------------------------------------------------------------------
if (mode %in% c("smoke", "run")) {
  for (model_name in keep("model", names(models)))
    for (type in keep("type", types))
      for (seed in keep("seed", seeds))
        for (arm in keep("arm", arms)) {
          key <- sprintf("%s-%s-seed%d-%s", model_name, type, seed, arm)
          cache <- file.path(outdir, paste0(key, ".rds"))
          old <- if (file.exists(cache)) readRDS(cache) else NULL
          if (!is.null(old) && identical(old$fingerprint, fingerprint)) next
          cat(format(Sys.time()), "start", key, "\n")
          res <- tryCatch(run_case(model_name, type, seed, arm),
                          error = function(e) list(error = conditionMessage(e)))
          res$fingerprint <- fingerprint
          if (is.null(res$error)) res$metrics <- case_rows(model_name, type, seed, arm, res)
          saveRDS(res, cache)
          if (!is.null(res$error)) { cat("  ERROR:", res$error, "\n"); next }
          g <- res$metrics[res$metrics$block %in% c("mu", "logsd"), ]
          cat(sprintf("  burn %d adapt %d | t burn %.0fs adapt %.0fs sample %.0fs | group min ESS %.0f, max Rhat %.3f\n",
                      res$iters[["burn"]], res$iters[["adapt"]], res$timing[["burn"]],
                      res$timing[["adapt"]], res$timing[["sample"]],
                      min(g$ess_bulk), max(g$rhat)))
        }
  if (smoke) mode <- "summarize"
}

# ---- summarize --------------------------------------------------------------
if (mode == "summarize") {
  files <- list.files(outdir, pattern = "\\.rds$", full.names = TRUE)
  res <- lapply(files, readRDS)
  errs <- vapply(res, function(r) !is.null(r$error), logical(1))
  if (any(errs)) {
    cat("FAILED CASES:\n")
    for (i in which(errs)) cat(" ", basename(files[i]), ":", res[[i]]$error, "\n")
  }
  stale <- vapply(res, function(r) !identical(r$fingerprint, fingerprint), logical(1))
  if (any(stale)) cat(sum(stale), "cached results are STALE (code changed) and excluded\n")
  m <- do.call(rbind, lapply(res[!errs & !stale], `[[`, "metrics"))
  if (is.null(m)) stop("No results")
  m$group <- m$block %in% c("mu", "logsd", "cor")

  # (1)+(2) per fit: cost to production and production efficiency
  per_fit <- do.call(rbind, lapply(split(m, list(m$model, m$type, m$seed, m$arm), drop = TRUE),
    function(x) {
      g <- x[x$group, ]; a <- x[x$block == "alpha", ]
      data.frame(model = x$model[1], type = x$type[1], seed = x$seed[1], arm = x$arm[1],
        burn_iter = x$burn_iter[1], adapt_iter = x$adapt_iter[1],
        t_to_sample = x$t_preburn[1] + x$t_burn[1] + x$t_adapt[1],
        t_sample = x$t_sample[1],
        group_max_rhat = max(g$rhat), group_min_ess = min(g$ess_bulk),
        group_min_tail_ess = min(g$ess_tail), group_median_ess = stats::median(g$ess_bulk),
        alpha_min_ess = min(a$ess_bulk), alpha_median_ess = stats::median(a$ess_bulk),
        group_min_ess_s = min(g$ess_bulk) / x$t_sample[1],
        alpha_median_ess_s = stats::median(a$ess_bulk) / x$t_sample[1],
        stringsAsFactors = FALSE)
    }))
  per_fit <- per_fit[order(per_fit$model, per_fit$type, per_fit$seed, per_fit$arm), ]
  write.csv(per_fit, file.path(outdir, "per-fit.csv"), row.names = FALSE)

  # Efficiency ratios vs original, paired within (model, type, seed)
  ratio <- do.call(rbind, lapply(split(per_fit, list(per_fit$model, per_fit$type, per_fit$seed), drop = TRUE),
    function(x) {
      o <- x[x$arm == "original", ]
      if (!nrow(o)) return(NULL)
      y <- x[x$arm != "original", ]
      data.frame(model = y$model, type = y$type, seed = y$seed, arm = y$arm,
        group_min_ess_s = y$group_min_ess_s / o$group_min_ess_s,
        alpha_median_ess_s = y$alpha_median_ess_s / o$alpha_median_ess_s,
        t_to_sample = y$t_to_sample / o$t_to_sample, stringsAsFactors = FALSE)
    }))
  write.csv(ratio, file.path(outdir, "ratios-vs-original.csv"), row.names = FALSE)

  # (3) correctness: pool seeds per arm (each seed's chains are more chains of
  # the same target) and compare each arm with original.
  pooled <- do.call(rbind, lapply(split(m, list(m$model, m$type, m$arm, m$variable), drop = TRUE),
    function(x) data.frame(model = x$model[1], type = x$type[1], arm = x$arm[1],
      block = x$block[1], variable = x$variable[1], mean = mean(x$mean),
      sd = sqrt(mean(x$sd^2)), mcse = sqrt(sum(x$mcse_mean^2)) / nrow(x),
      seed_spread = if (nrow(x) > 1) stats::sd(x$mean) / sqrt(nrow(x)) else NA_real_,
      max_rhat = max(x$rhat), ess = sum(x$ess_bulk), stringsAsFactors = FALSE)))
  # Error scale: the larger of the within-run MCSE and the between-seed spread,
  # so an arm that is stuck in different places per seed cannot pass by luck.
  pooled$err <- pmax(pooled$mcse, pooled$seed_spread, na.rm = TRUE)
  ref <- pooled[pooled$arm == "original", c("model", "type", "variable", "mean", "sd", "err", "max_rhat", "ess")]
  cmp <- merge(pooled[pooled$arm != "original", ], ref,
               by = c("model", "type", "variable"), suffixes = c("", "_orig"))
  cmp$z <- (cmp$mean - cmp$mean_orig) / sqrt(cmp$err^2 + cmp$err_orig^2)
  cmp$sd_ratio <- cmp$sd / cmp$sd_orig
  write.csv(cmp, file.path(outdir, "agreement-vs-original.csv"), row.names = FALSE)

  cat("\n== Per fit ==\n"); print(per_fit, row.names = FALSE, digits = 3)
  cat("\n== Efficiency relative to original (paired by seed; >1 = better) ==\n")
  agg <- aggregate(cbind(group_min_ess_s, alpha_median_ess_s, t_to_sample) ~ model + type + arm,
                   ratio, function(v) exp(mean(log(v))))
  print(agg, row.names = FALSE, digits = 3)
  cat("\n== Agreement with original (group + alpha; well-mixed variables only) ==\n")
  ok <- cmp$ess >= 400 & cmp$ess_orig >= 400
  agr <- do.call(rbind, lapply(split(cmp[ok, ], list(cmp$model[ok], cmp$type[ok], cmp$arm[ok]), drop = TRUE),
    function(x) data.frame(model = x$model[1], type = x$type[1], arm = x$arm[1],
      n = nrow(x), max_abs_z = max(abs(x$z)), n_abs_z_gt4 = sum(abs(x$z) > 4),
      sd_ratio_range = paste(round(range(x$sd_ratio), 3), collapse = "-"),
      stringsAsFactors = FALSE)))
  print(agr, row.names = FALSE, digits = 3)

  # ---- verdict for AM becoming the default ----
  am_r <- ratio[ratio$arm == "am", ]
  am_c <- cmp[cmp$arm == "am" & ok, ]
  cells <- unique(per_fit[, c("model", "type")])
  complete <- nrow(per_fit) == nrow(cells) * length(expected_arms) *
    length(expected_seeds) && !any(errs)
  orig_converged <- all(per_fit$group_max_rhat[per_fit$arm == "original"] < 1.05)
  gates <- c(
    "all cases ran, none errored" = complete,
    "AM production chains converged (group Rhat < 1.05, every fit)" =
      all(per_fit$group_max_rhat[per_fit$arm == "am"] < 1.05),
    "AM agrees with original (no |z| > 4; sd ratio in [0.9, 1.1])" =
      nrow(am_c) > 0 && all(abs(am_c$z) <= 4) && all(abs(am_c$sd_ratio - 1) <= 0.1),
    "AM group min ESS/s >= original (geo-mean over cells)" =
      exp(mean(log(am_r$group_min_ess_s))) >= 1,
    "AM group min ESS/s >= 0.8x original in EVERY fit" =
      all(am_r$group_min_ess_s >= 0.8),
    "AM subject (alpha) median ESS/s >= 0.9x original (geo-mean)" =
      exp(mean(log(am_r$alpha_median_ess_s))) >= 0.9,
    "AM time to production <= 1.25x original (geo-mean)" =
      exp(mean(log(am_r$t_to_sample))) <= 1.25)
  cat("\n== Default-change gates (AM vs original) ==\n")
  for (g in names(gates)) cat(if (isTRUE(gates[[g]])) " PASS  " else " FAIL  ", g, "\n")
  if (!orig_converged)
    cat(" NOTE  original did not converge in every fit: its agreement reference is weak there;",
        "report which cells, do not treat as an AM failure.\n")
  writeLines(capture.output(sessionInfo()), file.path(outdir, "session.txt"))
}
