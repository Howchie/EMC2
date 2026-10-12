# One default fit() arm of the optional-moves panel.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/fit_arm.R <case> <arm> <seed> <out.rds>
# case: T30 | corr15 | lbaA | lbaT | lbaE | forstmann | T30m | T30p (T30 with t0 drawn
#   independently of the other effects / shared by all subjects; read from
#   FIXTURE_DIR, default benchmarks/subject-hmc/fixtures)
# arm: base (ensemble only) | shear | hmc | full (ensemble + shear + subject HMC)
#   | fullc (full + collapsed group move).
# Chains 3 x CORES_PER_CHAIN cores (default 5). CPU is measured outside, from
# the run's cgroup.
# FEATURE=marginalise|pooled applies to FEATURE_PAR (default t0).
# Recovery is scored against the subjects the data were simulated from (their
# realised means and SDs), which is what the fit can see.
args <- commandArgs(TRUE)
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages({ library(EMC2, lib.loc = lib); library(posterior) })
cat("EMC2 from", find.package("EMC2"), "\n")
case <- args[1]; arm <- args[2]; seed <- as.integer(args[3]); out <- args[4]
options(emc2.subject_hmc = arm %in% c("hmc", "full", "fullc"),
        emc2.group_move = if (arm %in% c("shear", "full", "fullc")) "shear" else "ensemble",
        emc2.collapsed_move = arm == "fullc")
# PROFILE=1 records fit()'s per-block timings (prepare, proposals, sampling,
# concat, check_progress) into the result.
if (Sys.getenv("PROFILE") == "1") options(emc2.sampler_profile = TRUE)
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR
fx <- switch(case,
  T30 = readRDS("benchmarks/subject-hmc/fixtures/fixture_T30.rds"),
  lbaA = readRDS("benchmarks/subject-hmc/fixtures/fixture_lbaA.rds"),
  lbaT = readRDS("benchmarks/subject-hmc/fixtures/fixture_lbaT.rds"),
  corr15 = local({
    fx <- readRDS("benchmarks/mode-start/results_plausible/fixture_corr15.rds")
    fx$design <- q(design(data = fx$data, model = RDM, matchfun = matchfun,
                          formula = fx$settings$formula, constants = fx$settings$constants,
                          verbose = FALSE))
    fx
  }),
  forstmann = list(data = forstmann, design = q(design(data = forstmann, model = LNR,
    matchfun = matchfun, formula = list(m ~ lM, s ~ 1, t0 ~ 1), verbose = FALSE))),
  T30m = , T30p = , nback = , nback31 = , lbaE = readRDS(file.path(Sys.getenv("FIXTURE_DIR", "benchmarks/subject-hmc/fixtures"),
                                    paste0("fixture_", case, ".rds"))),
  stop("unknown case"))
feature <- Sys.getenv("FEATURE")
feature_par <- Sys.getenv("FEATURE_PAR", "t0")
prior_in <- NULL
if (nzchar(feature)) {
  model <- if (case == "forstmann") LNR else get(fx$settings$model, asNamespace("EMC2"))
  form <- if (case == "forstmann") list(m ~ lM, s ~ 1, t0 ~ 1) else fx$settings$formula
  consts <- if (case == "forstmann") NULL else fx$settings$constants
  extra <- switch(feature, marginalise = list(marginalise = feature_par),
                  pooled = list(pooled = feature_par), stop("unknown FEATURE"))
  fx$design <- q(do.call(design, c(list(data = fx$data, model = model, matchfun = matchfun,
    formula = form, constants = consts, verbose = FALSE), extra)))
}
# CPU of the whole scope (master and every worker) at each stage boundary,
# read from this process's cgroup around fit()'s own run_emc() calls.
cg_file <- file.path("/sys/fs/cgroup", sub("^0::", "", readLines("/proc/self/cgroup")[1]), "cpu.stat")
cg_cpu <- function() {
  x <- tryCatch(readLines(cg_file), error = function(e) character(0))
  u <- grep("^usage_usec", x, value = TRUE)
  if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_
}
stage_cpu <- list()
ns <- asNamespace("EMC2")
run_emc_orig <- ns$run_emc
timed_run_emc <- function(emc, stage, ...) {
  c0 <- cg_cpu(); t0 <- proc.time()[["elapsed"]]
  out <- run_emc_orig(emc, stage, ...)
  stage_cpu[[stage]] <<- c(cpu = cg_cpu() - c0, elapsed = proc.time()[["elapsed"]] - t0)
  out
}
unlockBinding("run_emc", ns); assign("run_emc", timed_run_emc, envir = ns); lockBinding("run_emc", ns)
set.seed(seed)
# A fixture may carry its own prior and sampler type (the N-back fixtures do).
if (is.null(prior_in) && !is.null(fx$prior)) prior_in <- fx$prior
emc <- if (is.null(fx$type)) q(make_emc(fx$data, fx$design, n_chains = 3, prior = prior_in)) else
  q(make_emc(fx$data, fx$design, n_chains = 3, prior = prior_in, type = fx$type))
t0 <- proc.time()[["elapsed"]]
cores_per_chain <- as.integer(Sys.getenv("CORES_PER_CHAIN", "5"))
# PARTICLE_FACTOR sets fit()'s particle_factor (default 50); VERBOSE=1 prints
# fit()'s progress; FIT_FILE=<path> autosaves the chains there as fit() runs.
particle_factor <- as.numeric(Sys.getenv("PARTICLE_FACTOR", "50"))
verbose <- identical(Sys.getenv("VERBOSE", "0"), "1")
fit_file <- Sys.getenv("FIT_FILE", "")
emc <- suppressWarnings(q(fit(emc, cores_per_chain = cores_per_chain, cores_for_chains = 3,
                              particle_factor = particle_factor, verbose = verbose,
                              fileName = if (nzchar(fit_file)) fit_file else NULL)))
elapsed <- proc.time()[["elapsed"]] - t0
emc <- EMC2:::restore_duplicates(emc)
pn <- emc[[1]]$par_names; p <- length(pn)
log <- emc[[1]]$convergence_log
last_of <- function(stage) { l <- Filter(function(e) identical(e$stage, stage), log); if (length(l)) l[[length(l)]] else NULL }
ls <- last_of("sample")
stage_iters <- vapply(c("preburn", "burn", "adapt", "sample"), function(stage) {
  l <- Filter(function(e) identical(e$stage, stage), log)
  if (length(l)) max(vapply(l, function(e) as.numeric(e$iter), numeric(1))) else 0
}, numeric(1))
stop_ok <- !is.null(ls) && !is.null(ls$gd$stats) && all(ls$gd$stats$value < ls$gd$stats$threshold)
draws <- function(f) {
  arr <- lapply(emc, function(ch) {
    w <- which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx)
    matrix(vapply(w, function(t) f(ch, t), numeric(length(f(ch, w[1])))), nrow = length(w), byrow = TRUE)
  })
  n <- min(vapply(arr, nrow, 1L))
  simplify2array(lapply(arr, function(a) a[seq_len(n), , drop = FALSE]))
}
hier <- seq_len(p)
ph <- pn[hier]
ut <- upper.tri(diag(length(hier)))
mu <- draws(function(ch, t) ch$samples$theta_mu[hier, t])
lsd <- draws(function(ch, t) 0.5 * log(diag(as.matrix(ch$samples$theta_var[, , t])))[hier])
cr <- if (emc[[1]]$type == "standard") draws(function(ch, t) cov2cor(as.matrix(ch$samples$theta_var[hier, hier, t]))[ut]) else NULL
al <- draws(function(ch, t) as.vector(ch$samples$alpha[, , t]))
pooled <- if (is.null(emc[[1]]$samples$pooled)) NULL else
  draws(function(ch, t) ch$samples$pooled[, t])
diag_of <- function(a) {
  a <- aperm(a, c(1, 3, 2))   # iterations x chains x quantities
  r <- apply(a, 3, rhat); eb <- apply(a, 3, ess_bulk); et <- apply(a, 3, ess_tail)
  list(max_rhat = max(r), min_ess_bulk = min(eb), med_ess_bulk = median(eb),
       min_ess_tail = min(et), n_over_1.01 = sum(r > 1.01), rhat = r, ess_bulk = eb)
}
d <- list(mu = diag_of(mu), lsd = diag_of(lsd), alpha = diag_of(al))
if (!is.null(cr)) d$cor <- diag_of(cr)
if (!is.null(pooled)) d$pooled <- diag_of(pooled)
rec <- NULL
if (!is.null(fx$subject_effects)) {
  se <- fx$subject_effects
  if (!is.null(rownames(se))) se <- se[dimnames(emc[[1]]$samples$alpha)[[2]], , drop = FALSE]
  zof <- function(a, truth) (apply(a, 2, mean) - truth) / apply(a, 2, sd)
  flat <- function(a) matrix(aperm(a, c(1, 3, 2)), ncol = dim(a)[2])
  z_mu <- stats::setNames(zof(flat(mu), colMeans(se[, ph, drop = FALSE])), ph)
  z_lsd <- stats::setNames(zof(flat(lsd), log(apply(se[, ph, drop = FALSE], 2, sd))), ph)
  z_alpha <- zof(flat(al), as.vector(t(se[, pn, drop = FALSE])))
  z_alpha <- matrix(z_alpha, length(pn), dimnames = list(pn, rownames(se)))
  z_pooled <- if (is.null(pooled)) NULL else {
    pnames <- rownames(emc[[1]]$samples$pooled)
    stats::setNames(zof(flat(pooled), colMeans(fx$subject_effects[, pnames, drop = FALSE])), pnames)
  }
  rec <- list(z_mu = z_mu, z_lsd = z_lsd, z_alpha = z_alpha, z_pooled = z_pooled)
}
moves <- lapply(emc, function(ch) list(
  hmc = ch$subject_hmc$totals,
  shear = if (is.null(ch$shear)) NULL else list(iterations = ch$shear$iterations, rounds = ch$shear$rounds),
  history = ch$alloc$history,
  collapsed = if (is.null(ch$collapsed)) NULL else list(iterations = ch$collapsed$iterations,
    rounds = ch$collapsed$rounds, accepted = ch$collapsed$accepted,
    last = ch$collapsed$plan$last[c("n", "acc", "moved")],
    eps = ch$collapsed$plan$eps)))
res <- list(case = case, arm = arm, seed = seed, feature = feature, feature_par = feature_par,
            lib = find.package("EMC2"), elapsed = elapsed,
            iters = chain_n(emc), stage_iters = stage_iters, stage_cpu = stage_cpu, stop_ok = stop_ok, sample_try = if (is.null(ls)) NA else ls$try,
            sample_gd = if (is.null(ls)) NULL else ls$gd$stats, diag = d, recovery = rec, moves = moves,
            blocks = if (Sys.getenv("PROFILE") == "1") EMC2:::.emc_block_profile() else NULL,
            iterations = if (Sys.getenv("PROFILE") == "1") EMC2:::.emc_iteration_profile() else NULL,
            n_draws = dim(mu)[1], n_subjects = emc[[1]]$n_subjects, p = p)
saveRDS(res, out)
if (Sys.getenv("SAVE_EMC") == "1") saveRDS(emc, sub("\\.rds$", "_emc.rds", out))
cat(sprintf("%s %s seed %d: iterations run %s (kept %d) | elapsed %.0f s | fit() stop rule met: %s (sample try %s)\n", case, arm, seed,
  paste(names(stage_iters), stage_iters, collapse = " "), res$n_draws, elapsed, stop_ok, res$sample_try))
cat("  CPU by stage (s):", paste(sprintf("%s %.0f", names(stage_cpu), vapply(stage_cpu, `[[`, 0, "cpu")), collapse = ", "), "\n")
for (nm in names(d)) cat(sprintf("  %-5s max R-hat %.3f  min ESS bulk %.0f tail %.0f  median %.0f\n",
  nm, d[[nm]]$max_rhat, d[[nm]]$min_ess_bulk, d[[nm]]$min_ess_tail, d[[nm]]$med_ess_bulk))
if (!is.null(rec)) {
  cat(sprintf("  recovery (z vs simulated subjects): mu max |z| %.2f, log sd max |z| %.2f mean z %+.2f | alpha within 2 sd %.0f%%, max |z| %.2f\n",
    max(abs(rec$z_mu)), max(abs(rec$z_lsd)), mean(rec$z_lsd), 100 * mean(abs(rec$z_alpha) < 2), max(abs(rec$z_alpha))))
  if (feature == "marginalise") cat(sprintf("  marginalised %s: group mean z %+.2f, log sd z %+.2f | by subject: within 2 sd %.0f%%, max |z| %.2f, mean z %+.2f\n", feature_par,
    rec$z_mu[[feature_par]], rec$z_lsd[[feature_par]],
    100 * mean(abs(rec$z_alpha[feature_par, ]) < 2), max(abs(rec$z_alpha[feature_par, ])), mean(rec$z_alpha[feature_par, ])))
  if (!is.null(rec$z_pooled)) cat(sprintf("  pooled %s: z %s\n", paste(names(rec$z_pooled), collapse = ","),
    paste(sprintf("%+.2f", rec$z_pooled), collapse = ",")))
}
for (j in seq_along(moves)) {
  m <- moves[[j]]
  if (!is.null(m$collapsed)) {
    l <- m$collapsed$last
    cat(sprintf("  chain %d: collapsed rounds %s over %s its, accepted iterations %s, eps %.3f | last measurement: %s rounds, mean acceptance %.3g, moved %s | last plan moves/it %.2f\n",
      j, format(m$collapsed$rounds %||% 0), format(m$collapsed$iterations %||% 0),
      format(m$collapsed$accepted %||% 0), m$collapsed$eps %||% NA_real_,
      format(l$n %||% 0), if (isTRUE(l$n > 0)) l$acc / l$n else NA_real_, format(l$moved %||% 0),
      if (length(m$history)) m$history[[length(m$history)]]$collapsed_rounds else NA))
  }
  h <- m$history
  last <- if (length(h)) h[[length(h)]] else NULL
  cat(sprintf("  chain %d: shear rounds %s over %s its | HMC attempts %s accepted %s cpu %.1f s calib %.1f s | last plan: rounds/it %.2f, rates>0 %d, tau %.2f -> %.2f\n",
    j, format(m$shear$rounds), format(m$shear$iterations), format(m$hmc$attempts), format(m$hmc$accepted),
    if (is.null(m$hmc$cpu)) 0 else m$hmc$cpu, if (is.null(m$hmc$calib_cpu)) 0 else m$hmc$calib_cpu,
    if (is.null(last)) 0 else sum(last$rounds), if (is.null(last)) 0L else sum(last$rate > 0),
    if (is.null(last)) NA else last$tau0, if (is.null(last)) NA else last$tau))
}
