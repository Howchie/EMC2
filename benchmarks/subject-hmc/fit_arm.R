# One default fit() arm of the optional-moves panel.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/fit_arm.R <case> <arm> <seed> <out.rds>
# case: T30 | corr15 | lbaA | forstmann
# arm: base (ensemble only) | shear | hmc | full (ensemble + shear + subject HMC)
#   | fullc (full + collapsed group move).
# Chains 3 x 5 cores. CPU is measured outside, from the run's cgroup.
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
  corr15 = local({
    fx <- readRDS("benchmarks/mode-start/results_plausible/fixture_corr15.rds")
    fx$design <- q(design(data = fx$data, model = RDM, matchfun = matchfun,
                          formula = fx$settings$formula, constants = fx$settings$constants,
                          verbose = FALSE))
    fx
  }),
  forstmann = list(data = forstmann, design = q(design(data = forstmann, model = LNR,
    matchfun = matchfun, formula = list(m ~ lM, s ~ 1, t0 ~ 1), verbose = FALSE))),
  stop("unknown case"))
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
emc <- q(make_emc(fx$data, fx$design, n_chains = 3))
t0 <- proc.time()[["elapsed"]]
emc <- suppressWarnings(q(fit(emc, cores_per_chain = 5, cores_for_chains = 3, verbose = FALSE)))
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
    t(vapply(w, function(t) f(ch, t), numeric(length(f(ch, w[1])))))
  })
  n <- min(vapply(arr, nrow, 1L))
  simplify2array(lapply(arr, function(a) a[seq_len(n), , drop = FALSE]))
}
ut <- upper.tri(diag(p))
mu <- draws(function(ch, t) ch$samples$theta_mu[, t])
lsd <- draws(function(ch, t) 0.5 * log(diag(as.matrix(ch$samples$theta_var[, , t]))))
cr <- if (emc[[1]]$type == "standard") draws(function(ch, t) cov2cor(as.matrix(ch$samples$theta_var[, , t]))[ut]) else NULL
al <- draws(function(ch, t) as.vector(ch$samples$alpha[, , t]))
diag_of <- function(a) {
  a <- aperm(a, c(1, 3, 2))   # iterations x chains x quantities
  r <- apply(a, 3, rhat); eb <- apply(a, 3, ess_bulk); et <- apply(a, 3, ess_tail)
  list(max_rhat = max(r), min_ess_bulk = min(eb), med_ess_bulk = median(eb),
       min_ess_tail = min(et), n_over_1.01 = sum(r > 1.01), rhat = r, ess_bulk = eb)
}
d <- list(mu = diag_of(mu), lsd = diag_of(lsd), alpha = diag_of(al))
if (!is.null(cr)) d$cor <- diag_of(cr)
rec <- NULL
if (!is.null(fx$group_means)) {
  z_mu <- (apply(mu, 2, mean) - fx$group_means[pn]) / apply(mu, 2, sd)
  z_lsd <- (apply(lsd, 2, mean) - log(sqrt(diag(fx$intended_covariance[pn, pn])))) / apply(lsd, 2, sd)
  rec <- list(z_mu = z_mu, z_lsd = z_lsd)
}
moves <- lapply(emc, function(ch) list(
  hmc = ch$subject_hmc$totals,
  shear = if (is.null(ch$shear)) NULL else list(iterations = ch$shear$iterations, rounds = ch$shear$rounds),
  history = ch$alloc$history,
  collapsed = if (is.null(ch$collapsed)) NULL else list(iterations = ch$collapsed$iterations,
    rounds = ch$collapsed$rounds, accepted = ch$collapsed$accepted, eps = ch$collapsed$plan$eps)))
res <- list(case = case, arm = arm, seed = seed, lib = find.package("EMC2"), elapsed = elapsed,
            iters = chain_n(emc), stage_iters = stage_iters, stage_cpu = stage_cpu, stop_ok = stop_ok, sample_try = if (is.null(ls)) NA else ls$try,
            sample_gd = if (is.null(ls)) NULL else ls$gd$stats, diag = d, recovery = rec, moves = moves,
            blocks = if (Sys.getenv("PROFILE") == "1") EMC2:::.emc_block_profile() else NULL,
            n_draws = dim(mu)[1], n_subjects = emc[[1]]$n_subjects, p = p)
saveRDS(res, out)
if (Sys.getenv("SAVE_EMC") == "1") saveRDS(emc, sub("\\.rds$", "_emc.rds", out))
cat(sprintf("%s %s seed %d: iterations run %s (kept %d) | elapsed %.0f s | fit() stop rule met: %s (sample try %s)\n", case, arm, seed,
  paste(names(stage_iters), stage_iters, collapse = " "), res$n_draws, elapsed, stop_ok, res$sample_try))
cat("  CPU by stage (s):", paste(sprintf("%s %.0f", names(stage_cpu), vapply(stage_cpu, `[[`, 0, "cpu")), collapse = ", "), "\n")
for (nm in names(d)) cat(sprintf("  %-5s max R-hat %.3f  min ESS bulk %.0f tail %.0f  median %.0f\n",
  nm, d[[nm]]$max_rhat, d[[nm]]$min_ess_bulk, d[[nm]]$min_ess_tail, d[[nm]]$med_ess_bulk))
if (!is.null(rec)) cat(sprintf("  recovery: max |z| mu %.2f, log sd %.2f\n", max(abs(rec$z_mu)), max(abs(rec$z_lsd))))
for (j in seq_along(moves)) {
  m <- moves[[j]]
  if (!is.null(m$collapsed)) cat(sprintf("  chain %d: collapsed moves with a likelihood round %s, accepted iterations %s over %s its, eps %.3f, last plan moves/it %.2f\n",
    j, format(m$collapsed$rounds), format(m$collapsed$accepted), format(m$collapsed$iterations), m$collapsed$eps,
    if (length(m$history)) m$history[[length(m$history)]]$collapsed_rounds else NA))
  h <- m$history
  last <- if (length(h)) h[[length(h)]] else NULL
  cat(sprintf("  chain %d: shear rounds %s over %s its | HMC attempts %s accepted %s cpu %.1f s calib %.1f s | last plan: rounds/it %.2f, rates>0 %d, tau %.2f -> %.2f\n",
    j, format(m$shear$rounds), format(m$shear$iterations), format(m$hmc$attempts), format(m$hmc$accepted),
    if (is.null(m$hmc$cpu)) 0 else m$hmc$cpu, if (is.null(m$hmc$calib_cpu)) 0 else m$hmc$calib_cpu,
    if (is.null(last)) 0 else sum(last$rounds), if (is.null(last)) 0L else sum(last$rate > 0),
    if (is.null(last)) NA else last$tau0, if (is.null(last)) NA else last$tau))
}
