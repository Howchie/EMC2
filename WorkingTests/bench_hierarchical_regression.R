# The audit's regression benchmark: matched hierarchical fits, reported as
# posterior per second rather than as iterations per second.
#
#   Rscript WorkingTests/bench_hierarchical_regression.R
#   EMC_CASES=cores Rscript WorkingTests/bench_hierarchical_regression.R
#   EMC_CASES=size EMC_OUT=after.rds EMC_BASELINE=before.rds Rscript ...
#
# Every later commit in the audit sequence that claims a speed-up is judged
# here.  Three numbers decide it, and they are deliberately not "time per
# iteration":
#
#   * worst-relevant-parameter bulk and tail ESS per elapsed second.  A change
#     that makes an iteration cheaper and the chain stickier is not an
#     improvement, and only an effective-sample measure can say so.  The
#     divisor is the *whole* wall clock -- preparation, start-up, every stage,
#     proposal updates, diagnostics, checkpoint writes and teardown -- because
#     that is what a user waits for.  The sample-stage-only rate is printed
#     next to it so a change in where the time goes is visible.  Worker pool
#     start-up and teardown happen inside a stage and are charged to it; the
#     `cleanup` column is only the checkpoint file's removal.
#   * peak private memory of the whole process tree, PSS where Linux offers it.
#     Summing RSS over a forked tree counts every shared page once per process,
#     so a pool of eight workers appears to hold eight times the data it has.
#   * full wall clock, split by stage.
#
# Environment:
#
#   EMC_CASES       case or group names, comma separated, or "all" (smoke)
#   EMC_MODEL       RDM (default), LBA, BAwL or DDM
#   EMC_LIB         library to load EMC2 from
#   EMC_ITER_SCALE  multiplies every stage's iteration count (1)
#   EMC_STEP        run_emc block size (100, capped at the shortest stage)
#   EMC_REPS        repeats of each case; >1 also checks determinism (1)
#   EMC_OUT         RDS to write   EMC_BASELINE  RDS to compare against
#   EMC_LABEL       name for this run in the output (the library's basename)
#   EMC_PROFILE     also print the preburn stage's per-iteration profile
#   EMC_ESS_ALPHA   include the N x P subject parameters in "worst"
#   EMC_SEED        (20260910)     EMC_POLL  memory poll interval, s (0.25)
#
# How the shapes are chosen.  `cond` is a synthetic factor whose level count
# sets the parameter vector's width, so P is exact rather than whatever a
# chosen interaction happened to produce: P = levels(cond) + 6 for RDM, LBA and
# DDM, + 7 for BAwL, and a case asking for less than a model's floor is raised
# to it and renamed.  Trial counts are uneven across subjects, and they scale
# with the level count so that every subject observes every condition.  P and
# dataset size are therefore *not* orthogonal across the size grid: comparing
# n32_p64 against n32_p250 measures both at once.  Only matched comparisons of
# the same case -- before against after, or one library against another -- mean
# anything, which is what this script is for.
#
# Traps, each of which has cost real time here:
#
#  * `nproc` reports 1 in some harnesses purely because `OMP_NUM_THREADS=1` is
#    exported.  `parallel::detectCores()` is the number that governs the pool.
#  * Core budgets are 1/2/4/8/16, never "all detected cores".  A run that takes
#    every core measures the machine's contention with the rest of the system.
#  * Do NOT A/B an R-only change by installing two temp libraries.  Separate
#    `R CMD INSTALL` runs differ in lazy-load database and bytecode layout by
#    several percent of whole-fit CPU, which is larger than most changes worth
#    measuring.  Patch the namespace in-process instead.
#  * The fit runs in a forked child so the parent can sample the tree's memory
#    while it works.  The child builds the data as well as running the fit, so
#    that the parent holds nothing the child shares: PSS divides a shared page
#    by the number of processes mapping it, so a parent holding the same data
#    would quietly halve the answer.
#  * A parent-link walk does NOT find a fit's workers.  The default `spawn`
#    backend starts each pool template with `system2(wait = FALSE)`, so init
#    adopts it and neither it nor its workers is a descendant of the fit any
#    more.  On the eight-subject case that is nine processes measured instead
#    of four and roughly twice the memory.  `.emc_profile_tree_pids` claims
#    them back through the `master=<pid>` record the pool writes into the
#    template's command line.
#  * The per-iteration profile record (R/profile_schema.R) survives only the
#    preburn stage; `concat_emc` rebuilds the sample arrays and drops the
#    attribute.  Per-block costs are C16's work, not this script's.

lib <- Sys.getenv("EMC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))

.env_int <- function(name, default) {
  as.integer(Sys.getenv(name, as.character(default)))
}
.env_num <- function(name, default) {
  as.numeric(Sys.getenv(name, as.character(default)))
}
.env_lgl <- function(name, default) {
  identical(tolower(Sys.getenv(name, tolower(as.character(default)))), "true")
}
.now <- function() proc.time()[["elapsed"]]
`%||%` <- function(x, y) if (is.null(x)) y else x

SEED       <- .env_int("EMC_SEED", 20260910L)
MODEL      <- Sys.getenv("EMC_MODEL", "RDM")
ITER_SCALE <- .env_num("EMC_ITER_SCALE", 1)
REPS       <- .env_int("EMC_REPS", 1L)
POLL       <- .env_num("EMC_POLL", 0.25)
OUT        <- Sys.getenv("EMC_OUT", file.path(tempdir(), "bench_hierarchical.rds"))
BASELINE   <- Sys.getenv("EMC_BASELINE", "")
LABEL      <- Sys.getenv("EMC_LABEL", if (nzchar(lib)) basename(lib) else "current")
PROFILE    <- .env_lgl("EMC_PROFILE", FALSE)
ESS_ALPHA  <- .env_lgl("EMC_ESS_ALPHA", FALSE)

# Fixed iteration counts, not adaptive stopping.  `fit()`'s defaults stop burn
# on Rhat and adapt on unique-sample counts, so two libraries would do different
# amounts of work and the comparison would be of luck rather than of speed.
# `min_unique = 2` is what pins adapt; without it `get_stop_criteria` supplies
# 600 and the stage runs until it is met.
#
# `setNames` is not decoration: `pmax` drops names, and a nameless ITERS would
# leave the stage loop with nothing to iterate over and report a fit that never
# ran as a very fast one.
ITERS <- setNames(pmax(1L, round(c(100, 100, 100, 300) * ITER_SCALE)),
                  c("preburn", "burn", "adapt", "sample"))
# One block per stage where possible: `run_emc` runs in chunks of `step_size`
# and rounds a stage up to a whole chunk, so a step larger than a stage's
# iteration count would overshoot it.
STEP <- min(.env_int("EMC_STEP", 100L), min(ITERS))

.stop_for <- function(stage) {
  out <- list(iter = unname(ITERS[[stage]]))
  if (identical(stage, "adapt")) out$min_unique <- 2L
  out
}

# --- model specifications ---------------------------------------------------

# `offset` is the number of parameters the formulas contribute besides `cond`,
# so a case can ask for a width and get it.
.model_spec <- function(model) {
  switch(model,
    RDM  = list(fun = RDM,  offset = 6L, race = TRUE,
                formula = list(v ~ lM + cond, B ~ E, A ~ 1, t0 ~ 1),
                constants = c(s = 0)),
    LBA  = list(fun = LBA,  offset = 6L, race = TRUE,
                formula = list(v ~ lM + cond, B ~ E, A ~ 1, t0 ~ 1),
                constants = c(sv = log(1))),
    BAwL = list(fun = BAwL, offset = 7L, race = TRUE,
                formula = list(v ~ lM + cond, B ~ E, A ~ 1, t0 ~ 1, k ~ 1),
                constants = c(sv = 0)),
    DDM  = list(fun = DDM,  offset = 6L, race = FALSE,
                formula = list(v ~ S + cond, a ~ E, t0 ~ 1, Z ~ 1),
                constants = c(s = log(1))),
    stop("no benchmark specification for model ", model))
}

# --- case registry ----------------------------------------------------------

.case <- function(name = NULL, group, subjects, pars, chains = 2L,
                  cores_per_chain = 2L, particle_factor = 50, trials = NULL) {
  spec <- .model_spec(MODEL)
  # `pars` is a target.  Each model's formulas contribute a fixed number of
  # parameters before `cond` says anything, and `cond` needs at least two
  # levels to be a factor at all, so a narrow width is raised to the model's
  # floor rather than refused: BAwL simply has no eight-parameter version of
  # this design, and the case is still worth running.
  pars <- max(as.integer(pars), spec$offset + 2L)
  n_cond <- pars - spec$offset
  # Enough trials that every subject sees every condition several times, and
  # never fewer than an ordinary experiment's worth.
  if (is.null(trials)) trials <- max(200L, 3L * n_cond)
  if (is.null(name)) name <- sprintf("n%d_p%d", as.integer(subjects), pars)
  list(name = name, group = group, model = MODEL,
       subjects = as.integer(subjects), pars = as.integer(pars),
       conditions = n_cond, trials = as.integer(trials),
       chains = as.integer(chains),
       cores_per_chain = as.integer(cores_per_chain),
       particle_factor = particle_factor)
}

.registry <- local({
  out <- list(
    # A single quick case, so that running the script with no arguments is a
    # useful thing to do rather than an hour's commitment.
    .case("smoke", "smoke", subjects = 8L, pars = 8L, chains = 2L,
          cores_per_chain = 2L))
  # N and P together: the audit asks for N = 8/32/140+ against P ~ 8/64/250+.
  for (n in c(8L, 32L, 140L)) {
    for (p in c(8L, 64L, 250L)) {
      # The name follows the width the model can actually reach, so a case
      # called n8_p250 always has 250 parameters in it.
      out[[length(out) + 1L]] <- .case(NULL, "size", subjects = n, pars = p,
                                       chains = 2L, cores_per_chain = 4L)
    }
  }
  # Core budget at a fixed shape.  One chain, so `cores_per_chain` is the only
  # thing varying and the chains are not competing for the same cores.
  for (cc in c(1L, 2L, 4L, 8L, 16L)) {
    out[[length(out) + 1L]] <- .case(sprintf("cores%d", cc), "cores",
                                     subjects = 32L, pars = 64L, chains = 1L,
                                     cores_per_chain = cc)
  }
  # Particle count: C17's experiment reads exactly this table, because more
  # particles buy a better proposal and cost time linearly.
  for (pf in c(25, 50, 100)) {
    out[[length(out) + 1L]] <- .case(sprintf("particles%d", pf), "particles",
                                     subjects = 32L, pars = 64L, chains = 1L,
                                     cores_per_chain = 4L, particle_factor = pf)
  }
  # Chains, with the per-chain budget held fixed: total cores therefore rise
  # with the chain count, which is how a real fit is run.
  for (ch in c(1L, 2L, 4L)) {
    out[[length(out) + 1L]] <- .case(sprintf("chains%d", ch), "chains",
                                     subjects = 32L, pars = 64L, chains = ch,
                                     cores_per_chain = 4L)
  }
  names(out) <- vapply(out, `[[`, character(1), "name")
  out
})

.select_cases <- function(spec) {
  want <- trimws(strsplit(spec, ",", fixed = TRUE)[[1L]])
  want <- want[nzchar(want)]
  if (!length(want) || "all" %in% want) return(.registry)
  groups <- vapply(.registry, `[[`, character(1), "group")
  keep <- names(.registry) %in% want | groups %in% want
  if (!any(keep)) {
    stop("no case or group matches ", shQuote(spec), ". Known cases: ",
         paste(names(.registry), collapse = ", "), "; groups: ",
         paste(unique(groups), collapse = ", "))
  }
  .registry[keep]
}

# --- data -------------------------------------------------------------------

.bench_data <- function(cs, seed) {
  set.seed(seed)
  base <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  levs <- sprintf("k%03d", seq_len(cs$conditions))
  # Uneven trial counts, deterministic given the seed.  Real fits do not hand
  # every subject the same number of usable trials, and an even split hides the
  # partition imbalance the worker metrics exist to measure.
  n_i <- as.integer(round(cs$trials * runif(cs$subjects, 0.7, 1.3)))
  n_i <- pmax(n_i, 2L * cs$conditions)
  dat <- do.call(rbind, lapply(seq_len(cs$subjects), function(s) {
    d <- base[rep(seq_len(nrow(base)), length.out = n_i[s]), , drop = FALSE]
    d$subjects <- sprintf("s%03d", s)
    # Balanced within subject and then permuted: assigning conditions in order
    # would alias `cond` with the E and S cycle the source rows already carry,
    # and the design would collapse.
    d$cond <- factor(sample(rep_len(levs, n_i[s])), levels = levs)
    d
  }))
  dat$subjects <- factor(dat$subjects, levels = sprintf("s%03d", seq_len(cs$subjects)))
  rownames(dat) <- NULL
  attr(dat, "trials_per_subject") <- n_i
  dat
}

# --- effective sample size --------------------------------------------------

# Rank-normalised bulk and split-tail ESS, from R/diagnostics.R, over every
# parameter of every requested selection.  "Relevant" is the group level by
# default: the group mean and variance are what a hierarchical fit is for, and
# including all N x P subject parameters would let one badly identified
# subject decide the headline for every case.
.bench_ess <- function(emc, selections) {
  rows <- list()
  for (sel in selections) {
    got <- tryCatch(get_pars(emc, selection = sel, stage = "sample",
                             merge_chains = FALSE),
                    error = function(e) NULL)
    if (is.null(got)) next
    for (nm in names(got)) {
      mats <- EMC2:::prep_mcmc_diagnostics(got[[nm]])
      for (par in names(mats)) {
        m <- mats[[par]]
        rows[[length(rows) + 1L]] <- data.frame(
          selection = sel,
          # For "mu" and "sigma2" the list is keyed by selection and the
          # columns are parameters; for "alpha" it is the other way round.
          parameter = if (identical(nm, sel)) par else paste0(nm, "[", par, "]"),
          bulk = EMC2:::ess_bulk(m), tail = EMC2:::ess_tail(m),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- one case, inside the forked child --------------------------------------

.run_case <- function(cs) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  spec <- .model_spec(cs$model)

  t0 <- .now()
  dat <- .bench_data(cs, SEED)
  des_args <- list(data = dat, model = spec$fun, formula = spec$formula,
                   constants = spec$constants, report_p_vector = FALSE)
  if (spec$race) des_args$matchfun <- function(d) d$S == d$lR
  des <- suppressMessages(do.call(design, des_args))
  n_pars <- length(sampled_pars(des))
  # A width that is not what the case asked for makes every cross-case reading
  # wrong, and it is the kind of thing a design change would do silently.
  if (n_pars != cs$pars) {
    warning(sprintf("case %s asked for %d parameters; the design has %d",
                    cs$name, cs$pars, n_pars))
  }
  emc <- suppressMessages(make_emc(dat, des, n_chains = cs$chains,
                                   compress = TRUE, type = "standard"))
  prepare <- .now() - t0

  ck <- tempfile(pattern = "bench_hier_", fileext = ".RData")
  stage_wall <- setNames(numeric(length(ITERS)), names(ITERS))
  profile <- NULL
  set.seed(SEED + 1L)
  for (stage in names(ITERS)) {
    t1 <- .now()
    emc <- run_emc(emc, stage = stage, stop_criteria = .stop_for(stage),
                   step_size = STEP, cores_per_chain = cs$cores_per_chain,
                   cores_for_chains = cs$chains,
                   particle_factor = cs$particle_factor,
                   fileName = ck, verbose = FALSE, verboseProgress = FALSE)
    stage_wall[[stage]] <- .now() - t1
    # Only preburn returns the sampler's own per-iteration record: every later
    # stage goes through `concat_emc`, which rebuilds the sample arrays and
    # leaves the attribute behind.
    if (PROFILE && identical(stage, "preburn")) {
      profile <- attr(emc[[1L]]$samples, "sampler_profile")
    }
  }
  t2 <- .now()
  unlink(ck)
  finalise <- .now() - t2

  t3 <- .now()
  ess <- .bench_ess(emc, if (ESS_ALPHA) c("mu", "sigma2", "alpha")
                    else c("mu", "sigma2"))
  diagnose <- .now() - t3

  # Three counts, because they answer different questions.  `n_rows` is the
  # experiment's size; `n_kernel` is how many rows the likelihood evaluates,
  # which for a race model is one per accumulator per *unique* trial; and the
  # compression factor is trials over unique trials.  `expand` maps each trial
  # to its unique row, so its length is the trial count and its maximum is the
  # unique-trial count.
  n_kernel <- sum(vapply(emc[[1L]]$data, nrow, numeric(1)))
  compress <- vapply(emc[[1L]]$data, function(d) {
    e <- attr(d, "expand")
    if (is.null(e) || !length(e)) c(nrow(d), nrow(d)) else c(length(e), max(e))
  }, numeric(2))
  n_trials <- sum(compress[1L, ])
  n_unique <- sum(compress[2L, ])
  # Matched runs must land in the same place, not merely in the same time.  A
  # candidate that is faster and has moved the posterior has not passed.
  checksum <- sum(vapply(emc, function(x) {
    idx <- x$samples$idx
    sum(x$samples$alpha[, , idx]) + sum(x$samples$subj_ll[, idx])
  }, numeric(1)))

  wall_fit <- sum(stage_wall)
  list(case = cs, label = LABEL, model = cs$model,
       n_pars = n_pars, n_rows = nrow(dat),
       n_kernel = n_kernel, n_trials = n_trials, n_unique = n_unique,
       trials = range(attr(dat, "trials_per_subject")),
       iters = ITERS, step = STEP, seed = SEED,
       chain_n = chain_n(emc),
       prepare = prepare, stage_wall = stage_wall, finalise = finalise,
       wall_fit = wall_fit, wall_total = prepare + wall_fit + finalise,
       diagnose = diagnose, ess = ess, checksum = checksum,
       profile = profile)
}

# --- running a case with the tree watched from outside ----------------------

# `mcparallel` forks; the child writes its result to a file rather than through
# the pipe, because a 140-subject sample store is not something to serialise
# twice for no reason.  While it works, the parent reads the private memory of
# the child's whole process tree -- the chain processes it forks, their worker
# templates and the workers those fork -- which is the only place a fit's real
# memory cost appears.
.run_measured <- function(cs) {
  if (.Platform$OS.type != "unix") {
    res <- .run_case(cs)
    res$peak_bytes <- NA_real_
    res$peak_method <- NA_character_
    res$memory_samples <- 0L
    return(res)
  }
  res_file <- tempfile(pattern = "bench_hier_res_", fileext = ".rds")
  on.exit(unlink(res_file), add = TRUE)
  # The child returns TRUE, not the result: `mccollect` reports "nothing yet"
  # and "the job returned NULL" the same way, so a NULL return would make the
  # poll loop unable to tell that the fit had finished.
  job <- parallel::mcparallel({ saveRDS(.run_case(cs), res_file); TRUE },
                              detached = FALSE)
  peak <- NA_real_
  method <- NA_character_
  n_samples <- 0L
  repeat {
    m <- tryCatch(EMC2:::.emc_profile_memory(job$pid, tree = TRUE,
                                             include_self = FALSE),
                  error = function(e) NULL)
    if (!is.null(m) && is.finite(m$bytes)) {
      n_samples <- n_samples + 1L
      # Record the method of the reading that produced the peak, not of the
      # last one: a poll that catches the tree mid-teardown reads short and
      # says so, and it must not relabel a complete maximum.
      if (is.na(peak) || m$bytes > peak) {
        peak <- m$bytes
        method <- m$method
      }
    }
    got <- parallel::mccollect(job, wait = FALSE, timeout = POLL)
    if (!is.null(got)) break
  }
  if (!file.exists(res_file)) {
    stop("case ", cs$name, ": the forked fit produced no result. Child said: ",
         paste(utils::capture.output(print(got)), collapse = " | "))
  }
  res <- readRDS(res_file)
  res$peak_bytes <- peak
  res$peak_method <- method
  res$memory_samples <- n_samples
  res
}

# --- reporting --------------------------------------------------------------

.fmt_bytes <- function(x) {
  if (!length(x) || all(is.na(x))) return("n/a")
  units <- c("B", "kB", "MB", "GB", "TB")
  i <- 1L
  while (x >= 1024 && i < length(units)) { x <- x / 1024; i <- i + 1L }
  sprintf("%.2f %s", x, units[i])
}

.worst <- function(ess, what) {
  if (is.null(ess) || !nrow(ess)) return(list(value = NA_real_, parameter = NA_character_))
  i <- which.min(ess[[what]])
  list(value = ess[[what]][i], parameter = ess$parameter[i])
}

.report_case <- function(res) {
  cs <- res$case
  say <- function(...) cat(..., sep = "")
  say("\ncase ", cs$name, "  (", cs$group, ", ", res$label, ")\n")
  say(sprintf("  shape      subjects %d   P %d   conditions %d   trials/subject %d-%d\n",
              cs$subjects, res$n_pars, cs$conditions, res$trials[1L], res$trials[2L]))
  say(sprintf("             %.0f trials (%.0f unique, compression %.2fx)   %.0f likelihood rows\n",
              res$n_trials, res$n_unique,
              res$n_trials / res$n_unique, res$n_kernel))
  say(sprintf("             chains %d   cores/chain %d   particle_factor %g\n",
              cs$chains, cs$cores_per_chain, cs$particle_factor))
  say(sprintf("  wall       prepare %.2f s", res$prepare))
  for (stage in names(res$stage_wall)) {
    say(sprintf("   %s %.2f s", stage, res$stage_wall[[stage]]))
  }
  say(sprintf("   cleanup %.2f s\n", res$finalise))
  say(sprintf("             total %.2f s   (diagnostics afterwards %.2f s, not counted)\n",
              res$wall_total, res$diagnose))
  wb <- .worst(res$ess, "bulk")
  wt <- .worst(res$ess, "tail")
  sample_wall <- res$stage_wall[["sample"]]
  say(sprintf("  ESS        worst bulk %8.1f  (%s)      worst tail %8.1f  (%s)\n",
              wb$value, wb$parameter, wt$value, wt$parameter))
  say(sprintf("             per total second   bulk %7.3f /s   tail %7.3f /s\n",
              wb$value / res$wall_total, wt$value / res$wall_total))
  say(sprintf("             per sample second  bulk %7.3f /s   tail %7.3f /s\n",
              wb$value / sample_wall, wt$value / sample_wall))
  say(sprintf("  memory     peak %s  (%s, %d samples of the process tree)\n",
              .fmt_bytes(res$peak_bytes),
              if (is.na(res$peak_method)) "unavailable" else res$peak_method,
              res$memory_samples))
  say(sprintf("  posterior  checksum %.17g\n", res$checksum))
  if (!is.null(res$profile)) {
    say("\n  preburn stage, per iteration:\n")
    EMC2:::.emc_profile_report(res$profile,
                               elapsed = res$stage_wall[["preburn"]])
  }
  invisible(NULL)
}

.summary_row <- function(res) {
  wb <- .worst(res$ess, "bulk")
  wt <- .worst(res$ess, "tail")
  data.frame(case = res$case$name, group = res$case$group, label = res$label,
             rep = res$rep %||% 1L,
             subjects = res$case$subjects, pars = res$n_pars,
             chains = res$case$chains, cores = res$case$cores_per_chain,
             particle_factor = res$case$particle_factor,
             wall = res$wall_total, sample_wall = unname(res$stage_wall[["sample"]]),
             bulk = wb$value, tail = wt$value,
             bulk_per_s = wb$value / res$wall_total,
             tail_per_s = wt$value / res$wall_total,
             peak_bytes = res$peak_bytes, checksum = res$checksum,
             stringsAsFactors = FALSE)
}


.print_summary <- function(rows) {
  cat("\nsummary\n")
  show <- data.frame(
    case = rows$case, rep = rows$rep, N = rows$subjects, P = rows$pars,
    chains = rows$chains, cores = rows$cores,
    wall = sprintf("%8.2f", rows$wall),
    bulk = sprintf("%7.1f", rows$bulk), tail = sprintf("%7.1f", rows$tail),
    `bulk/s` = sprintf("%7.3f", rows$bulk_per_s),
    `tail/s` = sprintf("%7.3f", rows$tail_per_s),
    peak = vapply(rows$peak_bytes, .fmt_bytes, character(1)),
    check.names = FALSE, stringsAsFactors = FALSE)
  print(show, row.names = FALSE)
}

# A candidate is accepted on reproducible gains, so the comparison is printed
# rather than left to the reader: same case, same rep, before against after.
# A moved checksum is called out first, because a speed-up that changed the
# posterior is not a speed-up.
.compare <- function(rows, baseline_file) {
  base <- tryCatch(readRDS(baseline_file), error = function(e) NULL)
  if (is.null(base) || is.null(base$summary)) {
    cat("\nbaseline ", baseline_file, " unreadable; no comparison\n", sep = "")
    return(invisible(NULL))
  }
  b <- base$summary
  key <- function(d) paste(d$case, d$rep, sep = "#")
  common <- intersect(key(rows), key(b))
  if (!length(common)) {
    cat("\nbaseline has no case in common with this run\n")
    return(invisible(NULL))
  }
  i <- match(common, key(rows))
  j <- match(common, key(b))
  pct <- function(new, old) 100 * (new - old) / old
  out <- data.frame(
    case = rows$case[i], rep = rows$rep[i],
    wall_base = sprintf("%8.2f", b$wall[j]),
    wall_new = sprintf("%8.2f", rows$wall[i]),
    wall_pct = sprintf("%+6.1f%%", pct(rows$wall[i], b$wall[j])),
    bulk_per_s_pct = sprintf("%+6.1f%%", pct(rows$bulk_per_s[i], b$bulk_per_s[j])),
    tail_per_s_pct = sprintf("%+6.1f%%", pct(rows$tail_per_s[i], b$tail_per_s[j])),
    peak_pct = sprintf("%+6.1f%%", pct(rows$peak_bytes[i], b$peak_bytes[j])),
    posterior = ifelse(rows$checksum[i] == b$checksum[j], "same", "MOVED"),
    stringsAsFactors = FALSE)
  cat("\nagainst baseline ", basename(baseline_file), " (label ",
      b$label[j[1L]], ")\n", sep = "")
  print(out, row.names = FALSE)
  if (any(out$posterior == "MOVED")) {
    cat("\n  the posterior moved: this is a different fit, not a faster one.\n",
        "  Check the C2 equivalence gates before reading any timing above.\n",
        sep = "")
  }
  invisible(out)
}

# --- main -------------------------------------------------------------------

cases <- .select_cases(Sys.getenv("EMC_CASES", "smoke"))
cat("namespace   : ", dirname(getNamespaceInfo("EMC2", "path")), "\n", sep = "")
cat("label       : ", LABEL, "\n", sep = "")
cat("detectCores : ", parallel::detectCores(), "\n", sep = "")
cat("model       : ", MODEL, "\n", sep = "")
cat("iterations  : ",
    paste(sprintf("%s %d", names(ITERS), ITERS), collapse = "  "),
    "   step ", STEP, "\n", sep = "")
cat("cases       : ", paste(names(cases), collapse = ", "), "\n", sep = "")
cat("reps        : ", REPS, "\n", sep = "")

if (PROFILE) {
  options(emc2.sampler_profile = TRUE, emc2.sampler_profile_expensive = TRUE)
}

results <- list()
for (rep in seq_len(REPS)) {
  for (cs in cases) {
    res <- .run_measured(cs)
    res$rep <- rep
    .report_case(res)
    results[[length(results) + 1L]] <- res
  }
}

rows <- do.call(rbind, lapply(results, .summary_row))
.print_summary(rows)

# Repeated runs of the same case share a seed, so an identical checksum is a
# statement about determinism and a differing one is a bug worth chasing before
# any timing here is believed.
if (REPS > 1L) {
  by_case <- split(rows$checksum, rows$case)
  unstable <- names(by_case)[vapply(by_case, function(x) length(unique(x)) > 1L,
                                    logical(1))]
  cat("\nreproducibility: ",
      if (length(unstable)) paste("checksum differs across reps for",
                                  paste(unstable, collapse = ", "))
      else "identical posterior checksums across reps", "\n", sep = "")
}

saveRDS(list(label = LABEL, model = MODEL, seed = SEED, iters = ITERS,
             step = STEP, summary = rows, results = results,
             detected_cores = parallel::detectCores(),
             when = Sys.time()),
        OUT)
cat("\nwritten ", OUT, "\n", sep = "")

if (nzchar(BASELINE)) .compare(rows, BASELINE)
