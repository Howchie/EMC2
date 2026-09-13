# Matched benchmark for the large-fit sampler architecture.
#
# Run this script in a fresh R process with EMC_LIB pointing at an installed
# baseline or candidate package.  It writes one RDS result selected by EMC_OUT.
# The benchmark intentionally uses public sample arrays and the installed
# clean-template worker backend.

lib <- Sys.getenv("EMC_LIB", .libPaths()[1L])
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))

mode <- Sys.getenv("EMC_BENCH_MODE", "fit")
out_file <- Sys.getenv("EMC_OUT", file.path(tempdir(), paste0("emc_", mode, ".rds")))
label <- Sys.getenv("EMC_LABEL", basename(lib))

if (identical(mode, "history")) {
  p <- as.integer(Sys.getenv("EMC_P", "40"))
  n <- as.integer(Sys.getenv("EMC_N", "140"))
  saved <- as.integer(Sys.getenv("EMC_SAVED", "1000"))
  writes <- as.integer(Sys.getenv("EMC_WRITES", "50"))
  samples <- list(
    alpha = array(0, c(p, n, saved + writes)),
    subj_ll = matrix(0, n, saved + writes),
    idx = saved
  )
  set.seed(20260806)
  proposals <- matrix(rnorm((p + 1L) * n), p + 1L, n)
  gc()
  timing <- system.time({
    for (j in saved + seq_len(writes)) {
      samples <- EMC2:::fill_samples_RE(samples, proposals, j, p)
    }
  })
  result <- list(
    label = label, mode = mode, p = p, n = n, saved = saved,
    writes = writes, elapsed = unname(timing[["elapsed"]]),
    user = unname(timing[["user.self"]]), system = unname(timing[["sys.self"]]),
    checksum = sum(samples$alpha[, , saved + writes]) +
      sum(samples$subj_ll[, saved + writes])
  )
  saveRDS(result, out_file)
  print(result)
  quit(save = "no")
}

if (identical(mode, "concat")) {
  # How the sample store's block joins scale.
  #
  # `run_emc` concatenates the accumulated history with each new block, so total
  # copying grows quadratically in blocks. This prints what that costs against
  # what appending into owned growable capacity would cost -- i.e. what C14's
  # representation change would be worth, which is the number that decides
  # whether it is worth putting an allocated/committed distinction in front of
  # every place that reads these arrays.
  P <- as.integer(Sys.getenv("EMC_P", "8"))
  N <- as.integer(Sys.getenv("EMC_N", "140"))
  step <- as.integer(Sys.getenv("EMC_STEP", "100"))
  set.seed(20260911)
  block <- function(n) list(
    alpha = array(rnorm(P * N * n), c(P, N, n)),
    subj_ll = matrix(rnorm(N * n), N, n),
    theta_mu = matrix(rnorm(P * n), P, n),
    theta_var = array(rnorm(P * P * n), c(P, P, n)))
  cat(sprintf("\nP=%d  subjects=%d  block=%d iterations\n", P, N, step))
  cat(sprintf("%8s %12s %12s %10s\n", "blocks", "joined_s", "one_copy_s", "overhead"))
  for (B in c(5L, 10L, 20L, 30L)) {
    acc <- block(step)
    joined <- system.time({
      for (i in seq_len(B - 1L)) {
        objs <- list(acc, block(step))
        keys <- unique(unlist(lapply(objs, names)))
        acc <- do.call(mapply, c(EMC2:::.emc_join_iteration,
                                 lapply(objs, "[", keys)))
      }
    })[["elapsed"]]
    # The floor: every block copied exactly once, which is what appending into
    # capacity would cost.
    one <- system.time({
      acc2 <- block(step)
      for (i in seq_len(B - 1L)) {
        nxt <- block(step)
        acc2$alpha <- array(c(acc2$alpha, nxt$alpha),
                            c(P, N, dim(acc2$alpha)[3] + step))
      }
    })[["elapsed"]]
    cat(sprintf("%8d %12.2f %12.2f %9.1fx\n", B, joined, one, joined / one))
  }
  quit(save = "no")
}

if (!identical(mode, "fit")) stop("unknown EMC_BENCH_MODE: ", mode)

n_subjects <- as.integer(Sys.getenv("EMC_N", "140"))
n_trials <- as.integer(Sys.getenv("EMC_TRIALS", "80"))
iterations <- as.integer(Sys.getenv("EMC_ITER", "20"))
workers <- as.integer(Sys.getenv("EMC_WORKERS", "4"))
particles <- as.integer(Sys.getenv("EMC_PARTICLES", "25"))
recycle_arg <- Sys.getenv("EMC_RECYCLE", "default")
recycle <- if (identical(recycle_arg, "default")) NULL else as.integer(recycle_arg)
profile_enabled <- identical(tolower(Sys.getenv("EMC_PROFILE", "false")), "true")
p_levels <- as.integer(Sys.getenv("EMC_P_LEVELS", "0"))

base <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
base <- base[rep(seq_len(nrow(base)), length.out = n_trials), , drop = FALSE]
dat <- do.call(rbind, lapply(seq_len(n_subjects), function(s) {
  d <- base
  d$subjects <- sprintf("bench%03d", s)
  d
}))
dat$subjects <- factor(dat$subjects, levels = sprintf("bench%03d", seq_len(n_subjects)))
if (p_levels > 0L) {
  dat$bench_condition <- factor(rep(seq_len(p_levels), length.out = nrow(dat)))
}
rownames(dat) <- NULL

matchfun <- function(d) d$S == d$lR
bench_formula <- if (p_levels > 0L) {
  list(v ~ 0 + bench_condition, B ~ 1, A ~ 1, t0 ~ 1)
} else {
  list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1)
}
des <- design(
  data = dat, model = LBA, matchfun = matchfun,
  formula = bench_formula,
  constants = c(sv = log(1))
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20260806)
sampler <- suppressMessages(make_emc(
  dat, des, n_chains = 1, compress = TRUE, rt_resolution = 0.02
))[[1L]]

options(emc2.sampler_profile = profile_enabled,
        emc2.sampler_profile_expensive = profile_enabled,
        emc2.worker_backend = "spawn")
if (!is.null(recycle)) options(emc2.worker_recycle = recycle)
pool_state <- get(".emc_pool_state", envir = asNamespace("EMC2"))
pool_state$last_error <- NULL
pool_state$warned <- FALSE
timing <- system.time({
  # `search_width = 1` is what fit() and run_emc() pass.  Without it, libraries
  # older than the check_tune_settings() default adapt every epsilon to empty at
  # iteration 26, and everything after that times caught errors, not sampling.
  sampler <- suppressMessages(EMC2:::run_stages(
    sampler, stage = "preburn", iter = iterations, search_width = 1,
    particle_factor = particles / sqrt(sampler$n_pars),
    n_cores = workers, verbose = FALSE, verboseProgress = FALSE, r_cores = 1
  ))
})
profile <- attr(sampler$samples, "sampler_profile")
result <- list(
  label = label, mode = mode, n = n_subjects, p = sampler$n_pars,
  trials = n_trials,
  iterations = iterations, workers = workers, particles = particles,
  recycle = if (is.null(recycle)) "default" else recycle,
  backend = EMC2:::.emc_wpool_backend(list(model = sampler$model)),
  pool_error = pool_state$last_error,
  elapsed = unname(timing[["elapsed"]]), user = unname(timing[["user.self"]]),
  system = unname(timing[["sys.self"]]), profile = profile,
  alpha_checksum = sum(sampler$samples$alpha[, , sampler$samples$idx]),
  ll_checksum = sum(sampler$samples$subj_ll[, sampler$samples$idx]),
  final_alpha = sampler$samples$alpha[, , sampler$samples$idx],
  final_ll = sampler$samples$subj_ll[, sampler$samples$idx]
)
saveRDS(result, out_file)
print(result[setdiff(names(result), c("profile", "final_alpha", "final_ll"))])
# Same schema and same formatter as WorkingTests/bench_worker_pool.R, so the
# two benchmarks cannot disagree about what a field is called or what it means.
if (!is.null(profile)) {
  EMC2:::.emc_profile_report(profile, elapsed = result$elapsed)
}
