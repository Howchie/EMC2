# Compare real cached Laplace-refresh jobs on the adaptive serial and fork
# routes. The first refresh warms subject-fit/check caches; subsequent parallel
# refreshes calibrate the observed job and startup costs before paired timing.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
input <- Sys.getenv("FIXTURE")
output <- Sys.getenv("OUT")
seed <- as.integer(Sys.getenv("SEED", "101"))
if (!nzchar(input) || !file.exists(input) || !nzchar(output))
  stop("set FIXTURE and OUT to existing/desired RDS paths")

ns <- asNamespace("EMC2")
getf <- function(name) get(name, envir = ns, inherits = FALSE)
emc <- getf("restore_custom_kernel_pointers")(
  readRDS(input)$fit, quiet = TRUE)
emc <- getf("restore_duplicates")(emc)
queue <- getf(".emc_rebuild_queue")
queue$costs <- list()
queue$last_decision <- NULL
getf(".emc_warm_psis")()
set.seed(seed)

# This first pass runs genuine mode and importance-check work. Further passes
# use the exact conditioning cache to measure the small residual dispatch cost.
emc <- getf(".emc_laplace_refresh")(emc, 3L)
options(emc2.rebuild_serial_margin = 0.01)
for (i in seq_len(8L)) emc <- getf(".emc_laplace_refresh")(emc, 3L)
prediction <- getf(".emc_rebuild_prediction")("laplace_refresh", 3L, 3L)
if (is.null(prediction) ||
    prediction$serial >= 0.9 * prediction$parallel) {
  stop("calibrated cached fixture does not predict a serial win")
}

base <- emc
run_arm <- function(route) {
  options(emc2.rebuild_serial_margin = if (route == "serial") 0.9 else 0.01)
  set.seed(987L)
  started <- proc.time()[["elapsed"]]
  out <- getf(".emc_laplace_refresh")(base, 3L)
  elapsed <- proc.time()[["elapsed"]] - started
  decision <- queue$last_decision
  if (!identical(decision$route, route))
    stop(sprintf("requested %s route, got %s", route, decision$route))
  list(fit = out, elapsed = elapsed, caller_seed = .Random.seed,
       decision = decision)
}

sample_fields <- c("alpha", "stage", "subj_ll", "idx")
same_state <- function(a, b) {
  samples <- all(vapply(seq_along(a), function(j) all(vapply(
    sample_fields, function(field)
      identical(a[[j]]$samples[[field]], b[[j]]$samples[[field]]),
    logical(1))), logical(1)))
  rng <- all(vapply(seq_along(a), function(j)
    identical(a[[j]]$rng, b[[j]]$rng), logical(1)))
  fields <- setdiff(union(names(a[[1L]]), names(b[[1L]])),
                    c("data", "prior", "model"))
  continuation <- all(vapply(fields, function(field) all(vapply(
    seq_along(a), function(j)
      identical(a[[j]][[field]], b[[j]][[field]]), logical(1))), logical(1)))
  list(samples = samples, rng = rng, continuation = continuation)
}

order <- if (seed %% 2L) c("parallel", "serial") else c("serial", "parallel")
pairs <- vector("list", 10L)
for (i in seq_along(pairs)) {
  arms <- if (i %% 2L) order else rev(order)
  runs <- setNames(lapply(arms, run_arm), arms)
  parity <- same_state(runs$parallel$fit, runs$serial$fit)
  exact <- all(unlist(parity)) &&
    identical(runs$parallel$caller_seed, runs$serial$caller_seed)
  if (!exact) stop("serial and parallel cached refreshes diverged")
  pairs[[i]] <- list(order = arms,
                     parallel_seconds = runs$parallel$elapsed,
                     serial_seconds = runs$serial$elapsed,
                     parallel_prediction = runs$parallel$decision,
                     serial_prediction = runs$serial$decision,
                     exact = exact)
}

parallel_seconds <- vapply(pairs, `[[`, numeric(1), "parallel_seconds")
serial_seconds <- vapply(pairs, `[[`, numeric(1), "serial_seconds")
saveRDS(list(seed = seed, input_md5 = unname(tools::md5sum(input)),
             prediction = prediction, pairs = pairs,
             median_parallel = stats::median(parallel_seconds),
             median_serial = stats::median(serial_seconds),
             all_exact = all(vapply(pairs, `[[`, logical(1), "exact")),
             session = sessionInfo()), output)
cat("seed", seed, "median parallel", stats::median(parallel_seconds),
    "median serial", stats::median(serial_seconds),
    "speedup", stats::median(parallel_seconds) / stats::median(serial_seconds),
    "all exact", all(vapply(pairs, `[[`, logical(1), "exact")), "\n")
