# Where does a sampler iteration go when the subject worker pool is running?
#
#   Rscript WorkingTests/bench_worker_pool.R
#   NSUBJ=32 CORES=32 WIDE=1 Rscript WorkingTests/bench_worker_pool.R
#
# Reports the per-iteration split from the sampler profile option, through the
# one reporting path in R/profile_schema.R.  The two numbers to read together
# are `response_wait` and `worker_cpu_max` -- the slowest *worker's* whole
# assigned partition, not the slowest single subject, which is what this script
# used to print under the name `worker_max`.  When they are close the pool is
# bound by the work, and when `request_send` is a large share it is bound by
# the wire instead.  The shared group state is written once to a file; each
# request carries only its path.
#
# Two traps, both of which have cost real time here:
#
#  * `nproc` reports 1 in some harnesses purely because `OMP_NUM_THREADS=1` is
#    exported. `parallel::detectCores()` is the number that governs this pool.
#    Check it (printed below) before concluding a pool "does not engage here".
#  * Do NOT A/B an R-only change by installing two temp libraries. Separate
#    `R CMD INSTALL` runs differ in lazy-load database and bytecode layout by
#    several percent of whole-fit CPU, which is larger than most changes worth
#    measuring. Patch the namespace in-process instead (`deparse` the function,
#    substitute, `eval`, `assignInNamespace`) so exactly one thing varies.

lib <- Sys.getenv("EMC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))
cat("namespace   :", dirname(getNamespaceInfo("EMC2", "path")), "\n")
cat("detectCores :", parallel::detectCores(), "\n")

NS    <- as.integer(Sys.getenv("NSUBJ",  "32"))
NT    <- as.integer(Sys.getenv("NTRIAL", "200"))
NIT   <- as.integer(Sys.getenv("NITER",  "40"))
CORES <- as.integer(Sys.getenv("CORES",  as.character(parallel::detectCores())))
WIDE  <- as.integer(Sys.getenv("WIDE",   "1"))

set.seed(11)
base <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
base <- base[rep(seq_len(nrow(base)), length.out = NT), , drop = FALSE]
dat <- do.call(rbind, lapply(seq_len(NS), function(s) {
  d <- base; d$subjects <- sprintf("s%02d", s); d
}))
dat$subjects <- factor(dat$subjects); rownames(dat) <- NULL
# Extra balanced factors exist only to widen the parameter vector: the shared
# broadcast is O(P^2), so P is what this benchmark is really varying.
dat$F1 <- factor(rep_len(paste0("f", 1:5), nrow(dat)))
dat$F2 <- factor(rep_len(paste0("g", 1:4), nrow(dat)))

matchfun <- function(d) d$S == d$lR
form <- if (WIDE) {
  list(v ~ lM * E * F1, B ~ E * lR * F2, A ~ 1, t0 ~ 1)
} else {
  list(v ~ lM * E, B ~ E * lR, A ~ 1, t0 ~ 1)
}
des <- design(data = dat, model = LBA, matchfun = matchfun, formula = form,
              constants = c(sv = log(1)))
emc <- make_emc(dat, des, n_chains = 1, compress = TRUE)
cat(sprintf("subjects=%d  trials/subj=%d  P=%d  iters=%d  cores=%d\n",
            NS, NT, length(sampled_pars(des)), NIT, CORES))

# `_expensive` additionally sizes each serialised request and reads the process
# tree's private memory; both are opt-in because measuring them costs an extra
# object walk and a /proc read per iteration.
options(emc2.sampler_profile = TRUE, emc2.sampler_profile_expensive = TRUE)
el <- system.time(
  emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = NIT),
                 cores_for_chains = 1, cores_per_chain = CORES,
                 verbose = FALSE, verboseProgress = FALSE))

# One schema, one formatter: everything printed below is declared in
# R/profile_schema.R, so a field added there appears here without this script
# learning its name, its units or where it nests.
EMC2:::.emc_profile_report(attr(emc[[1]]$samples, "sampler_profile"),
                           elapsed = el[["elapsed"]])
