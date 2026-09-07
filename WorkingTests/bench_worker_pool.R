# Where does a sampler iteration go when the subject worker pool is running?
#
#   Rscript WorkingTests/bench_worker_pool.R
#   NSUBJ=32 CORES=32 WIDE=1 Rscript WorkingTests/bench_worker_pool.R
#
# Reports the per-iteration split from `options(emc2.sampler_profile = TRUE)`.
# The two numbers to read together are `response_wait` and `worker_max`: when
# they are close the pool is bound by the work, and when `request_send` is a
# large share it is bound by the wire instead.  The shared group state is now
# written once to a file; each request carries only its path.
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

options(emc2.sampler_profile = TRUE)
el <- system.time(
  emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = NIT),
                 cores_for_chains = 1, cores_per_chain = CORES,
                 verbose = FALSE, verboseProgress = FALSE))

pr <- attr(emc[[1]]$samples, "sampler_profile")
if (is.null(pr) || !nrow(pr)) stop("no profile recorded")
pr <- pr[-1, , drop = FALSE]          # the first iteration pays pool start-up
f <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
tot <- f(pr$total)

cat(sprintf("\nworkers          %d\n", max(pr$workers)))
cat(sprintf("wall %.2f s   iteration %.2f ms\n", el[["elapsed"]], 1000 * tot))
for (nm in c("gibbs", "group_cache", "particle", "fill",
             "shared_serialize", "request_send", "response_wait")) {
  v <- f(pr[[nm]])
  cat(sprintf("  %-17s %8.2f ms  %5.1f%%\n", nm, 1000 * v, 100 * v / tot))
}
cat(sprintf("  %-17s %8.2f ms  (slowest worker's own work)\n",
            "worker_max", 1000 * f(pr$worker_max)))
cat(sprintf("  %-17s %8.2f ms  (summed over workers)\n",
            "worker_sum", 1000 * f(pr$worker_sum)))
cat(sprintf("\nshared file        %9.0f B\n", f(pr$shared_bytes)))
cat(sprintf("private per iter   %9.0f B\n", f(pr$private_bytes)))
cat(sprintf("wire per iteration %9.0f B  (%.2f MB)\n",
            f(pr$wire_bytes), f(pr$wire_bytes) / 1024^2))
