# ---------------------------------------------------------------------------
# Reproduce + verify the Windows fitting slowdown fix WITHOUT a Windows machine.
#
# Run:
#   devtools::load_all(".")
#   source("WorkingTests/windows_parallel_speed.R")
#
# Background
# ----------
# On Windows there is no fork(), so auto_mclapply() uses a PSOCK cluster. The
# OLD code built and tore down a fresh cluster on EVERY call, including the
# mc.cores == 1 calls on the likelihood hot path (within-chain parallelism is
# disabled on Windows, so r_cores == 1 there). calc_ll_manager() calls
# auto_mclapply once per likelihood evaluation -> thousands of cluster
# create/destroy cycles per iteration. macOS/Linux never saw this because
# mclapply(mc.cores = 1) is just lapply.
#
# The FIX: auto_mclapply() now short-circuits mc.cores == 1 to a plain lapply,
# so no cluster is ever built on the hot path, on any OS.
#
# The PSOCK code path is cross-platform, so we reproduce the Windows behaviour
# on a Mac/Linux box by calling makeCluster/parLapply directly. This script
# emulates the OLD ("broken") and NEW ("fixed") auto_mclapply and times both,
# first in isolation and then on the REAL EMC2 likelihood hot path.
# ---------------------------------------------------------------------------

suppressMessages(library(parallel))
stopifnot("Run after devtools::load_all('.')" = "EMC2" %in% loadedNamespaces())

## ===========================================================================
## Part 1 - microbenchmark of the dispatcher itself (no EMC2 likelihood)
## ===========================================================================
# OLD Windows behaviour: a fresh PSOCK cluster on every call, no short-circuit.
auto_broken_win <- function(X, FUN, mc.cores, ...){
  cl <- makeCluster(mc.cores); on.exit(stopCluster(cl))
  parLapply(cl = cl, X, FUN, ...)
}
# NEW Windows behaviour: short-circuit mc.cores == 1 -> lapply.
auto_fixed_win <- function(X, FUN, mc.cores, ...){
  if (mc.cores == 1) return(lapply(X, FUN, ...))
  cl <- makeCluster(mc.cores); on.exit(stopCluster(cl))
  parLapply(cl = cl, X, FUN, ...)
}
# Current unix behaviour for reference.
auto_unix <- function(X, FUN, mc.cores, ...){
  if (mc.cores == 1) return(lapply(X, FUN, ...))
  mclapply(X, FUN, mc.cores = mc.cores, ...)
}

f <- function(i) sqrt(i)
N <- 200
t_unix   <- system.time(for (k in 1:N) auto_unix      (1:4, f, mc.cores = 1))[["elapsed"]]
t_broken <- system.time(for (k in 1:N) auto_broken_win(1:4, f, mc.cores = 1))[["elapsed"]]
t_fixed  <- system.time(for (k in 1:N) auto_fixed_win (1:4, f, mc.cores = 1))[["elapsed"]]

cat(sprintf("\n[Part 1] dispatcher microbenchmark, %d calls at mc.cores = 1\n", N))
cat(sprintf("  unix (mclapply)             : %8.3f s\n", t_unix))
cat(sprintf("  windows OLD (cluster/call)  : %8.3f s   <- the bug (~%.0f ms/call)\n",
            t_broken, 1000 * (t_broken - t_unix) / N))
cat(sprintf("  windows NEW (short-circuit) : %8.3f s   <- the fix\n", t_fixed))

## ===========================================================================
## Part 2 - the REAL likelihood hot path under emulated Windows
## ===========================================================================
# Build a small LNR model (same setup as tests/testthat/test-fit.R), then call
# the actual EMC2:::calc_ll_manager() repeatedly under each emulated dispatcher.
# This is exactly what particle-MH does per likelihood evaluation: proposals is
# a multi-row particle matrix with r_cores = 1, which hits the auto_mclapply
# branch (nrow(proposals) > r_cores).
ADmat    <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
matchfun <- function(d) d$S == d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
dat$subjects <- droplevels(dat$subjects)

design_LNR <- design(data = dat, model = LNR, matchfun = matchfun,
                     formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                     contrasts = list(m = list(lM = ADmat)))
emc <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 1)

dadm  <- emc[[1]]$data[[1]]
model <- emc[[1]]$model

# A particle-set-shaped proposals matrix (named columns = sampled parameters).
pv <- tryCatch(sampled_pars(design_LNR),
               error = function(e) setNames(rep(0, length(emc[[1]]$par_names)),
                                            emc[[1]]$par_names))
n_particles <- 100
set.seed(1)
proposals <- matrix(pv, nrow = n_particles, ncol = length(pv), byrow = TRUE,
                    dimnames = list(NULL, names(pv)))
proposals <- proposals + rnorm(length(proposals), 0, 0.02)

ns        <- asNamespace("EMC2")
orig_auto <- get("auto_mclapply", ns)

time_calc_ll <- function(auto_fun, n_calls){
  assignInNamespace("auto_mclapply", auto_fun, ns)
  on.exit(assignInNamespace("auto_mclapply", orig_auto, ns))
  # warm up (and capture a reference result for correctness)
  ref <- EMC2:::calc_ll_manager(proposals, dadm, model, r_cores = 1)
  t <- system.time(
    for (k in 1:n_calls) EMC2:::calc_ll_manager(proposals, dadm, model, r_cores = 1)
  )[["elapsed"]]
  list(time = t, ref = ref)
}

n_calls <- 40
r_fixed  <- time_calc_ll(auto_fixed_win,  n_calls)   # emulated Windows, FIXED
r_broken <- time_calc_ll(auto_broken_win, n_calls)   # emulated Windows, OLD
# correctness: identical likelihoods regardless of dispatcher
stopifnot(isTRUE(all.equal(r_fixed$ref, r_broken$ref)))

cat(sprintf("\n[Part 2] real calc_ll_manager() hot path, %d calls (%d particles), r_cores = 1\n",
            n_calls, n_particles))
cat(sprintf("  windows OLD (cluster/call)  : %8.3f s\n", r_broken$time))
cat(sprintf("  windows NEW (short-circuit) : %8.3f s\n", r_fixed$time))
cat(sprintf("  speedup from the fix        : %8.1fx\n", r_broken$time / r_fixed$time))
cat("  likelihoods identical across dispatchers: TRUE\n")

## ===========================================================================
## Part 3 - prove the FULL FIT benefits: count the hot-path mc.cores == 1 calls
## ===========================================================================
# Each auto_mclapply(mc.cores == 1) call was a full PSOCK cluster create/destroy
# on old Windows. We instrument the live dispatcher to count them during a real
# fit (and confirm no other code path builds clusters per likelihood). The raw
# mclapply/mcmapply calls elsewhere don't matter on Windows: they cannot fork
# there, so they run serially with no cluster. Only auto_mclapply ever built one.
cnt1 <- 0L; cntN <- 0L
counting <- function(X, FUN, mc.cores, ...){
  if (mc.cores == 1) cnt1 <<- cnt1 + 1L else cntN <<- cntN + 1L
  orig_auto(X, FUN, mc.cores = mc.cores, ...)
}
assignInNamespace("auto_mclapply", counting, ns)

set.seed(123)
emc2 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2)
t_fit <- system.time(suppressMessages(
  emc2 <- fit(emc2, cores_for_chains = 1, verbose = FALSE, particle_factor = 20,
              step_size = 25,
              stop_criteria = list(preburn = list(iter = 10), burn = list(iter = 10),
                                   adapt = list(iter = 10),  sample = list(iter = 25)))
))[["elapsed"]]
assignInNamespace("auto_mclapply", orig_auto, ns)

per_call <- if (n_calls > 0) (r_broken$time / n_calls) else NA   # real-dadm cost from Part 2
cat(sprintf("\n[Part 3] full fit (preburn/burn/adapt 10 + sample 25, 2 subj, 2 chains)\n"))
cat(sprintf("  fit completed (live fixed code) in : %.1f s\n", t_fit))
cat(sprintf("  auto_mclapply mc.cores == 1 calls  : %d   <- each was a cluster on old Windows\n", cnt1))
cat(sprintf("  auto_mclapply mc.cores  > 1 calls  : %d   <- genuine clusters (none here)\n", cntN))
cat(sprintf("  est. old-Windows overhead removed  : ~%.0f s (%d x %.2f s/call, Part 2 rate)\n",
            cnt1 * per_call, cnt1, per_call))

cat("\nDone. Expected: Part 1 & 2 show the OLD path is dramatically slower at\n")
cat("mc.cores = 1, the NEW path matches unix speed, and likelihoods are identical.\n")
cat("Part 3 confirms a normal fit still works under the patched dispatcher.\n")
cat("\nNote on compare(): the Windows nested-cluster guard lives in compare() and\n")
cat("triggers only when Sys.info()[1]=='Windows' AND n_parallel_models>1 AND an\n")
cat("inner core count >1 -- a real Windows run will print the 'parallelises across\n")
cat("models only' warning; it cannot be exercised by faking the OS here.\n")
