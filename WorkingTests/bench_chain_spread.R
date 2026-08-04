# How unevenly do chains finish a block?  That spread is the whole budget any
# dynamic core-reallocation scheme has to reclaim.  Measured across consecutive
# blocks from a fresh, dispersed start -- the regime where stragglers are
# reported -- rather than from a converged state.
#
#   Rscript WorkingTests/bench_chain_spread.R

lib <- Sys.getenv("EMC_LIB", .libPaths()[1])
pf <- as.integer(Sys.getenv("EMC_PF", "100"))
n_subjects <- as.integer(Sys.getenv("EMC_SUBJ", "24"))
n_chains <- as.integer(Sys.getenv("EMC_CHAINS", "3"))
cpc <- as.integer(Sys.getenv("EMC_CPC", "4"))
step <- as.integer(Sys.getenv("EMC_STEP", "25"))

.libPaths(c(lib, .libPaths()))
suppressMessages(library(EMC2))
set.seed(1)
subj_trials <- round(seq(150, 450, length.out = n_subjects))
dat <- do.call(rbind, lapply(seq_len(n_subjects), function(i)
  data.frame(subjects = factor(i, levels = seq_len(n_subjects)),
             S = factor(sample(c("left", "right"), subj_trials[i], TRUE),
                        levels = c("left", "right")),
             R = factor(NA, levels = c("left", "right")), rt = NA_real_)))
ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
des <- design(data = dat, model = RDM, matchfun = function(d) d$S == d$lR,
              formula = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1),
              contrasts = list(v = list(lM = ADmat)), constants = c(s = log(1)))
pv <- c(v = 1.5, v_lMd = 0.8, B = log(1), t0 = log(0.2), A = log(0.3))
pmat <- matrix(rep(pv, n_subjects), nrow = n_subjects, byrow = TRUE,
               dimnames = list(levels(dat$subjects), names(pv)))
sim <- make_data(pmat, design = des, data = dat)

stamp <- file.path(tempdir(), "spread"); dir.create(stamp, showWarnings = FALSE)
options(emc2.chain_timer_dir = stamp)

RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
emc <- suppressMessages(make_emc(sim, des, n_chains = n_chains))
cat(sprintf("%d subjects, %d chains, cpc=%d, pf=%d, step=%d\n\n",
            n_subjects, n_chains, cpc, pf, step))

block <- function(stage, iter) {
  for (b in seq_len(iter %/% step)) {
    unlink(list.files(stamp, full.names = TRUE))
    t <- system.time(emc <<- suppressMessages(run_emc(emc, stage = stage,
      stop_criteria = list(iter = step, max_gd = Inf, min_unique = 0, min_es = 0),
      cores_for_chains = n_chains, cores_per_chain = cpc, verbose = FALSE,
      particle_factor = pf, step_size = step, max_tries = 1)))
    ct <- sort(vapply(list.files(stamp, full.names = TRUE), readRDS, numeric(1)))
    cat(sprintf("%-7s block %d  wall=%5.1fs  chains=[%s]  spread=%4.1f%%  idle_core_s=%5.1f  master=%4.1fs\n",
                stage, b, t[["elapsed"]], paste(sprintf("%.1f", ct), collapse = " "),
                100 * (max(ct) - min(ct)) / max(ct),
                cpc * sum(max(ct) - ct), t[["elapsed"]] - max(ct)))
  }
}
block("preburn", 75)
block("burn", 75)
