# Two measurements that decide whether dynamic cores_per_chain reallocation can
# pay off at all:
#
#   (1) does cores_per_chain actually buy speed?  If a chain cannot use extra
#       cores, donating a finished chain's cores to it gains nothing.
#   (2) how far apart do chains finish within a block?  That spread is the
#       entire budget a dynamic scheme has to work with.
#
#   Rscript WorkingTests/bench_cpc_scaling.R
#
# EMC_PF   particle factor (work per subject-iteration)
# EMC_SUBJ number of subjects
# EMC_CPC  comma-separated cores_per_chain grid

lib <- Sys.getenv("EMC_LIB", .libPaths()[1])
self <- "WorkingTests/bench_cpc_scaling.R"
setup_file <- Sys.getenv("EMC_SETUP", file.path(tempdir(), "emc_cpc_data.rds"))
pf <- as.integer(Sys.getenv("EMC_PF", "30"))
n_subjects <- as.integer(Sys.getenv("EMC_SUBJ", "24"))

if (identical(Sys.getenv("EMC_CHILD"), "")) {
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
  saveRDS(list(sim = sim, des = des), setup_file)
  message(sprintf("%d subjects, %d-%d trials, particle_factor=%d",
                  n_subjects, min(table(sim$subjects)), max(table(sim$subjects)), pf))

  for (cpc in as.integer(strsplit(Sys.getenv("EMC_CPC", "1,2,4,8"), ",")[[1]]))
    system2("Rscript", c(self, cpc),
            env = c("EMC_CHILD=1", paste0("EMC_LIB=", lib),
                    paste0("EMC_SETUP=", setup_file), paste0("EMC_PF=", pf),
                    paste0("EMC_SUBJ=", n_subjects)))
  quit(save = "no")
}

# --- child: one cores_per_chain setting --------------------------------------
.libPaths(c(lib, .libPaths()))
suppressMessages(library(EMC2))
cpc <- as.integer(commandArgs(TRUE)[1])
d <- readRDS(setup_file)
RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
emc <- suppressMessages(make_emc(d$sim, d$des, n_chains = 3))
emc <- suppressMessages(run_emc(emc, stage = "preburn",
  stop_criteria = list(iter = 10, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = 3, cores_per_chain = cpc, verbose = FALSE,
  particle_factor = pf, step_size = 10, max_tries = 1))

# Per-chain wall time for one block, to expose the straggler spread.  Each
# chain writes its own duration; the block's cost is the slowest of them.
stamp <- file.path(tempdir(), "chain_times")
unlink(list.files(stamp, full.names = TRUE)); dir.create(stamp, showWarnings = FALSE)
options(emc2.chain_timer_dir = stamp)

t <- system.time(emc <- suppressMessages(run_emc(emc, stage = "burn",
  stop_criteria = list(iter = 40, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = 3, cores_per_chain = cpc, verbose = FALSE,
  particle_factor = pf, step_size = 40, max_tries = 1)))
cpu <- sum(t[c("user.self", "sys.self", "user.child", "sys.child")])
ct <- sort(vapply(list.files(stamp, full.names = TRUE), readRDS, numeric(1)))
cat(sprintf("cpc=%d total=%2d elapsed=%6.1fs cpu=%7.1fs util=%.2f",
            cpc, 3 * cpc, t[["elapsed"]], cpu, cpu / (t[["elapsed"]] * 3 * cpc)))
if (length(ct)) cat(sprintf("  chains=[%s] spread=%.0f%%",
                            paste(sprintf("%.1f", ct), collapse = " "),
                            100 * (max(ct) - min(ct)) / max(ct)))
cat("\n")
