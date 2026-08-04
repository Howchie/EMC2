# Reproduces the nested-vs-flat scheduling benchmark quoted in R/parallel_pool.R
# and NEWS.md.
#
#   Rscript WorkingTests/bench_flat_parallel.R
#
# Each configuration runs in its own R process, so no configuration inherits
# another's RNG state, JIT state or heap.  Set EMC_PF to change the particle
# factor (work per task) and EMC_CORES to change the core grid.

lib <- Sys.getenv("EMC_LIB", .libPaths()[1])
self <- "WorkingTests/bench_flat_parallel.R"   # run from the package root
setup_file <- Sys.getenv("EMC_SETUP", file.path(tempdir(), "emc_bench_data.rds"))

if (identical(Sys.getenv("EMC_CHILD"), "")) {
  .libPaths(c(lib, .libPaths()))
  suppressMessages(library(EMC2))
  set.seed(1)
  n_subjects <- 24L
  # Deliberately uneven: equal-sized subjects hide the load imbalance that
  # flattened scheduling is supposed to remove.
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
  message(sprintf("simulated %d rows, %d subjects, %d-%d trials each",
                  nrow(sim), n_subjects, min(table(sim$subjects)), max(table(sim$subjects))))

  cores_grid <- as.integer(strsplit(Sys.getenv("EMC_CORES", "3,6,12"), ",")[[1]])
  n_chains <- 3L
  for (cores in cores_grid) {
    cfc <- min(n_chains, cores); cpc <- max(1L, cores %/% cfc)
    for (flat in c("FALSE", "TRUE")) {
      system2("Rscript", c(self, flat, cfc, cpc),
              env = c("EMC_CHILD=1", paste0("EMC_LIB=", lib),
                      paste0("EMC_SETUP=", setup_file),
                      paste0("EMC_PF=", Sys.getenv("EMC_PF", "30"))))
    }
  }
  quit(save = "no")
}

# --- child: one configuration ------------------------------------------------
.libPaths(c(lib, .libPaths()))
suppressMessages(library(EMC2))
a <- commandArgs(TRUE)
flat <- as.logical(a[1]); cfc <- as.integer(a[2]); cpc <- as.integer(a[3])
pf <- as.integer(Sys.getenv("EMC_PF", "30"))
d <- readRDS(setup_file)
options(emc2.flat_parallel = flat)
RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
emc <- suppressMessages(make_emc(d$sim, d$des, n_chains = 3))
emc <- suppressMessages(run_emc(emc, stage = "preburn",
  stop_criteria = list(iter = 10, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = cfc, cores_per_chain = cpc, verbose = FALSE,
  particle_factor = pf, step_size = 10, max_tries = 1))
t <- system.time(emc <- suppressMessages(run_emc(emc, stage = "burn",
  stop_criteria = list(iter = 40, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = cfc, cores_per_chain = cpc, verbose = FALSE,
  particle_factor = pf, step_size = 40, max_tries = 1)))
cpu <- t[["user.self"]] + t[["sys.self"]] + t[["user.child"]] + t[["sys.child"]]
cat(sprintf("flat=%-5s cores=%2d (cfc=%d,cpc=%d) elapsed=%7.1fs cpu=%8.1fs util=%.2f\n",
            flat, cfc * cpc, cfc, cpc, t[["elapsed"]], cpu, cpu / (t[["elapsed"]] * cfc * cpc)))
