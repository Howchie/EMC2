suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
set.seed(123)
des <- design(model = LBA, matchfun = function(d) d$S == d$lR,
              formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
              constants = c(sv = log(1)), pooled = 't0', report_p_vector = FALSE,
              factors = list(subjects = 1:6, S = c('left', 'right')),
              Rlevels = c('left', 'right'))
truth <- c(v = 1, v_lMTRUE = 1.5, B = log(1), A = log(.4), t0 = log(.25))
re <- make_random_effects(des, truth, n_subj = 6)
dat <- make_data(re, des, n_trials = 80)
file <- tempfile(); saveRDS(dat, file); hash <- unname(tools::md5sum(file)); unlink(file)
set.seed(as.integer(Sys.getenv('SEED', '101')))
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE)
sample_iter <- as.integer(Sys.getenv('SAMPLE_ITER', '400'))
emc <- make_emc(dat, des, n_chains = 3)
start <- proc.time()
out <- suppressWarnings(fit(emc, verbose = FALSE, cores_per_chain = 2, cores_for_chains = 2,
  step_size = 100, stop_criteria = list(preburn = list(iter = 50), burn = list(iter = 200),
    adapt = list(iter = 200, min_unique = 1), sample = list(iter = sample_iter))))
timing <- proc.time() - start
profile <- function(name) {
  ns <- asNamespace('EMC2')
  if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
  get(name, envir = ns, inherits = FALSE)()
}
saveRDS(list(fit = out, seed = .Random.seed, timing = timing, data_md5 = hash,
             sample_iter = sample_iter,
             blocks = profile('.emc_block_profile'), iterations = profile('.emc_iteration_profile'),
             session = sessionInfo()), Sys.getenv('OUT'))
cat('timing', timing, '\n')
