# Analytic bimodal target with exact 50/50 mode mass, using custom likelihood.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
custom_ll <- function(pars, dadm, ...) {
  x <- pars[1]
  a <- dnorm(x, -3, .5, log = TRUE); b <- dnorm(x, 3, .5, log = TRUE)
  m <- max(a, b)
  m + log((exp(a - m) + exp(b - m)) / 2) - dnorm(x, log = TRUE)
}
des <- design(model = custom_ll, custom_p_vector = 'theta')
dat <- get_data(samples_LNR); dat <- dat[dat$subjects == levels(dat$subjects)[1], , drop = FALSE]
set.seed(as.integer(Sys.getenv('SEED', '101')))
sample_draws <- as.integer(Sys.getenv('DRAWS', '400'))
if (length(sample_draws) != 1L || is.na(sample_draws) || sample_draws < 1L)
  stop('DRAWS must be a positive integer')
cache_max_reuses <- Sys.getenv('CACHE_MAX_REUSES', 'default')
if (cache_max_reuses != 'default')
  options(emc2.subject_fit_cache_max_reuses = as.numeric(cache_max_reuses))
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE)
emc <- make_emc(dat, des, type = 'single', n_chains = 3)
start <- proc.time()
out <- suppressWarnings(fit(emc, verbose = FALSE, cores_for_chains = 3, cores_per_chain = 1,
  stop_criteria = list(preburn = list(iter = 50), burn = list(iter = 200),
                       adapt = list(iter = 200, min_unique = 1), sample = list(iter = sample_draws)), step_size = 100))
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, sampler_seed = as.integer(Sys.getenv('SEED', '101')),
             sample_draws = sample_draws,
             cache_max_reuses = getOption('emc2.subject_fit_cache_max_reuses'),
             timing = timing,
             blocks = EMC2:::.emc_block_profile(), session = sessionInfo()), Sys.getenv('OUT'))
cat('time', timing, 'fraction right mode', mean(unlist(lapply(out, function(ch) ch$samples$alpha)) > 0), '\n')
