# Fresh-process fixed-transition paper fixtures. Environment selects model,
# library, host, concurrency, seed and output. Data generation has its own seed.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
kind <- Sys.getenv('CASE', 'EXG')
seed <- as.integer(Sys.getenv('SEED', '101'))
set.seed(123)
if (kind == 'EXG') {
  des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1,
                model = REXG, formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  dat <- make_data(log(c(mu = .4, sigma = .1, tau = .2)), des,
                   n_trials = 1e4, rt_resolution = .001)
  dat <- make_missing(dat, UC = gamlss.dist::qexGAUS(.9, mu = .4, sigma = .1, nu = .2),
                      rt_resolution = .001, verbose = FALSE)
} else if (kind == 'LNR') {
  Dfun <- function(dadm) factor(abs(as.numeric(dadm$S) - as.numeric(dadm$lR)))
  des <- design(factors = list(subjects = 1, S = 1:3), Rlevels = 1:3,
                functions = list(D = Dfun), formula = list(m ~ D, s ~ 1, t0 ~ 1), model = LNR)
  par <- sampled_pars(des); par[] <- c(-.5, .5, 1, log(.75), log(.2))
  dat <- make_data(par, des, n_trials = 1e4)
  dat <- make_missing(dat, LT = .41,
                      UC = function(d) quantile(d$rt[is.finite(d$rt)], probs = .9),
                      rt_resolution = .001, verbose = FALSE)
} else stop('unknown CASE')
dat <- readRDS(Sys.getenv('DATA_IN'))
data_file <- tempfile(); saveRDS(dat, data_file); data_md5 <- unname(tools::md5sum(data_file)); unlink(data_file)
if (Sys.getenv('ARCHIVE_COMPAT', '0') == '1') {
  for (nm in c('LT', 'UT', 'LC', 'UC')) {
    values <- unique(dat[[nm]])
    if (length(values) != 1L) stop('expected a single bound for ', nm)
    value <- values[[1L]]
    # Legacy check uses strict inequalities after millisecond rounding.
    # Likelihood preparation reads the unchanged bound columns first.
    if (is.finite(value)) value <- if (nm %in% c('UC', 'UT')) value + 1e-8 else max(0, value - 1e-8)
    attr(dat, nm) <- value
  }
}
set.seed(seed)
options(emc2.fit_host = Sys.getenv('HOST', '1') == '1',
        emc2.sampler_profile = Sys.getenv('PROFILE', '0') == '1')
emc <- make_emc(dat, des, type = 'single', n_chains = 3, rt_resolution = .001)
initial_seed <- .Random.seed
start <- proc.time()
fit_result <- suppressWarnings(fit(emc, verbose = FALSE, cores_per_chain = 1,
  cores_for_chains = as.integer(Sys.getenv('CHAINS', '3'))))
timing <- proc.time() - start
saveRDS(list(fit = fit_result, initial_seed = initial_seed, final_seed = .Random.seed,
             timing = timing, case = kind, seed = seed, data_md5 = data_md5,
             session = sessionInfo(), options = options()[c('emc2.fit_host','emc2.sampler_profile')]),
        Sys.getenv('OUT'))
cat('case', kind, 'seed', seed, 'timing', timing, '\n')
