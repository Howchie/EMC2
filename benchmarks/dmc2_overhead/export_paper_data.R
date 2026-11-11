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
saveRDS(dat, Sys.getenv('OUT'))
cat('rows', nrow(dat), '\n')
