# Full native likelihood call, including setup, per 704-row particle workload.
# DMC2 REXG uses an untruncated distribution; timings are not a correctness comparison.
pkg <- Sys.getenv('PKG', 'EMC2')
suppressPackageStartupMessages(library(pkg, character.only = TRUE, lib.loc = Sys.getenv('EMC_LIB')))
ns <- asNamespace(pkg)
set.seed(123)
Dfun <- function(dadm) factor(abs(as.numeric(dadm$S) - as.numeric(dadm$lR)))
des <- design(factors = list(subjects = 1, S = 1:3), Rlevels = 1:3,
              functions = list(D = Dfun), formula = list(m ~ D, s ~ 1, t0 ~ 1), model = LNR)
p <- sampled_pars(des); p[] <- c(-.5, .5, 1, log(.75), log(.2))
dat <- make_data(p, des, n_trials = 1e4)
dat <- make_missing(dat, LT = .41, UC = function(d) quantile(d$rt[is.finite(d$rt)], .9),
                    rt_resolution = .001, verbose = FALSE)
emc <- make_emc(dat, des, type = 'single', n_chains = 1, rt_resolution = .001)
dadm <- emc[[1]]$data[[1]]; model <- emc[[1]]$model()
p_types <- names(model$p_types)
designs <- lapply(p_types, function(nm) { d <- attr(dadm, 'designs')[[nm]]; d[attr(d, 'expand'), , drop = FALSE] })
names(designs) <- p_types
constants <- attr(dadm, 'constants'); if (is.null(constants)) constants <- NA
particles <- matrix(rep(p, each = 256), nrow = 256, dimnames = list(NULL, names(p)))
fun <- get(if (pkg == 'EMC2') 'calc_ll_oo' else 'calc_ll', ns)
call <- function() fun(particles, dadm, constants, designs, model$c_name, model$bound,
                       model$transform, model$pre_transform, p_types, log(1e-10), model$trend)
invisible(call())
timing <- replicate(10L, system.time(invisible(call()))[['elapsed']])
cat('package', pkg, 'rows', nrow(dadm), 'particles', nrow(particles), 'times', timing,
    'seconds per row-particle', median(timing) / (nrow(dadm) * nrow(particles)), '\n')
