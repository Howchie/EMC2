# REXG raw-kernel fixture with four exact parameter tuples across design cells.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
set.seed(123)
des <- design(factors = list(subjects = 1, S = 1:4), Rlevels = 1:4,
              model = REXG,
              formula = list(mu ~ S, sigma ~ 1, tau ~ 1))
p <- c(mu = log(.4), mu_S2 = log(1.15), mu_S3 = log(.85),
       mu_S4 = log(1.3), sigma = log(.1), tau = log(.2))
dat <- make_data(p, des, n_trials = 4000, rt_resolution = .001)
dat <- make_missing(dat, UC = function(d) quantile(d$rt[is.finite(d$rt)], .9),
                    rt_resolution = .001, verbose = FALSE)
emc <- make_emc(dat, des, type = 'single', n_chains = 1, rt_resolution = .001)
dadm <- emc[[1]]$data[[1]]; model <- emc[[1]]$model()
p_types <- names(model$p_types)
designs <- lapply(p_types, function(nm) {
  d <- attr(dadm, 'designs')[[nm]]
  d[attr(d, 'expand'), , drop = FALSE]
})
names(designs) <- p_types
constants <- attr(dadm, 'constants'); if (is.null(constants)) constants <- NA
particles <- matrix(rep(p, each = 256), nrow = 256, dimnames = list(NULL, names(p)))
fun <- get('calc_ll_oo', asNamespace('EMC2'))
call <- function() fun(particles, dadm, constants, designs, model$c_name, model$bound,
                        model$transform, model$pre_transform, p_types, log(1e-10), model$trend)
result <- call()
timing <- replicate(10L, system.time(invisible(call()))[['elapsed']])
out <- Sys.getenv('OUT')
if (nzchar(out)) {
  data_file <- tempfile(); saveRDS(dadm, data_file)
  data_md5 <- unname(tools::md5sum(data_file)); unlink(data_file)
  saveRDS(list(result = result, timing = timing, data_md5 = data_md5,
               rows = nrow(dadm), session = sessionInfo()), out)
}
cat('rows', nrow(dadm), 'particles', nrow(particles), 'times', timing, '\n')
