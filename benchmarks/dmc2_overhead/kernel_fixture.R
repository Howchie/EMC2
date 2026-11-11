# Full native likelihood call, including setup, per 704-row particle workload.
# DMC2 REXG uses an untruncated distribution; timings are not a correctness comparison.
pkg <- Sys.getenv('PKG', 'EMC2')
suppressPackageStartupMessages(library(pkg, character.only = TRUE, lib.loc = Sys.getenv('EMC_LIB')))
ns <- asNamespace(pkg)
set.seed(123)
des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1,
              model = REXG, formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
p <- log(c(mu = .4, sigma = .1, tau = .2))
dat <- make_data(p, des, n_trials = 1e4, rt_resolution = .001)
dat <- make_missing(dat, UC = gamlss.dist::qexGAUS(.9, mu = .4, sigma = .1, nu = .2),
                    rt_resolution = .001, verbose = FALSE)
emc <- make_emc(dat, des, type = 'single', n_chains = 1, rt_resolution = .001)
dadm <- emc[[1]]$data[[1]]; model <- emc[[1]]$model()
p_types <- names(model$p_types)
designs <- lapply(p_types, function(nm) { d <- attr(dadm, 'designs')[[nm]]; d[attr(d, 'expand'), , drop = FALSE] })
names(designs) <- p_types
constants <- attr(dadm, 'constants'); if (is.null(constants)) constants <- NA
particle_count <- as.integer(Sys.getenv('N_PARTICLES', '256'))
repeat_count <- as.integer(Sys.getenv('REPS', '10'))
stopifnot(length(particle_count) == 1L, is.finite(particle_count),
          particle_count > 0L, length(repeat_count) == 1L,
          is.finite(repeat_count), repeat_count > 0L)
particles <- matrix(rep(p, each = particle_count), nrow = particle_count,
                    dimnames = list(NULL, names(p)))
fun <- get(if (pkg == 'EMC2') 'calc_ll_oo' else 'calc_ll', ns)
call <- function() fun(particles, dadm, constants, designs, model$c_name, model$bound,
                       model$transform, model$pre_transform, p_types, log(1e-10), model$trend)
result <- call()
timing <- replicate(repeat_count, system.time(invisible(call()))[['elapsed']])
out <- Sys.getenv('OUT')
if (nzchar(out)) {
  data_file <- tempfile(); saveRDS(dadm, data_file)
  data_md5 <- unname(tools::md5sum(data_file)); unlink(data_file)
  saveRDS(list(result = result, timing = timing, data_md5 = data_md5,
               session = sessionInfo()), out)
}
cat('package', pkg, 'rows', nrow(dadm), 'particles', nrow(particles), 'times', timing,
    'seconds per row-particle', median(timing) / (nrow(dadm) * nrow(particles)), '\n')
