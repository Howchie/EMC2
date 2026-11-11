# Two-choice REXG native likelihood fixture for `rexg_logS_at_t` caching.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
data('forstmann', package = 'EMC2')

lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
matchfun <- function(d) d$S == d$lR
subjects <- unique(forstmann$subjects)
dat <- forstmann[forstmann$subjects %in% subjects[[1L]], , drop = FALSE]
dat$subjects <- droplevels(dat$subjects)
replicate_data <- as.integer(Sys.getenv('DATA_REPS', '8'))
if (!is.finite(replicate_data) || replicate_data < 1L)
  stop('DATA_REPS must be a positive integer')
dat <- do.call(rbind, replicate(replicate_data, dat, simplify = FALSE))
rownames(dat) <- NULL

des <- design(
  data = dat, model = REXG, matchfun = matchfun,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
  functions = list(lI = lI_fun), UT = 2)
emc <- make_emc(dat, des, type = 'single', n_chains = 1, compress = FALSE)
dadm <- emc[[1L]]$data[[1L]]
model <- emc[[1L]]$model()
p_types <- names(model$p_types)
designs <- lapply(p_types, function(nm) {
  d <- attr(dadm, 'designs')[[nm]]
  d[attr(d, 'expand'), , drop = FALSE]
})
names(designs) <- p_types

# Use one shared parameter tuple across the two racers so the endpoint callback
# can reuse the same zero-truncation normalization within each particle.
p <- stats::setNames(rep(-5, length(p_types)), p_types)
sampled_names <- rownames(emc[[1L]]$samples$alpha)
p[sampled_names] <- as.numeric(emc[[1L]]$samples$alpha[, 1L, 1L])
particle_count <- as.integer(Sys.getenv('N_PARTICLES', '256'))
repeat_count <- as.integer(Sys.getenv('REPS', '12'))
particles <- matrix(rep(p, each = particle_count), nrow = particle_count,
                    dimnames = list(NULL, p_types))
call <- function() EMC2:::calc_ll_oo(
  particles, dadm, constants = attr(dadm, 'constants'), designs = designs,
  type = model$c_name, bounds = model$bound, transforms = model$transform,
  pretransforms = model$pre_transform, p_types = p_types,
  min_ll = log(1e-10), trend = model$trend)

result <- call()
timing <- replicate(repeat_count, system.time(invisible(call()))[['elapsed']])
output <- Sys.getenv('OUT')
if (nzchar(output)) {
  data_file <- tempfile()
  saveRDS(dadm, data_file)
  data_md5 <- unname(tools::md5sum(data_file))
  unlink(data_file)
  saveRDS(list(result = result, timing = timing, data_md5 = data_md5,
               rows = nrow(dadm), particles = particle_count,
               session = sessionInfo()), output)
}
cat('rows', nrow(dadm), 'particles', particle_count,
    'median_seconds', median(timing), 'timings', paste(timing, collapse = ','),
    'data_md5', data_md5, '\n')
