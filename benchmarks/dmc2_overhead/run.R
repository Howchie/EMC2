# Run a fixture in a fresh R process: Rscript run.R fixture.R output.rds seed
# The fixture defines run_fixture(seed), returning the fit and any measurements.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('usage: run.R fixture.R output.rds seed')
seed <- suppressWarnings(as.integer(args[[3L]]))
if (is.na(seed)) stop('seed must be an integer')
fixture <- normalizePath(args[[1L]], mustWork = TRUE)
source(fixture, local = .GlobalEnv)
if (!exists('run_fixture', mode = 'function', inherits = FALSE))
  stop('fixture must define run_fixture(seed)')
set.seed(seed)
initial_rng <- .Random.seed
started <- Sys.time()
timing <- system.time(result <- run_fixture(seed))
record <- list(fixture = fixture, fixture_md5 = unname(tools::md5sum(fixture)),
               seed = seed, initial_rng = initial_rng,
               final_rng = if (exists('.Random.seed', .GlobalEnv, inherits = FALSE))
                 get('.Random.seed', .GlobalEnv) else NULL,
               started = started, timing = timing, result = result,
               session = sessionInfo(),
               environment = Sys.getenv(c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                          'MKL_NUM_THREADS', 'R_LIBS', 'R_LIBS_USER')),
               options = options()[c('emc2.sampler_profile', 'emc2.fit_host')])
saveRDS(record, args[[2L]])
