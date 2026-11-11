# Fixed-transition one-chain REXG race fit with a finite upper truncation.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
data('forstmann', package = 'EMC2')

lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
matchfun <- function(d) d$S == d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ,
                 drop = FALSE]
dat$subjects <- droplevels(dat$subjects)
replicate_data <- as.integer(Sys.getenv('DATA_REPS', '2'))
dat <- do.call(rbind, replicate(replicate_data, dat, simplify = FALSE))
rownames(dat) <- NULL
n_racers <- as.integer(Sys.getenv('N_RACERS', '2'))
response_levels <- levels(dat$R)
if (length(response_levels) < n_racers) {
  response_levels <- c(response_levels,
    paste0('extra', seq_len(n_racers - length(response_levels))))
  dat$R <- factor(as.character(dat$R), levels = response_levels)
  dat$S <- factor(as.character(dat$S), levels = response_levels)
}
upper_bound <- if (identical(Sys.getenv('NONUNIFORM_UT', '0'), '1'))
  rep(c(1.5, 2), length.out = nrow(dat)) else 2
des <- design(
  data = dat, model = REXG, matchfun = matchfun, Rlevels = response_levels,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
  functions = list(lI = lI_fun), UT = upper_bound
)

seed <- as.integer(Sys.getenv('SEED', '101'))
set.seed(seed)
emc <- make_emc(dat, des, type = 'single', n_chains = 1, compress = FALSE)
initial_seed <- .Random.seed
data_file <- tempfile()
saveRDS(dat, data_file)
data_md5 <- unname(tools::md5sum(data_file))
unlink(data_file)
options(emc2.fit_host = FALSE, emc2.sampler_profile = FALSE,
        emc2.persistent_chains = FALSE)
started <- proc.time()
result <- suppressWarnings(fit(
  emc, verbose = FALSE, cores_per_chain = 1, cores_for_chains = 1,
  step_size = 10,
  stop_criteria = list(
    preburn = list(iter = 10), burn = list(iter = 10),
    adapt = list(iter = 10, min_unique = 1), sample = list(iter = 40)
  )
))
timing <- proc.time() - started
out <- Sys.getenv('OUT')
saveRDS(list(fit = result, initial_seed = initial_seed,
             final_seed = .Random.seed, timing = timing, seed = seed,
             data_md5 = data_md5, session = sessionInfo()), out)
cat('seed', seed, 'rows', nrow(dat), 'elapsed', timing[['elapsed']],
    'data_md5', data_md5, '\n')
