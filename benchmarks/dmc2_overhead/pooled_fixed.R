suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
assignInNamespace('.emc_alloc_solve_units', function(u) {
  hmc <- vapply(u$units, function(unit) identical(unit$kind, 'subject') && !is.null(unit$h), logical(1))
  x <- vapply(u$units, function(unit) if (identical(unit$kind, 'shear')) 1 else 0, numeric(1))
  x[hmc] <- .5
  list(x = x, q = ifelse(hmc, .5, 1), joint = if (!is.null(u$joint)) 1,
       rho = NA_real_, tau = NA_real_, cost = u$C_rest,
       objective = NA_real_, objective0 = NA_real_, tau0 = NA_real_)
}, asNamespace('EMC2'))
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
emc <- make_emc(dat, des, n_chains = 3)
start <- proc.time()
out <- suppressWarnings(fit(emc, verbose = FALSE, cores_per_chain = 2, cores_for_chains = 2,
  step_size = 100, stop_criteria = list(preburn = list(iter = 50), burn = list(iter = 200),
    adapt = list(iter = 200, min_unique = 1), sample = list(iter = 400))))
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, timing = timing, data_md5 = hash,
             blocks = EMC2:::.emc_block_profile(), iterations = EMC2:::.emc_iteration_profile(),
             session = sessionInfo()), Sys.getenv('OUT'))
cat('timing', timing, '\n')
