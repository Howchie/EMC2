# The native iteration steps (src/iteration_step_*.cpp) against the R steps
# (R/iteration_steps.R): each stage block run both ways from the same state
# and random stream.

native_steps_match <- function(cores_per_chain) {
  old <- options(emc2.fit_host = FALSE, emc2.native_steps_trace = TRUE)
  on.exit(options(old), add = TRUE)
  assign("native_runs", NULL, envir = EMC2:::.emc_pool_state)
  des <- design(factors = list(subjects = 1:8, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(47)
  pars <- cbind(mu = log(0.4) + rnorm(8, 0, 0.15), sigma = log(0.05) + rnorm(8, 0, 0.3),
                tau = log(0.2) + rnorm(8, 0, 0.3))
  rownames(pars) <- 1:8
  dat <- make_data(pars, des, n_trials = 60)
  emc <- make_emc(dat, des, rt_resolution = 0.001, n_chains = 2)
  stage_r <- EMC2:::run_stage
  blocks <- differing <- 0L
  runs <- NULL
  both <- function(pmwgs, ...) {
    seed <- get(".Random.seed", envir = globalenv())
    options(emc2.native_steps = FALSE)
    a <- stage_r(pmwgs, ...)
    seed_a <- get(".Random.seed", envir = globalenv())
    assign(".Random.seed", seed, envir = globalenv())
    options(emc2.native_steps = TRUE)
    b <- stage_r(pmwgs, ...)
    blocks <<- blocks + 1L
    if (!identical(a$samples, b$samples) || !identical(a$rng, b$rng) ||
        !identical(a$shear$plan, b$shear$plan) ||
        !identical(a$collapsed$plan, b$collapsed$plan) ||
        !identical(seed_a, get(".Random.seed", envir = globalenv())))
      differing <<- differing + 1L
    b
  }
  local_mocked_bindings(run_stage = both, .package = "EMC2")
  fitted <- suppressWarnings(fit(emc, verbose = FALSE, cores_for_chains = 1,
                                 cores_per_chain = cores_per_chain,
                                 stop_criteria = list(preburn = list(iter = 20),
                                                      burn = list(iter = 60),
                                                      adapt = list(iter = 60),
                                                      sample = list(iter = 100))))
  options(emc2.native_steps = NULL)
  runs <- colSums(EMC2:::.emc_pool_state$native_runs)
  expect_gt(blocks, 3L)
  expect_identical(differing, 0L)
  # Every step ran natively in nearly every iteration of the native runs
  # (half of the iterations counted).
  steps <- c("group", "particle", "ensemble", "shear", "collapsed", "hmc", "store")
  expect_true(all(runs[steps] > 0))
  expect_gt(min(runs[steps]) / runs[["iter"]], 0.4)
}

test_that("the native iteration steps are the R steps on the same random stream", {
  skip_on_cran()
  native_steps_match(1L)
})

test_that("the native iteration steps are the R steps with a worker pool", {
  skip_on_cran()
  skip_on_os("windows")
  native_steps_match(2L)
})

# The workers' native claim of the subject queue (emc2_native_particle_claim)
# against their per-subject path: each stage block run both ways from the
# same state and random stream, on forked workers so that each run's workers
# take the option set for it.
test_that("workers claiming subjects in the particle kernel give the per-subject replies", {
  skip_on_cran()
  skip_on_os("windows")
  old <- options(emc2.fit_host = FALSE, emc2.worker_backend = "fork")
  on.exit(options(old), add = TRUE)
  des <- design(factors = list(subjects = 1:8, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(53)
  pars <- cbind(mu = log(0.4) + rnorm(8, 0, 0.15), sigma = log(0.05) + rnorm(8, 0, 0.3),
                tau = log(0.2) + rnorm(8, 0, 0.3))
  rownames(pars) <- 1:8
  dat <- make_data(pars, des, n_trials = 60)
  emc <- make_emc(dat, des, rt_resolution = 0.001, n_chains = 1)
  stage_r <- EMC2:::run_stage
  blocks <- differing <- 0L
  both <- function(pmwgs, ...) {
    seed <- get(".Random.seed", envir = globalenv())
    options(emc2.native_claim = FALSE)
    a <- stage_r(pmwgs, ...)
    seed_a <- get(".Random.seed", envir = globalenv())
    assign(".Random.seed", seed, envir = globalenv())
    options(emc2.native_claim = TRUE)
    b <- stage_r(pmwgs, ...)
    blocks <<- blocks + 1L
    if (!identical(a$samples, b$samples) || !identical(a$rng, b$rng) ||
        !identical(seed_a, get(".Random.seed", envir = globalenv())))
      differing <<- differing + 1L
    b
  }
  local_mocked_bindings(run_stage = both, .package = "EMC2")
  fitted <- suppressWarnings(fit(emc, verbose = FALSE, cores_for_chains = 1,
                                 cores_per_chain = 3,
                                 stop_criteria = list(preburn = list(iter = 20),
                                                      burn = list(iter = 40),
                                                      adapt = list(iter = 40),
                                                      sample = list(iter = 40))))
  options(emc2.native_claim = NULL)
  expect_gt(blocks, 3L)
  expect_identical(differing, 0L)
})
