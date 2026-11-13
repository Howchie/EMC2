# Proposal refreshes use transition counts, independent of stored history.

local_rng_guard()

test_that("production refreshes its proposals on a quadrupling schedule", {
  expect_identical(EMC2:::.emc_next_sample_rebuild(0), 100L)
  expect_identical(EMC2:::.emc_next_sample_rebuild(100), 400L)
  expect_identical(EMC2:::.emc_next_sample_rebuild(400), 1600L)
  expect_identical(EMC2:::.emc_next_sample_rebuild(1600), 6400L)
})

test_that("production proposals are not rebuilt before their next refresh", {
  calls <- 0L
  local_mocked_bindings(
    create_chain_proposals = function(emc, ...) { calls <<- calls + 1L; emc },
    create_eff_proposals = function(emc, ...) { calls <<- calls + 1L; emc },
    .package = "EMC2")
  emc <- list(list(samples = list(stage = rep("adapt", 20), idx = 20L)))
  emc <- EMC2:::add_proposals(emc, "sample", 1, NULL)
  expect_equal(calls, 2L)
  expect_identical(emc[[1]]$sample_proposals_next, 100L)
  emc[[1]]$sample_transitions <- 50L
  emc[[1]]$samples <- list(stage = c(rep("adapt", 20), rep("sample", 50)), idx = 70L)
  later <- EMC2:::add_proposals(emc, "sample", 1, NULL)
  expect_equal(calls, 2L)
  expect_identical(later, emc)
  expect_false(EMC2:::.emc_sample_rebuild_due(emc))
})

test_that("production refreshes use every stored production draw", {
  calls <- list()
  local_mocked_bindings(
    create_chain_proposals = function(emc, samples_idx = NULL, ...) {
      calls[[length(calls) + 1L]] <<- list(kind = "chain", window = samples_idx)
      emc
    },
    create_eff_proposals = function(emc, n_cores, ...) {
      calls[[length(calls) + 1L]] <<- list(kind = "eff")
      emc
    },
    .package = "EMC2")
  make_chain <- function(n_adapt, n_sample) {
    st <- c(rep("adapt", n_adapt), rep("sample", n_sample))
    list(samples = list(stage = st, idx = length(st)))
  }
  # run_emc() advances sample_transitions by the transitions it ran; thinning
  # stores fewer draws than that.
  step <- function(emc, n_sample, n_stored = n_sample) {
    for (j in seq_along(emc)) {
      nxt <- emc[[j]]$sample_proposals_next
      emc[[j]]$samples <- make_chain(50, n_stored)$samples
      emc[[j]]$sample_proposals_next <- nxt
      if (!is.null(emc[[j]]$sample_transitions))
        emc[[j]]$sample_transitions <- n_sample
    }
    calls <<- list()
    EMC2:::add_proposals(emc, "sample", 1, NULL)
  }
  emc <- list(make_chain(50, 0), make_chain(50, 0))
  emc <- step(emc, 0)
  # Production start: warm-up window (NULL = the default recent window).
  expect_length(calls, 2L)
  expect_null(calls[[1]]$window)
  expect_identical(emc[[1]]$sample_proposals_next, 100L)
  emc <- step(emc, 50)
  expect_length(calls, 0L)                   # between rebuilds: unchanged
  emc <- step(emc, 100)
  expect_identical(calls[[1]]$window, 51:150) # every production draw so far
  expect_identical(emc[[2]]$sample_proposals_next, 400L)
  for (n in c(150, 200, 300)) { emc <- step(emc, n); expect_length(calls, 0L) }
  # Rebuild times follow transitions, not stored draws: thinned 4x, the stage
  # holds 100 draws after 400 transitions and still rebuilds on time.
  emc <- step(emc, 399, n_stored = 99)
  expect_length(calls, 0L)
  emc <- step(emc, 400, n_stored = 100)
  expect_identical(calls[[1]]$window, 51:150)
  expect_identical(emc[[1]]$sample_proposals_next, 1600L)
  # Warm-up stages still rebuild every time. Adapt tunes the full production
  # mixture, so it also builds the efficient proposal.
  calls <- list()
  EMC2:::add_proposals(list(make_chain(10, 0)), "adapt", 1, NULL)
  expect_identical(vapply(calls, `[[`, "", "kind"), c("chain", "eff"))
})

# End to end: a real fit on a weak-data Gaussian mock whose exact posterior is
# known, with production long enough for two rebuilds (100, 400).
test_that("rebuilding sample proposals keeps an exact Gibbs posterior", {
  skip_model_validation()
  skip_on_os("windows")
  data("forstmann", package = "EMC2")
  subject_names <- as.character(unique(forstmann$subjects)[1:4])
  dat <- droplevels(forstmann[forstmann$subjects %in% subject_names, ])
  y <- stats::setNames(c(-0.7, 0.1, 0.3, 1.2), subject_names)
  observation_var <- 25
  des <- design(data = dat, model = LNR,
                formula = list(m ~ 1, s ~ 1), constants = c(t0 = log(0.2)))
  fake_likelihood <- function(proposals, dadm, model, r_cores = 1L,
                              marginalise = NULL, ...) {
    subject <- as.character(unique(dadm$subjects)[1L])
    ll <- -0.5 * (proposals[, "m"] - y[[subject]])^2 / observation_var
    matrix(ll, ncol = 1L, dimnames = list(NULL, "mock_gaussian"))
  }
  batch_means <- function(x, n_batch) {
    n_each <- floor(length(x) / n_batch)
    colMeans(matrix(x[seq_len(n_each * n_batch)], nrow = n_each))
  }
  # Exact Gibbs reference for the group mean of m (diagonal-gamma hierarchy).
  reference <- function(prior) {
    set.seed(20260926L)
    n_ref <- 70000L; sigma2 <- 0.3
    mu <- prior$theta_mu_mean[1L]; out <- numeric(n_ref)
    pv <- prior$theta_mu_var[1L, 1L]
    for (it in seq_len(n_ref)) {
      av <- 1 / (1 / sigma2 + 1 / observation_var)
      alpha <- stats::rnorm(length(y), av * (mu / sigma2 + y / observation_var), sqrt(av))
      mv <- 1 / (1 / pv + length(y) / sigma2)
      mu <- stats::rnorm(1L, mv * (prior$theta_mu_mean[1L] / pv + sum(alpha) / sigma2), sqrt(mv))
      sigma2 <- 1 / stats::rgamma(1L, shape = prior$shape[1L] + length(y) / 2,
                                  rate = prior$rate[1L] + sum((alpha - mu)^2) / 2)
      out[it] <- mu
    }
    out[-seq_len(10000L)]
  }
  # A small ensemble pool keeps the test fast.
  small_pool <- EMC2:::.emc_ensemble_options()
  small_pool$pool <- 8L
  local_mocked_bindings(.emc_ensemble_options = function() small_pool, .package = "EMC2")
  set.seed(20260925L)
  fitted <- testthat::with_mocked_bindings({
    emc <- suppressMessages(make_emc(dat, des, type = "diagonal-gamma",
                                     compress = FALSE, n_chains = 2L))
    fit(emc, cores_for_chains = 1L,
        stop_criteria = list(preburn = list(iter = 20L), burn = list(iter = 50L),
                             adapt = list(iter = 50L, min_unique = 0L),
                             sample = list(iter = 1200L)),
        verbose = FALSE, particle_factor = 10, step_size = 25L,
        max_tries = 1L, trim = FALSE)
  }, calc_ll_manager = fake_likelihood, .package = "EMC2")
  fitted <- EMC2:::restore_duplicates(fitted)
  expect_identical(fitted[[1L]]$sample_proposals_next, 1600L)
  ref <- reference(fitted[[1L]]$prior)
  draws <- lapply(fitted, function(ch) ch$samples$theta_mu["m", ch$samples$stage == "sample"])
  se <- stats::sd(unlist(lapply(draws, batch_means, n_batch = 20L))) / sqrt(20L * length(draws))
  ref_se <- stats::sd(batch_means(ref, 50L)) / sqrt(50L)
  expect_lt(abs(mean(unlist(draws)) - mean(ref)) / sqrt(se^2 + ref_se^2), 5,
            label = "group mean")
  m2 <- unlist(lapply(draws, function(x) (x - mean(ref))^2))
  se2 <- stats::sd(unlist(lapply(draws, function(x) batch_means((x - mean(ref))^2, 20L)))) /
    sqrt(20L * length(draws))
  ref_se2 <- stats::sd(batch_means((ref - mean(ref))^2, 50L)) / sqrt(50L)
  expect_lt(abs(mean(m2) - mean((ref - mean(ref))^2)) / sqrt(se2^2 + ref_se2^2), 5,
            label = "group variance")
})
