local_rng_guard()

test_that("diagonal and blocked fits start and sample priors without cross-group covariance", {
  set.seed(4)
  design_in <- get_design(samples_LNR)[[1]]
  data_in <- get_data(samples_LNR)
  pars <- names(sampled_pars(design_in))
  off_block <- function(groups) outer(groups, groups, "!=")
  for (type in c("diagonal", "blocked")) {
    groups <- if (type == "diagonal") seq_along(pars) else
      c(1L, 1L, rep(2L, length(pars) - 2L))
    emc <- suppressMessages(make_emc(
      data_in, design_in, n_chains = 2, type = type,
      par_groups = if (type == "blocked") pars[1:2]))
    emc <- suppressMessages(init_chains(emc, particles = 10,
                                        cores_for_chains = 1))
    for (chain in emc) {
      start <- chain$samples$theta_var[, , 1L]
      expect_true(all(start[off_block(groups)] == 0))
      expect_true(all(eigen(start, only.values = TRUE)$values > 0))
    }
    # The prior overlays in plots and savage_dickey() sample it this way.
    draws <- EMC2:::get_objects(
      design = get_design(emc), type = emc[[1]]$type, sample_prior = TRUE,
      prior = get_prior(emc), selection = "Sigma", N = 50,
      sampler = emc)[[1]]$samples$theta_var
    expect_true(all(apply(draws, 3L, function(v) all(v[off_block(groups)] == 0))))
  }
})

test_that("a diagonal type stays diagonal when a prior_list is supplied", {
  design_in <- get_design(samples_LNR)[[1]]
  data_in <- get_data(samples_LNR)
  pars <- names(sampled_pars(design_in))
  standard_prior <- prior(design_in, type = "standard")
  diagonal_prior <- prior(design_in, type = "diagonal")
  expect_identical(attr(diagonal_prior, "sampler_type"), "diagonal")

  emc <- suppressMessages(make_emc(
    data_in, design_in, n_chains = 1, type = "diagonal",
    prior_list = standard_prior))
  expect_identical(emc[[1]]$type, "standard")
  expect_identical(emc[[1]]$par_group, seq_along(pars))
  expect_false(any(emc[[1]]$is_blocked))

  emc <- suppressMessages(init_chains(emc, particles = 10,
                                      cores_for_chains = 1))
  start <- emc[[1]]$samples$theta_var[, , 1L]
  expect_true(all(start[row(start) != col(start)] == 0))

  inferred <- suppressMessages(make_emc(
    data_in, design_in, n_chains = 1, prior_list = diagonal_prior))[[1]]
  expect_identical(inferred$par_group, seq_along(pars))
  expect_false(any(inferred$is_blocked))
})
