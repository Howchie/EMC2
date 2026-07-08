# Regression tests for the group_design information-criterion path
# (group__IC_standard / standard_subj_ll), which drives compare() / IC() /
# DIC / BPIC. Before the fix this crashed ("invalid 'length' argument")
# because the standard_subj_ll call passed arguments in the wrong order, fetched
# the mapped p-vector "mu" instead of the stacked "beta" coefficients, and took
# the parameter dimension from theta_mu instead of alpha.

RNGkind("L'Ecuyer-CMRG")

test_that("standard_subj_ll matches the dmvnorm reference (intercept-only design)", {
  set.seed(1)
  p <- 3; n_subj <- 5; N <- 4
  par_names <- c("v", "a", "t0")
  # Intercept-only group designs: each subject's mean is just the intercept
  # vector, so M == p and standard_subj_ll must reduce to the plain MVN density.
  group_designs <- EMC2:::add_group_design(par_names, NULL, n_subj)

  theta_mu <- matrix(rnorm(p * N), p, N)
  alpha <- array(rnorm(p * n_subj * N), c(p, n_subj, N))
  theta_var <- array(0, c(p, p, N))
  for (i in seq_len(N)) {
    A <- matrix(rnorm(p * p), p, p)
    theta_var[, , i] <- crossprod(A) + diag(p)
  }

  lls <- EMC2:::standard_subj_ll(theta_var, theta_mu, alpha, n_subj, group_designs)

  ref <- vapply(seq_len(N), function(i) {
    sum(mvtnorm::dmvnorm(t(alpha[, , i]), theta_mu[, i], theta_var[, , i], log = TRUE))
  }, numeric(1))

  expect_equal(lls, ref, tolerance = 1e-10)
})

test_that("compare() returns finite DIC/BPIC for a group_design fit", {
  skip_on_cran()
  set.seed(123)

  subs <- unique(forstmann$subjects)[1:4]
  dat <- droplevels(forstmann[forstmann$subjects %in% subs, ])
  subj_design <- design(data = dat, model = LNR,
                        formula = list(m ~ 1, s ~ 1, t0 ~ 1))
  gd <- group_design(formula = list(m ~ 1), data = dat,
                     subject_design = subj_design)
  pri <- prior(subj_design, group_design = gd)
  emc <- make_emc(dat, subj_design, group_design = gd, prior_list = pri,
                  compress = FALSE, n_chains = 2)
  emc <- fit(emc, cores_for_chains = 1, stop_criteria = list(
    preburn = list(iter = 10), burn = list(iter = 10), adapt = list(iter = 10),
    sample = list(iter = 25)), verbose = FALSE, particle_factor = 20,
    step_size = 10)

  # BayesFactor = FALSE: the DIC/BPIC (IC) path is what commit 4 fixes; bridge
  # sampling needs far more samples than this short fit and is covered elsewhere.
  cmp <- compare(list(m = emc), BayesFactor = FALSE, print_summary = FALSE)
  expect_true(all(is.finite(cmp$DIC)))
  expect_true(all(is.finite(cmp$BPIC)))
})
