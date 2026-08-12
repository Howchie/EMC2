# Commit 7: bridge sampling already supports group designs (the handling lives in
# the variant functions, not bridge_sampling.R). These tests (a) pin the
# vectorized group-likelihood micro-cleanup in bridge_group_and_prior_and_jac_standard
# against the explicit per-subject loop it replaced, and (b) encode the
# equivalence property that an intercept-only group_design reproduces the
# no-group-design log marginal likelihood.

RNGkind("L'Ecuyer-CMRG")

test_that("bridge group indices follow separate covariance blocks", {
  p <- 5; n_subj <- 3; n_iter <- 4
  samples <- list(
    n_pars = p,
    n_subjects = n_subj,
    par_group = c(1, 1, 2, 2, 3),
    is_blocked = rep(TRUE, p),
    group_designs = NULL,
    samples = list(
      theta_mu = matrix(0, p, n_iter),
      a_half = matrix(1, p, n_iter),
      theta_var = array(rep(diag(p), n_iter), c(p, p, n_iter))
    )
  )
  idx <- rep(TRUE, n_iter)
  base <- matrix(0, n_iter, p * n_subj)
  complete <- EMC2:::bridge_add_group(base, samples, idx, "standard")
  info <- EMC2:::bridge_add_info_standard(list(), samples)

  expect_equal(max(info$group_idx), ncol(complete))
  expect_equal(length(info$group_idx), ncol(complete) - ncol(base))
})

test_that("vectorized group likelihood equals the per-subject loop", {
  set.seed(3)
  p <- 3; n_subj <- 6
  # covariate design on param 1 (M = 2 + 1 + 1 = 4 stacked coefficients)
  X1 <- cbind(1, as.numeric(scale(rnorm(n_subj))))
  group_designs <- list(m = X1,
                        s = matrix(1, n_subj, 1),
                        t0 = matrix(1, n_subj, 1))
  M <- sum(vapply(group_designs, ncol, integer(1)))
  theta_mu_row <- rnorm(M)
  A <- matrix(rnorm(p * p), p, p); var_curr <- crossprod(A) + diag(p)
  proposals_list <- lapply(seq_len(n_subj), function(s) matrix(rnorm(p), nrow = 1))
  i <- 1

  # OLD explicit loop (verbatim from the pre-cleanup implementation)
  group_ll_loop <- 0
  for (s in seq_len(n_subj)) {
    alpha_s <- proposals_list[[s]][i, ]
    mu_s <- numeric(p); par_idx <- 0
    for (k in seq_len(p)) {
      x_sk <- group_designs[[k]][s, , drop = FALSE]
      mu_s[k] <- x_sk %*% theta_mu_row[par_idx + 1:ncol(group_designs[[k]])]
      par_idx <- par_idx + ncol(group_designs[[k]])
    }
    group_ll_loop <- group_ll_loop +
      mvtnorm::dmvnorm(alpha_s, mu_s, var_curr, log = TRUE)
  }

  # NEW vectorized form
  subj_mu <- EMC2:::calculate_subject_means(group_designs, theta_mu_row)
  alpha_i <- vapply(proposals_list, function(pr) pr[i, ], numeric(p))
  group_ll_vec <- sum(mvtnorm::dmvnorm(t(alpha_i - subj_mu),
                                       sigma = var_curr, log = TRUE))

  expect_equal(group_ll_vec, group_ll_loop, tolerance = 1e-10)
})

test_that("bridge logMLL: intercept-only group_design ~ no group_design", {
  skip_on_cran()
  set.seed(2024)
  subs <- unique(forstmann$subjects)[1:5]
  dat <- droplevels(forstmann[forstmann$subjects %in% subs, ])
  sd <- design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))

  fit_it <- function(gd = NULL, pri = NULL) {
    emc <- make_emc(dat, sd, group_design = gd, prior_list = pri,
                    compress = FALSE, n_chains = 3)
    fit(emc, cores_for_chains = 1, stop_criteria = list(
      preburn = list(iter = 20), burn = list(max_gd = 1.2),
      adapt = list(min_unique = 40), sample = list(iter = 150)),
      verbose = FALSE, particle_factor = 30, step_size = 25)
  }

  emc0 <- fit_it()
  gd <- group_design(formula = list(m ~ 1), data = dat, subject_design = sd)
  pri <- prior(sd, group_design = gd)
  emcg <- fit_it(gd, pri)

  b0 <- run_bridge_sampling(emc0, cores_for_props = 1, both_splits = FALSE)
  bg <- run_bridge_sampling(emcg, cores_for_props = 1, both_splits = FALSE)

  expect_true(is.finite(b0) && is.finite(bg))
  # An intercept-only group design is the same model as the plain hierarchical
  # one, so the log marginal likelihoods must agree up to bridge MC error. This
  # pins the plumbing (that a group_design fit is bridge-able at all), not the
  # point estimate; a real plumbing break gives a large gap or an error.
  expect_equal(bg, b0, tolerance = 1.0)
})
