# refit_alpha(): subject-level parameters for new data with the group level
# frozen at the training posterior, and the calibrated score built on it.
#
# The thing that has to be true throughout is that the group level does not
# move: a "single" sampler cannot change it, but it can still be installed
# wrongly (in the wrong parameter order, or before make_emc overwrites it), and
# the alphas would then be shrunk towards the wrong place without any error.

# Two subjects, trimmed, so the refits in this file stay cheap.
held_out_pair <- function(n = 120) {
  dat  <- get_data(samples_LNR)
  subs <- levels(factor(dat$subjects))[1:2]
  out <- do.call(rbind, lapply(subs, function(s) utils::head(dat[dat$subjects == s, ], n)))
  out$subjects <- factor(as.character(out$subjects))
  out
}

refit_fast <- function(data, ...) {
  refit_alpha(samples_LNR, data, rt_resolution = NULL, n_chains = 2, iter = 50,
              verbose = FALSE, ...)
}

test_that("the frozen prior is the group level's predictive moments", {
  e <- samples_LNR
  mu    <- get_pars(e, selection = "mu", merge_chains = TRUE, return_mcmc = FALSE)
  Sigma <- get_pars(e, selection = "Sigma", merge_chains = TRUE, return_mcmc = FALSE)

  fp <- EMC2:::.frozen_group_prior(e, group_uncertainty = TRUE)
  expect_equal(fp$theta_mu_mean, rowMeans(mu))
  # law of total covariance: E[Sigma] + Var[mu]
  expect_equal(unname(fp$theta_mu_var),
               unname(apply(Sigma, c(1, 2), mean) + var(t(mu))))

  plug <- EMC2:::.frozen_group_prior(e, group_uncertainty = FALSE)
  expect_equal(unname(plug$theta_mu_var), unname(apply(Sigma, c(1, 2), mean)))
  # ignoring group-level uncertainty can only narrow the prior
  expect_true(all(diag(plug$theta_mu_var) <= diag(fp$theta_mu_var)))
})

test_that("refit_alpha installs the frozen group level and never moves it", {
  d  <- held_out_pair()
  rf <- refit_fast(d)
  fp <- EMC2:::.frozen_group_prior(samples_LNR)

  expect_s3_class(rf, "emc")
  expect_equal(rf[[1]]$type, "single")
  expect_equal(rf[[1]]$par_names, samples_LNR[[1]]$par_names)
  expect_equal(rf[[1]]$n_subjects, 2L)

  # fit() strips information duplicated across chains, so ask for it back: the
  # prior has to be right on every chain, not just the one that keeps a copy.
  for (chain in EMC2:::restore_duplicates(rf)) {
    # same order as the sampler's parameter vector, not the group level's
    expect_equal(chain$prior$theta_mu_mean, fp$theta_mu_mean[chain$par_names])
    expect_equal(dimnames(chain$prior$theta_mu_var)[[1]], chain$par_names)
    expect_true(attr(chain$prior, "frozen_group"))
    # a single sampler stores no group level at all, so there is nothing that
    # the new data could have updated
    expect_null(chain$samples$theta_mu)
  }

  # Every chain sampled against the same frozen prior, so they must agree; a
  # chain left on make_emc's default N(0, 1) prior would not.
  ch_means <- lapply(1:2, function(i) {
    a <- get_pars(rf, selection = "alpha", chain = i, by_subject = TRUE,
                  merge_chains = TRUE)
    vapply(a, function(x) colMeans(do.call(rbind, x)),
           numeric(length(rf[[1]]$par_names)))
  })
  expect_lt(max(abs(ch_means[[1]] - ch_means[[2]])), 0.3)
})

test_that("refitted alphas track the subjects' own data, not just the prior", {
  d  <- held_out_pair(200)
  rf <- refit_fast(d)
  a  <- get_pars(rf, selection = "alpha", by_subject = TRUE, merge_chains = TRUE)
  means <- t(vapply(a, function(x) colMeans(do.call(rbind, x)),
                    numeric(length(rf[[1]]$par_names))))

  # The two subjects differ in the training fit; refitting from their own data
  # with a shared prior must keep them apart rather than collapsing them onto
  # the frozen mean.
  train <- get_pars(samples_LNR, selection = "alpha", by_subject = TRUE,
                    merge_chains = TRUE)
  tm <- t(vapply(train[rownames(means)], function(x) colMeans(do.call(rbind, x)),
                 numeric(ncol(means))))
  expect_gt(cor(as.vector(means), as.vector(tm)), 0.9)
  expect_gt(max(abs(means[1, ] - means[2, ])), 0.1)
})

test_that("alpha_group_distance measures distance from the frozen group level", {
  d  <- held_out_pair()
  rf <- refit_fast(d)
  gd <- alpha_group_distance(rf)

  expect_equal(nrow(gd), 2L)
  expect_setequal(gd$subjects, levels(d$subjects))
  expect_true(all(c("d", "p", rf[[1]]$par_names) %in% names(gd)))
  expect_true(all(gd$p >= 0 & gd$p <= 1))

  # d is the Mahalanobis distance of the posterior mean under the frozen prior
  a <- get_pars(rf, selection = "alpha", by_subject = TRUE, merge_chains = TRUE)
  means <- t(vapply(a, function(x) colMeans(do.call(rbind, x)),
                    numeric(length(rf[[1]]$par_names))))
  ref <- sqrt(stats::mahalanobis(means[gd$subjects, , drop = FALSE],
                                 rf[[1]]$prior$theta_mu_mean,
                                 rf[[1]]$prior$theta_mu_var))
  expect_equal(gd$d, unname(ref))

  # an ordinary fit has no frozen group level to be distant from
  expect_error(alpha_group_distance(samples_LNR), "frozen group level")
})

test_that("the calibration split is disjoint, covering, and per subject", {
  d <- held_out_pair(100)
  sp <- EMC2:::.calibration_split(d, 0.5, "first")
  expect_equal(nrow(sp$calibration) + nrow(sp$scored), nrow(d))
  expect_length(intersect(rownames(sp$calibration), rownames(sp$scored)), 0)
  expect_equal(as.vector(table(sp$calibration$subjects)), c(50, 50))
  expect_equal(as.vector(table(sp$scored$subjects)), c(50, 50))
  # "first" means the subject's first trials, not the data frame's first rows
  expect_equal(rownames(sp$calibration),
               unlist(lapply(split(rownames(d), d$subjects), utils::head, 50),
                      use.names = FALSE))

  # a count rather than a proportion
  sp_n <- EMC2:::.calibration_split(d, 30, "first")
  expect_equal(as.vector(table(sp_n$calibration$subjects)), c(30, 30))

  set.seed(1)
  sp_r <- EMC2:::.calibration_split(d, 0.5, "random")
  expect_equal(as.vector(table(sp_r$calibration$subjects)), c(50, 50))
  expect_false(identical(rownames(sp_r$calibration), rownames(sp$calibration)))

  expect_error(EMC2:::.calibration_split(d, 500, "first"), "nothing.*to score")
  expect_error(EMC2:::.calibration_split(d, 0, "first"), "positive number")
})

test_that("the calibrated score is conditional on refitted alphas", {
  d <- held_out_pair(150)
  set.seed(3)
  sc <- score_held_out(samples_LNR, d, type = "calibrated", rt_resolution = NULL,
                       n_calibration = 0.5,
                       refit = list(n_chains = 2, iter = 50))

  expect_s3_class(sc, "emc.score")
  expect_equal(sc$type, "calibrated")
  expect_equal(sc$n_points, nrow(d) / 2)          # only the scored half
  expect_true(all(is.finite(sc$pointwise)))
  expect_equal(sc$n_subjects, 2L)
  expect_s3_class(sc$refit, "emc")
  expect_equal(as.vector(sc$n_calibration), c(75L, 75L))

  # It must be the refit's alphas doing the work: scoring the same rows with the
  # training fit's alphas is a different number.
  sp <- EMC2:::.calibration_split(d, 0.5, "first")
  cond <- score_held_out(samples_LNR, sp$scored, type = "trial", rt_resolution = NULL)
  expect_equal(sc$point_id, cond$point_id)
  expect_false(isTRUE(all.equal(sc$elpd, cond$elpd)))
  # ... but not a wildly different one, since both alphas describe the same person
  expect_lt(abs(sc$elpd - cond$elpd) / sc$n_points, 0.5)

  # pooling calibrated folds works like any other trial-level score, but a
  # calibrated point and a conditional one are not the same quantity
  cv <- combine_scores(list(sc))
  expect_equal(cv$elpd, sc$elpd)
  expect_error(combine_scores(list(sc, cond)), "different types")
  expect_error(combine_scores(list(sc, sc)), "overlap")
})
