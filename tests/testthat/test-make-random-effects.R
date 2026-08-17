test_that("make_random_effects uses a bounded, transform-aware default spread", {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design_lba <- design(
    factors = list(
      subjects = 1:20,
      S = c("left", "right"),
      L = c("Low", "Med", "High")
    ),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL,
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ L, t0 ~ 1, A ~ 1),
    report_p_vector = FALSE
  )
  means <- sampled_pars(design_lba)
  means[] <- c(2, 2.5, log(.7), log(1.2), log(1.5), log(.2), log(.5))

  set.seed(42)
  effects <- make_random_effects(design_lba, means, n_subj = 20)
  simulated <- make_data(effects, design_lba, n_trials = 1)

  expect_equal(dim(effects), c(20, 7))
  expect_false(anyNA(simulated$R))
  expect_true(all(is.finite(simulated$rt)))
})

test_that("make_random_effects successfully redraws when a subset of subjects falls out of bounds", {
  matchfun <- function(d) as.numeric(factor(d$S)) == as.numeric(factor(d$lR))
  design_lba <- design(
    factors = list(
      subjects = 1:10,
      S = c("left", "right"),
      L = c("Low", "Med", "High")
    ),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL,
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ L, t0 ~ 1, A ~ 1),
    report_p_vector = FALSE
  )
  means <- sampled_pars(design_lba)
  # t0 close to bound (0.05) with moderate variance so some fail on attempt 1
  means[] <- c(2, 2.5, log(.7), log(1.2), log(1.5), log(.06), log(.5))

  for (seed in 1:5) {
    set.seed(seed)
    effects <- make_random_effects(design_lba, means, n_subj = 10, variance_proportion = 0.5)
    expect_equal(dim(effects), c(10, 7))
    expect_true(all(is.finite(effects)))
  }
})

test_that("make_random_effects properly aligns permuted named group_means and covariances", {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design_lba <- design(
    factors = list(subjects = 1:5, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL,
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1),
    report_p_vector = FALSE
  )
  canonical_names <- names(sampled_pars(design_lba))
  canonical_means <- c(2, 2.5, log(.7), log(.2), log(.5))
  names(canonical_means) <- canonical_names

  # Scramble the names
  rev_means <- rev(canonical_means)
  cov_mat <- diag(0.01, length(canonical_names))
  dimnames(cov_mat) <- list(canonical_names, canonical_names)
  rev_cov <- cov_mat[names(rev_means), names(rev_means)]

  set.seed(123)
  eff1 <- make_random_effects(design_lba, canonical_means, n_subj = 5, covariances = cov_mat)
  set.seed(123)
  eff2 <- make_random_effects(design_lba, rev_means, n_subj = 5, covariances = rev_cov)

  expect_equal(eff1, eff2)
})

test_that("make_random_effects errors when bounded draws cannot be resampled", {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design_lba <- design(
    factors = list(subjects = 1:2, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL,
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1),
    report_p_vector = FALSE
  )
  means <- sampled_pars(design_lba)
  means[] <- c(2, 2.5, log(.7), log(.01), log(.5))
  zero_covariance <- matrix(0, nrow = length(means), ncol = length(means))

  expect_error(
    make_random_effects(
      design_lba, means, n_subj = 2,
      covariances = zero_covariance, max_tries = 2
    ),
    "Could not draw in-bound random effects"
  )
})
