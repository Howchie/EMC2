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
    model = plain_link(BAwL),
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

test_that("the default spread is set per design cell, not per coefficient", {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design_lba <- design(
    factors = list(subjects = 1:5, S = c("left", "right"), L = c("Low", "Med", "High")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = BAwL,
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ L, t0 ~ 1, A ~ 1),
    report_p_vector = FALSE
  )
  means <- sampled_pars(design_lba)
  means[] <- c(2, 2.5, log(.7), log(1.2), log(1.5), log(.2), log(.5))
  cv <- 0.3
  set.seed(7)
  eff <- make_random_effects(design_lba, means, n_subj = 4000, variance_proportion = cv)

  # Every B cell (log scale) has the log-normal spread of coefficient of
  # variation cv, whether it is the intercept or intercept + contrast.
  B_cells <- cbind(eff[, "B"], eff[, "B"] + eff[, "B_LMed"], eff[, "B"] + eff[, "B_LHigh"])
  expect_equal(apply(B_cells, 2, var), rep(log1p(cv^2), 3), tolerance = 0.08)
  expect_lt(cor(eff[, "B"], eff[, "B_LMed"]), -0.5)
  # Identity cells: sd = cv * |cell mean|.
  v_cells <- cbind(eff[, "v"], eff[, "v"] + eff[, "v_lMTRUE"])
  expect_equal(apply(v_cells, 2, sd), cv * c(2, 4.5), tolerance = 0.08)
  # Types stay independent.
  expect_lt(abs(cor(eff[, "B"], eff[, "v"])), 0.06)
})

test_that("make_random_effects applies the likelihood's strict bounds", {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design_lba <- design(
    factors = list(subjects = 1:2, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = plain_link(BAwL),
    constants = c(sv = log(1)),
    formula = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1, pContaminant ~ 1),
    report_p_vector = FALSE
  )
  means <- sampled_pars(design_lba)
  # pContaminant exactly on its upper bound (0.999): accepted by an inclusive
  # check, but the likelihood rejects it.
  means[] <- c(2, 2.5, log(.7), log(.2), log(.5), qnorm(0.999))
  expect_identical(pnorm(qnorm(0.999)), 0.999)
  zero_covariance <- matrix(0, nrow = length(means), ncol = length(means))
  expect_error(
    make_random_effects(design_lba, means, n_subj = 2,
                        covariances = zero_covariance, max_tries = 2),
    "Could not draw in-bound random effects"
  )
  means["pContaminant"] <- qnorm(0.05)
  expect_equal(dim(make_random_effects(design_lba, means, n_subj = 2,
                                       covariances = zero_covariance, max_tries = 2)),
               c(2, length(means)))
})
