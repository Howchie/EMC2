# Regression test for group_design() subject-order hardening. The design
# matrices are built via factor-level ordering, but the "data" attribute
# (consumed by minimal_group_design and summary()) used to be appearance-ordered,
# so with unsorted input data the two diverged and covariates misaligned. Both
# are now keyed on levels(factor(subjects)).

RNGkind("L'Ecuyer-CMRG")

make_base <- function() {
  set.seed(42)
  base <- droplevels(forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:6], ])
  subs <- levels(factor(base$subjects))
  covs <- stats::setNames(round(rnorm(length(subs)), 3), subs)
  base$agecov <- covs[as.character(base$subjects)]
  list(base = base, covs = covs)
}

subj_design_lnr <- function(base) {
  design(data = base, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1))
}

test_that("group_design is invariant to input row order", {
  bc <- make_base(); base <- bc$base; covs <- bc$covs
  sd <- subj_design_lnr(base)
  mk <- function(d) group_design(formula = list(m ~ agecov), data = d,
                                 subject_design = sd)

  set.seed(7)
  gd_sorted <- mk(base[order(base$subjects), ])
  gd_shuf   <- mk(base[sample(nrow(base)), ])

  # (a) design matrices identical regardless of input order
  expect_equal(unname(gd_sorted[["m"]]), unname(gd_shuf[["m"]]))
  # rows are labelled and ordered by subject level
  expect_equal(rownames(gd_shuf[["m"]]), levels(factor(base$subjects)))

  # (b) the "data" attribute is level-ordered and identical across orderings
  expect_equal(attr(gd_sorted, "data"), attr(gd_shuf, "data"))

  # (c) each subject is paired with its OWN covariate value
  sd_shuf <- attr(gd_shuf, "data")
  expect_equal(sd_shuf[["agecov"]], unname(covs[rownames(sd_shuf)]))
})

test_that("a shuffled-input group_design is consumable by the sampler", {
  skip_on_cran()
  set.seed(11)
  bc <- make_base(); base <- bc$base
  shuffled <- base[sample(nrow(base)), ]
  sd <- subj_design_lnr(shuffled)
  gd <- group_design(formula = list(m ~ agecov), data = shuffled, subject_design = sd)
  pri <- prior(sd, group_design = gd)
  emc <- make_emc(shuffled, sd, group_design = gd, prior_list = pri,
                  compress = FALSE, n_chains = 1)
  # init_chains drives make_emc's minimal_group_design mapping and the sampler's
  # calculate_subject_means(group_designs, ...) start step -- both consumers of
  # the (now level-ordered) design matrices and "data" attribute.
  expect_no_error(init_chains(emc, particles = 10, cores_per_chain = 1))
})
