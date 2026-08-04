# Canary regression test for the build-flag bug fixed in the "restore IEEE NaN
# semantics" commit. -ffast-math implies -ffinite-math-only, which lets the
# compiler assume no NaN/Inf operands and fold `NaN == x` to TRUE. NA_REAL is a
# NaN, so under the buggy flags an NA-rt (both-sided censored) trial silently
# takes the `rt == R_NegInf` left-censor branch in c_log_likelihood_race and is
# scored as if it were only left-censored. -fno-finite-math-only restores correct
# Inf/NaN comparisons. If someone "simplifies" the Makevars and drops that flag,
# this test fails.

lI_single <- function(d) factor(rep(1, nrow(d)), levels = 1)

design_rexg <- design(
  factors = list(subjects = 1, S = 1),
  Rlevels = 1,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1),
  functions = list(lI = lI_single),
  model = REXG
)

p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18))
p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))

# Common censoring bounds. The observable window is [LC, UC]; a trial whose RT
# falls in the truncation-but-censored bands [LT, LC] or [UC, UT] is recorded
# with rt = -Inf (below LC), rt = +Inf (above UC), or rt = NA (censored, side
# unknown -> union of both bands).
LT <- 0.20; LC <- 0.35; UC <- 0.90; UT <- 1.50

single_trial <- function(rt_value) {
  data.frame(
    subjects = factor(1),
    S = factor(1),
    R = factor(1, levels = 1),
    rt = rt_value,
    LT = LT, LC = LC, UC = UC, UT = UT
  )
}

trial_ll <- function(rt_value) {
  dat <- single_trial(rt_value)
  emc <- make_emc(dat, design_rexg, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  as.numeric(EMC2:::calc_ll_oo(
    p_mat, dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = p_types,
    min_ll = log(1e-10),
    trend = model$trend
  ))
}

test_that("NA-rt (both-sided censored) differs from -Inf-rt (left-censored)", {
  ll_left  <- trial_ll(-Inf)  # mass over [LT, LC]
  ll_right <- trial_ll(Inf)   # mass over [UC, UT]
  ll_na    <- trial_ll(NA_real_)

  # All three intervals carry non-trivial, distinct probability mass.
  expect_true(is.finite(ll_left))
  expect_true(is.finite(ll_right))
  expect_true(is.finite(ll_na))

  # The NA-rt likelihood is the union of the two censoring bands: with the same
  # winner-specific integrator over [LT, LC] and [UC, UT], it must equal
  # log(exp(ll_left) + exp(ll_right)). This is the R reference formula.
  ll_union <- log(exp(ll_left) + exp(ll_right))
  expect_equal(ll_na, ll_union, tolerance = 1e-8)

  # Canary: under -ffinite-math-only the NA trial collapses onto the left-censor
  # branch and ll_na == ll_left. The right band adds strictly positive mass, so
  # the correct answer must be strictly larger.
  expect_gt(ll_na, ll_left + 1e-6)
})
