# Regression test for the pContaminant compressed-path fix. A contaminant is an
# intrinsic omission that only manifests as a never-responded (rt == +Inf)
# trial, so only +Inf trials pick up the log(pC) mass. The expand branch of
# c_log_likelihood_race already used `rt == R_PosInf`; the compressed branch used
# `!R_FINITE(rt)`, which wrongly also credited left-censored (-Inf) and missing
# (NA) trials. Both branches now share one lambda. On a dataset whose trials are
# all unique (identity expand) the two branches must produce the same total ll.

lI_single <- function(d) factor(rep(1, nrow(d)), levels = 1)

design_rexg_pc <- design(
  factors = list(subjects = 1, S = 1),
  Rlevels = 1,
  formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1, pContaminant ~ 1),
  functions = list(lI = lI_single),
  model = REXG
)

# pContaminant is pnorm-transformed: qnorm(0.05) -> pC = 0.05.
p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0),
       pContaminant = qnorm(0.05))
p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))

# A handful of distinct finite RTs plus one left-censored (-Inf) trial, so the
# dadm rows are all unique and the expand vector is the identity map.
dat <- data.frame(
  subjects = factor(1),
  S = factor(1),
  R = factor(1, levels = 1),
  rt = c(0.45, 0.55, 0.70, 0.95, -Inf),
  LT = 0.20, LC = 0.35, UC = Inf, UT = Inf
)

emc <- make_emc(dat, design_rexg_pc, type = "single", n_chains = 1,
                compress = FALSE, verbose = FALSE)
dadm <- emc[[1]]$data[[1]]
model <- emc[[1]]$model()
p_types <- names(model$p_types)

designs <- list()
for (nm in p_types) {
  designs[[nm]] <- attr(dadm, "designs")[[nm]][
    attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
}

run_ll <- function(d) {
  constants <- attr(d, "constants")
  if (is.null(constants)) constants <- NA
  as.numeric(EMC2:::calc_ll_oo(
    p_mat, d,
    constants = constants,
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

test_that("pContaminant expand and compressed branches agree on a -Inf trial", {
  # Sanity: the expand attribute is present and identity (all trials unique).
  expect_true(length(attr(dadm, "expand")) > 0)

  ll_expand <- run_ll(dadm)

  # Strip the expand attribute to force the compressed summation branch on the
  # exact same underlying data.
  dadm_stripped <- dadm
  attr(dadm_stripped, "expand") <- integer(0)
  ll_compressed <- run_ll(dadm_stripped)

  expect_true(is.finite(ll_expand))
  expect_equal(ll_compressed, ll_expand, tolerance = 1e-10)
})
