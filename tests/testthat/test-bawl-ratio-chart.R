skip_model_validation()

# BAwL(parameterization = "ratio"): the sampled leak coordinate is
# r = k / (mean launch strength) instead of the leak rate k itself, so the LBA
# limit sits at r = 0 rather than at a pole.  The map is a bijection, so every
# test below is an equivalence at k = r * mean -- the chart must not change the
# model, only the coordinates it is sampled in.
#
# The leak is rebuilt in two independent places that must agree exactly: the
# Ttransform of BAwL() (R/model_LBA.R), which feeds the R kernels and the
# simulators, and bawl_leak() (src/model_BAwL.cpp), which feeds the compiled
# likelihood.  The compiled path never runs Ttransform, so a drift between the
# two would be silent.

ratio_data <- function(n = 40, seed = 11) {
  set.seed(seed)
  data.frame(
    subjects = factor(rep(1, n)), trials = seq_len(n),
    S = factor(sample(c("r1", "r2"), n, TRUE), levels = c("r1", "r2")),
    R = factor(sample(c("r1", "r2"), n, TRUE), levels = c("r1", "r2")),
    rt = round(runif(n, 0.35, 1.4), 3))
}

ratio_design <- function(dat, model, formula, constants = NULL) {
  utils::capture.output(des <- suppressMessages(design(
    data = dat, model = model, matchfun = function(d) d$S == d$lR,
    formula = formula, constants = constants)))
  des
}

# The compiled likelihood is the authoritative path for these models, so the
# charts are compared through it rather than through calc_ll_R.
ratio_ll <- function(dat, des, vals) {
  p <- sampled_pars(des)
  p[] <- vals[names(p)]
  dadm <- design_model(dat, des)
  as.numeric(calc_ll_manager(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    dadm = dadm, model = des$model))
}

geom <- c(B = log(0.9), A = log(0.4), t0 = log(0.2))

# ---------------------------------------------------------------------------
# The column contract
# ---------------------------------------------------------------------------

test_that("the ratio chart swaps r for k in place and keeps k's bounds", {
  rate <- BAwL()
  ratio <- BAwL(parameterization = "ratio")

  # r takes k's POSITION: the leading p_types names are what
  # validate_col_prefix() checks against emc2col::bawl::spec_ratio().
  expect_identical(names(ratio$p_types)[1:6], c("v", "sv", "B", "A", "t0", "r"))
  expect_identical(names(rate$p_types)[1:6], c("v", "sv", "B", "A", "t0", "k"))
  expect_identical(ratio$p_types_canonical, c("v", "sv", "B", "A", "t0", "r"))
  expect_false("k" %in% names(ratio$p_types))

  # Same defaults, bounds and zero exception as the leak rate it replaces.
  expect_identical(unname(ratio$p_types[["r"]]), unname(rate$p_types[["k"]]))
  expect_identical(ratio$transform$func[["r"]], "exp")
  expect_identical(unname(ratio$bound$minmax[, "r"]),
                   unname(rate$bound$minmax[, "k"]))
  expect_identical(unname(ratio$bound$exception[["r"]]), 0)

  expect_identical(names(BAwL(drift_distribution = "lognormal",
                              parameterization = "ratio")$p_types)[1:6],
                   c("mean", "cv", "B", "A", "t0", "r"))
  expect_identical(names(BAwL(drift_distribution = "weibull",
                              parameterization = "ratio")$p_types)[1:6],
                   c("shape", "mean", "B", "A", "t0", "r"))
})

test_that("the c_name infix is RAT and never collides with the IO flag", {
  # resolve_race_model_adapter() detects the untruncated launch with an
  # unanchored find("IO"), and so do R/likelihood.R and R/make_data.R.  An
  # infix of "RATIO" would match all three and silently turn posdrift off.
  expect_identical(BAwL(parameterization = "ratio")$c_name, "BAwLRAT")
  expect_false(grepl("IO", BAwL(parameterization = "ratio")$c_name, fixed = TRUE))
  expect_identical(BAwL(posdrift = FALSE, parameterization = "ratio")$c_name,
                   "BAwLIORAT")
  expect_identical(BAwL(drift_distribution = "lognormal",
                        parameterization = "ratio")$c_name, "BAwLRAT_LOGN")
  expect_identical(BAwL(drift_distribution = "weibull", erlang_type = "local_kill",
                        parameterization = "ratio")$c_name,
                   "BAwLRAT_WEIB_LOCAL_KILL")
})

test_that("the ratio chart refuses the launches and paths it cannot support", {
  # mu is the log-scale median, not a natural-scale mean, so there is no column
  # for r to divide by.
  expect_error(BAwL(drift_distribution = "splitlognormal",
                    parameterization = "ratio"), "natural-scale mean")
  # The correlated kernels still read the leak straight out of k's column.
  expect_error(BAwL(correlated = TRUE, parameterization = "ratio"),
               "correlated = TRUE")
  expect_error(BAwLcorr(), NA)
})

# ---------------------------------------------------------------------------
# Equivalence at the map, through the compiled likelihood
# ---------------------------------------------------------------------------

test_that("the ratio chart reproduces the rate chart at k = r * mean", {
  dat <- ratio_data()

  # normal launch: the mean is v, the first kernel column.
  v <- 1.6; k <- 0.6
  d_k <- ratio_design(dat, function() BAwL(),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), c(sv = log(1)))
  d_r <- ratio_design(dat, function() BAwL(parameterization = "ratio"),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
  expect_equal(ratio_ll(dat, d_r, c(geom, v = v, r = log(k / v))),
               ratio_ll(dat, d_k, c(geom, v = v, k = log(k))),
               tolerance = 1e-12)

  # lognormal launch: the mean is the `mean` column, still first.
  mn <- 2.0; k2 <- 0.5
  bl <- c(geom, mean = log(mn), cv = log(0.4))
  d_k2 <- ratio_design(dat, function() BAwL(drift_distribution = "lognormal"),
                       list(mean ~ 1, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1))
  d_r2 <- ratio_design(dat, function() BAwL(drift_distribution = "lognormal",
                                            parameterization = "ratio"),
                       list(mean ~ 1, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1))
  expect_equal(ratio_ll(dat, d_r2, c(bl, r = log(k2 / mn))),
               ratio_ll(dat, d_k2, c(bl, k = log(k2))),
               tolerance = 1e-12)
})

test_that("the Weibull launch divides by its SECOND column, not its first", {
  # (shape, mean) is the one launch pair whose mean is not column 0.  Reading
  # the mean from column 0 would silently divide by the shape here and nowhere
  # else, so this case is what pins down bawl_leak()'s index.
  dat <- ratio_data()
  shape <- 2; mn <- 1.8; k <- 0.45
  bw <- c(geom, shape = log(shape), mean = log(mn))
  d_k <- ratio_design(dat, function() BAwL(drift_distribution = "weibull"),
                      list(shape ~ 1, mean ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1))
  d_r <- ratio_design(dat, function() BAwL(drift_distribution = "weibull",
                                           parameterization = "ratio"),
                      list(shape ~ 1, mean ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1))
  expect_equal(ratio_ll(dat, d_r, c(bw, r = log(k / mn))),
               ratio_ll(dat, d_k, c(bw, k = log(k))),
               tolerance = 1e-12)
  # Dividing by the shape instead would land somewhere else entirely.
  expect_false(isTRUE(all.equal(
    ratio_ll(dat, d_r, c(bw, r = log(k / shape))),
    ratio_ll(dat, d_k, c(bw, k = log(k))))))
})

test_that("the ratio chart matches the trend engine's independent k = r * v", {
  # make_trend(par_input) imposes the same coupling inside the C++ mapping and
  # shares no code with the chart, so it is an oracle rather than a mirror.
  dat <- ratio_data()
  v <- 1.6; k <- 0.6
  tr <- make_trend(par_names = "k", cov_names = NULL, par_input = list("v"),
                   kernels = "poly2", bases = "identity",
                   phase = "posttransform", at = NULL)
  utils::capture.output(d_tr <- suppressMessages(design(
    data = dat, model = function() BAwL(), matchfun = function(d) d$S == d$lR,
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k.d1 ~ 1),
    constants = c(sv = log(1), k = 0, k.d2 = 0),
    transform = list(func = c(k.d1 = "exp")), trend = tr)))
  d_r <- ratio_design(dat, function() BAwL(parameterization = "ratio"),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
  expect_equal(ratio_ll(dat, d_r, c(geom, v = v, r = log(k / v))),
               ratio_ll(dat, d_tr, c(geom, v = v, k.d1 = log(k / v))),
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# The limits: r = 0 and a non-positive mean drift
# ---------------------------------------------------------------------------

test_that("r = 0 is exactly the LBA limit", {
  dat <- ratio_data()
  d_r <- ratio_design(dat, function() BAwL(parameterization = "ratio"),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
  d_lba <- ratio_design(dat, function() LBA(),
                        list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1), c(sv = log(1)))
  expect_equal(ratio_ll(dat, d_r, c(geom, v = 1.6, r = log(0))),
               ratio_ll(dat, d_lba, c(geom, v = 1.6)), tolerance = 1e-12)

  # At the LBA limit a positive-drift model finishes with probability one, so
  # this also proves the "RAT" infix did not trip the IO/posdrift detection:
  # an untruncated launch leaves defective mass behind.
  pars <- cbind(v = 1.6, sv = 1, B = 0.9, A = 0.4, t0 = 0.2, r = 0,
                mG = 1, mK = 1, lambda_g = 0, lambda_k = 0)
  pars <- cbind(pars, k = pars[, "r"] * pmax(pars[, "v"], 0),
                b = pars[, "B"] + pars[, "A"])
  expect_equal(BAwL(parameterization = "ratio")$pfun(Inf, pars), 1,
               tolerance = 1e-10)
  expect_lt(BAwL(posdrift = FALSE, parameterization = "ratio")$pfun(Inf, pars), 1)
})

test_that("a non-positive mean drift clamps the leak to zero", {
  # k = r * max(v, 0).  The clamp is on the mean-drift PARAMETER and is
  # independent of posdrift, which governs the per-trial draws.
  dat <- ratio_data()
  for (pd in c(TRUE, FALSE)) {
    d_k <- ratio_design(dat, function() BAwL(posdrift = pd),
                        list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), c(sv = log(1)))
    d_r <- ratio_design(dat, function() BAwL(posdrift = pd,
                                             parameterization = "ratio"),
                        list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
    expect_equal(ratio_ll(dat, d_r, c(geom, v = -0.5, r = log(0.8))),
                 ratio_ll(dat, d_k, c(geom, v = -0.5, k = log(0))),
                 tolerance = 1e-12)
  }
})

test_that("a positive mean keeps its full leak under posdrift = FALSE", {
  # The clamp must not leak into the IO model's ordinary regime.
  dat <- ratio_data()
  v <- 1.6; k <- 0.6
  d_k <- ratio_design(dat, function() BAwL(posdrift = FALSE),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), c(sv = log(1)))
  d_r <- ratio_design(dat, function() BAwL(posdrift = FALSE,
                                           parameterization = "ratio"),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
  leaky <- ratio_ll(dat, d_r, c(geom, v = v, r = log(k / v)))
  expect_equal(leaky, ratio_ll(dat, d_k, c(geom, v = v, k = log(k))),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    leaky, ratio_ll(dat, d_k, c(geom, v = v, k = log(0))))))
})

# ---------------------------------------------------------------------------
# The simulator reads the leak Ttransform derives
# ---------------------------------------------------------------------------

test_that("the ratio chart simulates the same data as the rate chart", {
  dat <- ratio_data()
  v <- 1.6; k <- 0.6
  d_k <- ratio_design(dat, function() BAwL(),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), c(sv = log(1)))
  d_r <- ratio_design(dat, function() BAwL(parameterization = "ratio"),
                      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, r ~ 1), c(sv = log(1)))
  p_k <- sampled_pars(d_k); p_k[] <- c(geom, v = v, k = log(k))[names(p_k)]
  p_r <- sampled_pars(d_r); p_r[] <- c(geom, v = v, r = log(k / v))[names(p_r)]

  set.seed(99); sim_k <- make_data(p_k, design = d_k, n_trials = 60)
  set.seed(99); sim_r <- make_data(p_r, design = d_r, n_trials = 60)
  expect_identical(sim_r$R, sim_k$R)
  expect_identical(sim_r$rt, sim_k$rt)
})
