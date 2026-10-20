library(EMC2)

# Reference values are the drift integral of the Wald survivor over the
# positive-truncated normal drift law, computed with mpmath at 30 digits
# (mu = 0.4, b = 1.25, A = 0, s = 0.446).
test_that("RDMSWTN sv > 0 CDF and survivor stay accurate at huge sv / s", {
  ref <- data.frame(
    t  = c(0.15, 1.5, 0.15, 1.5, 0.15, 1.5),
    sv = c(1e6, 1e6, 1e8, 1e8, 1e12, 1e12),
    logS = c(-11.930632535655933, -14.319818452621854, -16.535802405671372,
             -18.924988322647653, -25.746142774456334, -28.135328691432615),
    logF = c(-6.5855727875363817e-06, -6.0392359179768584e-07,
             -6.5855534004014481e-08, -6.0392360207654340e-09,
             -6.5855532045915198e-12, -6.0392360218035972e-13))
  for (i in seq_len(nrow(ref))) {
    lF <- EMC2:::prdmswtn(ref$t[i], 0.4, 1.25, 0, 0.446, 0, ref$sv[i], log_out = TRUE)
    expect_equal(lF, ref$logF[i], tolerance = 1e-8)
    # A negligible guess clock makes prdmswtn return log(1 - F) with the
    # survivor resolved below the double resolution of F.
    lG <- EMC2:::prdmswtn(ref$t[i], 0.4, 1.25, 0, 0.446, 0, ref$sv[i], lambda_g = 1e-9,
                          guess = TRUE, log_out = TRUE)
    expect_equal(log(-expm1(lG)) + 1e-9 * ref$t[i], ref$logS[i], tolerance = 1e-8)
  }
})

xsv_ll <- function(sv, A, LT = 0.25, UT = 3.5, rt = c(1, 2), R = c("left", "right")) {
  dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                    R = factor(R, levels = c("left", "right")), rt = rt, LT = LT, UT = UT)
  des <- design(data = dat, Rlevels = c("left", "right"),
                matchfun = function(d) d$S == d$lR,
                formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ lM),
                model = RDMSWTN, report_p_vector = FALSE,
                constants = if (A == 0) c(A = log(0)) else NULL)
  p <- sampled_pars(des)
  p[] <- 0
  p[c("v", "v_lMTRUE", "t0", "s")] <- log(c(1.5, 2, 0.15, 0.7))
  if (A > 0) p["A"] <- log(A)
  p["sv"] <- log(sv[1]); p["sv_lMTRUE"] <- log(sv[2] / sv[1])
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = NULL)
  pm <- matrix(p, 1, dimnames = list(NULL, names(p)))
  dadm <- emc[[1]]$data[[1]]
  ll <- as.numeric(EMC2:::calc_ll_pw(pm, dadm = dadm, model = emc[[1]]$model))
  # Exact truncated value from the scalar kernels: f_w(t) S_l(t) / Z.
  mp <- EMC2:::get_pars_matrix_oo(pm, dadm, emc[[1]]$model)
  lS <- function(t, i) with(as.list(mp[i, ]), log(-expm1(
    EMC2:::prdmswtn(t, v, B + A, A, s, t0, sv, log_out = TRUE))))
  lf <- function(t, i) with(as.list(mp[i, ]),
    EMC2:::drdmswtn(t, v, B + A, A, s, t0, sv, log_out = TRUE))
  logZ <- log(exp(lS(LT, 1) + lS(LT, 2)) - exp(lS(UT, 1) + lS(UT, 2)))
  exact <- vapply(seq_along(rt), function(j) {
    rows <- 2 * j - 1:0
    w <- rows[dadm$winner[rows]]
    lf(rt[j], w) + lS(rt[j], setdiff(rows, w))
  }, numeric(1)) - logZ
  list(ll = ll, exact = unname(exact))
}

test_that("truncated RDMSWTN normalises a near-deterministic drift (sv below its cutoff)", {
  for (A in c(0, 0.3)) {
    r <- xsv_ll(c(1e-12, 1e-12), A)
    expect_equal(r$ll, r$exact, tolerance = 1e-6)
  }
})

test_that("truncated RDMSWTN floors after normalising, not per accumulator", {
  # The mismatch accumulator's huge drift spread leaves it ~1e-10 survival at
  # LT, so each trial has a component below min_ll while f / Z is well above it.
  r <- xsv_ll(c(6e10, 0.3), 0)
  expect_true(all(r$exact > -20))
  expect_equal(r$ll, r$exact, tolerance = 1e-6)
})

test_that("a narrow truncation window integrates the density instead of differencing survivors", {
  narrow_ll <- function(w) {
    dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                      R = factor("left", levels = c("left", "right")), rt = 1 + w / 2,
                      LT = 1, UT = 1 + w)
    des <- design(data = dat, Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1), model = RDMSWTN,
                  report_p_vector = FALSE, constants = c(A = log(0)))
    p <- sampled_pars(des)
    p[] <- 0
    p[c("v", "v_lMTRUE", "t0", "s", "sv")] <- log(c(1.5, 2, 0.15, 0.7, 0.3))
    emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                    verbose = FALSE, rt_resolution = NULL)
    ll <- as.numeric(EMC2:::calc_ll_pw(matrix(p, 1, dimnames = list(NULL, names(p))),
                                       dadm = emc[[1]]$data[[1]], model = emc[[1]]$model))
    # ll + log(width) tends to log P(left wins | rt) at the window as it narrows.
    ll + log(dat$UT - dat$LT)
  }
  expect_equal(narrow_ll(1e-14), narrow_ll(1e-6), tolerance = 1e-5)
})
