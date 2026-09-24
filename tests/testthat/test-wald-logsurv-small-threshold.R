local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# Regression: a Wald accumulator whose threshold sits far below its diffusion
# scale, i.e. b/s small.
#
# There the accumulator has all but certainly crossed already, so its log
# survivor is large and negative.  wald_pt_log_surv() reaches that through the
# Mills-ratio difference R(-z1) - R(-z2); the difference cancels, and the
# cancellation fallback used to be the LARGE-argument asymptotic
# -R'(x) ~ 1/(x1 x2).  For b/s small both arguments approach ZERO, where
# -R'(x) -> 1 instead, so the fallback overstated the survivor by ~1/(x1 x2) and
# the fmin(., 0) clamp turned that into log S = 0: an accumulator certain to
# have finished was scored as certain to still be running.
#
# In a race that silently drops the accumulator from the likelihood, while the
# truncation normaliser reaches the same survivor by a natural-scale 1 - F and
# keeps it.  The two then disagree by tens of nats and the truncated likelihood
# log f - log Z runs away POSITIVE and unbounded, which is what this pins.

small_threshold_design <- function() {
  design(
    factors  = list(S = c("left", "right"), subjects = 1),
    Rlevels  = c("left", "right"),
    matchfun = function(d) d$S == d$lR,
    formula  = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 0 + lM, sv ~ 1),
    constants = c(A = log(0), sv = log(0)),
    model    = RDMSWTN,
    report_p_vector = FALSE
  )
}

small_threshold_pars <- function(des) {
  p <- sampled_pars(des)
  p["v"] <- log(3)
  p["B"] <- log(1)
  p["t0"] <- log(0.05)
  p["s_lMTRUE"] <- log(1)
  p["s_lMFALSE"] <- log(1)
  p
}

# Race log-likelihood of one trial from the natural-scale kernels: the winner's
# density times every loser's survivor.  Kept on the natural scale on purpose --
# it is the independent check on the log-space branch, and it stays accurate
# while 1 - F is still resolvable in double precision (s up to ~1e6 here).
race_ll_reference <- function(dadm, p) {
  is_match <- as.logical(as.character(dadm$lM))
  s_row <- ifelse(is_match, exp(p[["s_lMTRUE"]]), exp(p[["s_lMFALSE"]]))
  pars <- cbind(v = exp(p[["v"]]), b = exp(p[["B"]]), B = exp(p[["B"]]), A = 0,
                t0 = exp(p[["t0"]]), s = as.numeric(s_row), sv = 0,
                mG = 1, mK = 1, lambda_g = 0, lambda_k = 0)
  dens <- EMC2:::dRDMSWTN(dadm$rt, pars)
  surv <- 1 - EMC2:::pRDMSWTN(dadm$rt, pars)
  trial <- rep(seq_len(nrow(dadm) / 2L), each = 2L)
  vapply(split(seq_len(nrow(dadm)), trial), function(rows) {
    w <- rows[dadm$winner[rows]]
    l <- rows[!dadm$winner[rows]]
    log(dens[w]) + sum(log(surv[l]))
  }, numeric(1))
}

make_small_threshold_dadm <- function(des, p, truncate) {
  # Pin the stream: other files in the suite set RNGkind(), and the plateau
  # below only holds while no trial's loser survivor has fallen through the
  # min_ll floor, which depends on the simulated RTs.
  old <- RNGkind("Mersenne-Twister")
  on.exit(RNGkind(old[1]), add = TRUE)
  set.seed(1)
  dat <- make_data(p, design = des, n_trials = 40)
  dat <- dat[is.finite(dat$rt) & dat$rt > 0.25 & dat$rt < 4, ]
  if (truncate) {
    dat$LT <- 0.25
    dat$UT <- 4
  }
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dadm = emc[[1]]$data[[1]], model = emc[[1]]$model)
}

test_that("a saturated loser keeps its survivor in the race likelihood", {
  des <- small_threshold_design()
  p <- small_threshold_pars(des)
  fit <- make_small_threshold_dadm(des, p, truncate = FALSE)

  # s = 1e9 puts b/s at 1e-9.  That is past EMC2_CDF_SAT_MARGIN, so the natural
  # CDF is rejected and the log-space branch -- the broken one -- has to answer,
  # while 1 - F still carries ~7 digits and so remains a sound reference.
  for (s_false in c(1, 1e2, 1e4, 1e6, 1e9)) {
    q <- p
    q["s_lMFALSE"] <- log(s_false)
    got <- as.numeric(EMC2:::calc_ll_pw(
      matrix(q, nrow = 1, dimnames = list(NULL, names(q))),
      dadm = fit$dadm, model = fit$model))
    want <- unname(race_ll_reference(fit$dadm, q))
    keep <- want > log(1e-10)   # the likelihood floors below this
    expect_gt(sum(keep), 0)
    expect_equal(got[keep], want[keep], tolerance = 1e-5,
                 info = paste("s_lMFALSE =", s_false))
  }
})

test_that("truncated likelihood does not run away when an accumulator saturates", {
  des <- small_threshold_design()
  p <- small_threshold_pars(des)
  fit <- make_small_threshold_dadm(des, p, truncate = TRUE)

  ll_of <- function(s_false) {
    q <- p
    q["s_lMFALSE"] <- log(s_false)
    as.numeric(EMC2:::calc_ll_manager(
      matrix(q, nrow = 1, dimnames = list(NULL, names(q))),
      dadm = fit$dadm, model = fit$model))
  }

  # A total log-likelihood may legitimately be positive -- these are RT
  # densities -- so the contract is not a sign but a shape: b/s -> 0 is a
  # BOUNDED ridge, because conditioning on a response inside [LT, UT] cancels
  # the vanishing survivor between numerator and normaliser.  The likelihood
  # must therefore settle onto a plateau no better than the sane value, never
  # climb past it.  Before the fix the numerator lost the survivor while the
  # normaliser kept it, and this ladder spiked ~125 nats above ll_sane.
  ll_sane <- ll_of(1)
  expect_true(is.finite(ll_sane))

  ladder <- vapply(10^c(1, 2, 3, 4, 6, 8, 10, 12, 14), ll_of, numeric(1))
  expect_false(anyNA(ladder))
  expect_true(all(ladder <= ll_sane + 1e-8))      # no rung beats the sane fit
  expect_true(all(diff(ladder) <= 1e-4))          # never climbs back up

  # The plateau itself, read below 1e8.  Past there a loser row falls through
  # the min_ll floor and the trial is scored at the floor, which is a
  # deliberate cliff rather than part of the ridge.
  expect_equal(ladder[6], ladder[5], tolerance = 1e-3)
})
