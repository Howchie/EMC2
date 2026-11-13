# calc_ll_oo keeps the data-fixed race state (response coding, RT bins,
# winner/loser masks, the mixed path's shared state) across calls, keyed by a
# content hash of the fields it is built from.  Cached and uncached likelihoods
# must be identical, and any change to those fields must rebuild the state.

race_cache_fixture <- function(censor = TRUE) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = matchfun, model = LNR,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1))
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- c(-0.5, -0.5, log(0.8), log(0.2))
  set.seed(11)
  dat <- make_data(p, des, n_trials = 150)
  if (censor) {
    dat <- make_missing(dat, UC = quantile(dat$rt, 0.85), rt_resolution = 0.001,
                        verbose = FALSE)
  }
  emc <- suppressMessages(make_emc(dat, des, type = "single", rt_resolution = 0.001))
  dadm <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]], force_rebuild = TRUE)
  set.seed(12)
  P <- matrix(rep(p, each = 9) + rnorm(9 * length(p), sd = 0.05), 9,
              dimnames = list(NULL, names(p)))
  list(dadm = dadm, model = emc[[1]]$model, P = P)
}

race_ll <- function(fx, dadm = fx$dadm, P = fx$P) {
  EMC2:::calc_ll_manager(P, dadm, fx$model, ll_cache_validated = TRUE)
}

with_cache_off <- function(expr) {
  old <- Sys.getenv("EMC2_LL_DATA_CACHE_OFF", unset = NA)
  Sys.setenv(EMC2_LL_DATA_CACHE_OFF = "1")
  on.exit(if (is.na(old)) Sys.unsetenv("EMC2_LL_DATA_CACHE_OFF") else
    Sys.setenv(EMC2_LL_DATA_CACHE_OFF = old), add = TRUE)
  force(expr)
}

test_that("cached race likelihoods equal uncached ones", {
  for (censor in c(TRUE, FALSE)) {
    fx <- race_cache_fixture(censor)
    EMC2:::emc2_ll_data_cache_clear()
    first <- race_ll(fx)
    again <- race_ll(fx)
    one <- race_ll(fx, P = fx$P[1, , drop = FALSE])
    off <- with_cache_off(race_ll(fx))
    expect_identical(again, first)
    expect_identical(off, first)
    expect_identical(one, first[1])
    expect_true(all(is.finite(first)))
  }
})

test_that("changed data rebuilds the cached race state", {
  fx <- race_cache_fixture(TRUE)
  EMC2:::emc2_ll_data_cache_clear()
  base <- race_ll(fx)
  changed <- fx$dadm
  finite <- is.finite(changed$rt)
  changed$rt[finite] <- changed$rt[finite] + 0.002
  expect_identical(race_ll(fx, dadm = changed),
                   with_cache_off(race_ll(fx, dadm = changed)))
  expect_false(identical(race_ll(fx, dadm = changed), base))
  # A changed censoring bound is part of the cached state as well.
  bound <- fx$dadm
  bound$UC <- bound$UC + 0.01
  expect_identical(race_ll(fx, dadm = bound),
                   with_cache_off(race_ll(fx, dadm = bound)))
  # The original data still finds its own entry.
  expect_identical(race_ll(fx), base)
})

test_that("modifying data already seen rebuilds the cached race state", {
  fx <- race_cache_fixture(TRUE)
  EMC2:::emc2_ll_data_cache_clear()
  d <- fx$dadm
  d$rt <- d$rt + 0
  first <- race_ll(fx, dadm = d)
  finite <- is.finite(d$rt)
  d$rt[finite] <- d$rt[finite] + 0.002
  attr(d, "expand") <- attr(d, "expand")
  expect_identical(race_ll(fx, dadm = d), with_cache_off(race_ll(fx, dadm = d)))
  expect_false(identical(race_ll(fx, dadm = d), first))
})

test_that("the race state cache stays bounded", {
  fx <- race_cache_fixture(FALSE)
  EMC2:::emc2_ll_data_cache_clear()
  for (k in 1:4200) {
    d <- fx$dadm
    d$rt <- d$rt + k * 1e-6
    race_ll(fx, dadm = d, P = fx$P[1, , drop = FALSE])
  }
  # 4096 entries at most; clearing removes them all.
  n_cleared <- EMC2:::emc2_ll_data_cache_clear()
  expect_gt(n_cleared, 64L)
  expect_lte(n_cleared, 4096L)
  expect_identical(EMC2:::emc2_ll_data_cache_clear(), 0L)
})
