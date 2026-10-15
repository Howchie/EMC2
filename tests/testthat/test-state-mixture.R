mf_sm <- function(d) d$S == d$lR
fac_sm <- list(subjects = 1, S = c("left", "right"))
form_sm <- list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1)

design_sm <- function(model, formula, constants = NULL) {
  suppressMessages(design(factors = fac_sm, Rlevels = c("left", "right"),
                          matchfun = mf_sm, formula = formula,
                          constants = c(sv = log(1), constants), model = model,
                          report_p_vector = FALSE))
}

pw_sm <- function(d, p, data) {
  emc <- suppressMessages(make_emc(data, d, type = "single", compress = FALSE,
                                   rt_resolution = NULL))
  as.vector(calc_ll_pw(matrix(p, 1, dimnames = list(NULL, names(p))),
                       emc[[1]]$data[[1]], emc[[1]]$model))
}

mix_sm <- function(l0, l1, p) log((1 - p) * exp(l0) + p * exp(l1))

plain_sm <- design_sm(LBA, form_sm)
target_sm <- design_sm(mix_states(LBA, "v"),
                       c(form_sm, list(v_mix ~ 1, pmix ~ 0 + lM)),
                       c(pmix_lMFALSE = qnorm(0)))
p0_sm <- c(v = 1, v_lMTRUE = 2, B = log(1), A = log(.5), t0 = log(.2))
p1_sm <- replace(p0_sm, "v_lMTRUE", 0.4)                 # target v = 1.4
pm_sm <- c(p0_sm, v_mix = 1.4, pmix_lMTRUE = qnorm(.3))

dat_sm <- local({
  set.seed(11)
  suppressMessages(make_data(pm_sm, target_sm, n_trials = 100))
})

test_that("an infinite constant on a cell-means column maps to its link limit", {
  mp <- mapped_pars(target_sm, pm_sm)
  expect_equal(mp$pmix[mp$lM == "TRUE"], c(.3, .3))
  expect_equal(mp$pmix[mp$lM == "FALSE"], c(0, 0))
})

test_that("pmix = 0 reproduces the parent model exactly", {
  p_off <- replace(pm_sm, "pmix_lMTRUE", -40)
  expect_identical(pw_sm(target_sm, p_off, dat_sm), pw_sm(plain_sm, p0_sm, dat_sm))
})

test_that("a mixing target racer gives the trial-level mixture", {
  ref <- mix_sm(pw_sm(plain_sm, p0_sm, dat_sm), pw_sm(plain_sm, p1_sm, dat_sm), .3)
  expect_equal(pw_sm(target_sm, pm_sm, dat_sm), ref, tolerance = 1e-12)
})

test_that("independently mixing racers give the four-state mixture", {
  both <- design_sm(mix_states(LBA, "v"), c(form_sm, list(v_mix ~ lM, pmix ~ 1)))
  pb <- c(p0_sm, v_mix = .5, v_mix_lMTRUE = 1.2, pmix = qnorm(.4))
  comp <- sapply(0:3, function(s) {
    tgt <- s %% 2; ntg <- s %/% 2
    vN <- if (ntg) .5 else 1
    vT <- if (tgt) 1.7 else 3
    q <- replace(p0_sm, c("v", "v_lMTRUE"), c(vN, vT - vN))
    log(ifelse(tgt, .4, .6) * ifelse(ntg, .4, .6)) + pw_sm(plain_sm, q, dat_sm)
  })
  ref <- apply(comp, 1, function(x) max(x) + log(sum(exp(x - max(x)))))
  expect_equal(pw_sm(both, pb, dat_sm), ref, tolerance = 1e-12)
})

test_that("alt = 'below' places the alternative state under the base value", {
  below <- design_sm(mix_states(LBA, "v", alt = "below"),
                     c(form_sm, list(v_below ~ 1, pmix ~ 0 + lM)),
                     c(pmix_lMFALSE = qnorm(0)))
  pb <- c(p0_sm, v_below = log(3 - 1.4), pmix_lMTRUE = qnorm(.3))
  ref <- mix_sm(pw_sm(plain_sm, p0_sm, dat_sm), pw_sm(plain_sm, p1_sm, dat_sm), .3)
  expect_equal(pw_sm(below, pb, dat_sm), ref, tolerance = 1e-12)
})

test_that("truncation normalises the mixed joint by the mixed window mass", {
  dt <- make_missing(dat_sm[dat_sm$rt > .25 & dat_sm$rt < 1.6, ], LT = .25, UT = 1.6)
  log_z <- function(p) {
    one <- function(rt, ...) {
      d <- dat_sm[1:2, ]
      d$S <- factor(c("left", "right"), levels = c("left", "right"))
      d$R <- factor(NA, levels = levels(dat_sm$R))
      d$rt <- rt
      make_missing(d, ...)
    }
    up <- exp(pw_sm(plain_sm, p, one(Inf, UC = 1.6)))
    lo <- exp(pw_sm(plain_sm, p, one(-Inf, LC = .25)))
    stats::setNames(log(1 - up - lo), c("left", "right"))[as.character(dt$S)]
  }
  z0 <- log_z(p0_sm); z1 <- log_z(p1_sm)
  ref <- mix_sm(pw_sm(plain_sm, p0_sm, dt) + z0, pw_sm(plain_sm, p1_sm, dt) + z1, .3) -
    mix_sm(z0, z1, .3)
  expect_equal(pw_sm(target_sm, pm_sm, dt), unname(ref), tolerance = 1e-9)
})

test_that("mix_states refuses models without independent accumulators", {
  expect_error(mix_states(DDM, "v")(), "independent race")
  expect_error(mix_states(LBA, "nope")(), "not parameters")
})
