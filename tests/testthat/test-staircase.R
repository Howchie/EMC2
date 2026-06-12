# Staircase SSD simulation via the `staircase=` argument of make_data(),
# through both the conditional (block rfun) and unconditional (trial-by-trial)
# simulation paths, for SSEXG and SSRDEX.

stair_list <- list(SSD0 = .25, stairstep = .05, stairmin = 0, stairmax = Inf)

stair_design <- function(model, label) {
  design(model = model,
    factors  = list(subjects = 1, S = c("left", "right")),
    Rlevels  = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(
      lI  = function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2),
      SSD = function(d) SSD_function(d, SSD = NA, pSSD = .25)),
    formula  = if (identical(label, "SSRDEX"))
        list(v ~ 0 + lM, B ~ 1, A ~ 1, t0 ~ 1,
             muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
      else
        list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
             muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
  )
}

stair_pvector <- function(des, label) {
  p <- sampled_pars(des, doMap = FALSE)
  if (identical(label, "SSRDEX")) {
    p["v_lMFALSE"] <- log(0.2); p["v_lMTRUE"] <- log(1.5)
    p["B"] <- log(1); p["A"] <- log(0.4); p["t0"] <- log(0.15)
  } else {
    p["mu_lMFALSE"] <- log(0.6); p["mu_lMTRUE"] <- log(0.45)
    p["sigma"] <- log(0.05); p["tau"] <- log(0.1)
  }
  p["muS"] <- log(0.2); p["sigmaS"] <- log(0.08); p["tauS"] <- log(0.1)
  p["gf"] <- qnorm(0.08); p["tf"] <- qnorm(0.08)
  p
}

expect_staircase_ok <- function(dat) {
  ssd <- dat$SSD[is.finite(dat$SSD)]
  expect_gt(length(ssd), 0)
  # stop-trial proportion near pSSD = .25
  expect_gt(length(ssd) / nrow(dat), .15)
  expect_lt(length(ssd) / nrow(dat), .35)
  # staircase actually moves and starts at SSD0
  expect_gt(length(unique(round(ssd, 3))), 1)
  expect_equal(ssd[1], stair_list$SSD0)
  # consecutive staircase SSDs step by exactly +/- stairstep
  expect_true(all(round(abs(diff(ssd)), 3) %in% c(0, stair_list$stairstep)))
}

for (label in c("SSEXG", "SSRDEX")) {
  model <- if (label == "SSRDEX") SSRDEX else SSEXG
  des <- stair_design(model, label)
  p   <- stair_pvector(des, label)

  test_that(paste0(label, ": conditional staircase simulation works"), {
    set.seed(11)
    dat <- suppressMessages(make_data(p, des, n_trials = 200, staircase = stair_list))
    expect_staircase_ok(dat)
  })

  test_that(paste0(label, ": unconditional staircase simulation works"), {
    set.seed(12)
    dat <- suppressMessages(make_data(p, des, n_trials = 100, staircase = stair_list,
                                      conditional_on_data = FALSE))
    expect_staircase_ok(dat)
  })

  test_that(paste0(label, ": unconditional vectorised staircase simulation works"), {
    set.seed(13)
    dat <- suppressMessages(make_data(p, des, n_trials = 100, staircase = stair_list,
                                      conditional_on_data = FALSE, use_vectorised = TRUE))
    expect_staircase_ok(dat)
  })
}

test_that("unconditional simulation with NA SSDs and no staircase errors clearly", {
  des <- stair_design(SSRDEX, "SSRDEX")
  p   <- stair_pvector(des, "SSRDEX")
  expect_error(
    suppressMessages(make_data(p, des, n_trials = 10, conditional_on_data = FALSE)),
    "staircase list must be supplied")
})

# --- Unconditional branch: censoring (make_missing) ------------------------

test_that("unconditional simulation applies UC censoring and staircase reacts", {
  des <- stair_design(SSEXG, "SSEXG")
  p   <- stair_pvector(des, "SSEXG")
  set.seed(21)
  dat <- suppressMessages(make_data(p, des, n_trials = 100, staircase = stair_list,
                                    conditional_on_data = FALSE,
                                    TC = list(UC = .7, UCresponse = TRUE)))
  # censoring applied: no finite rt above UC, censored trials coded Inf/NA
  expect_true(all(dat$rt[is.finite(dat$rt)] <= .7))
  expect_true(any(is.infinite(dat$rt)))
  expect_true(all(is.na(dat$R[is.infinite(dat$rt)])))
  # ladder mechanics still intact
  expect_staircase_ok(dat)
})

test_that("harsh UC drives the unconditional staircase up (censored = stop success)", {
  des <- stair_design(SSEXG, "SSEXG")
  p   <- stair_pvector(des, "SSEXG")
  set.seed(22)
  dat_uc <- suppressMessages(make_data(p, des, n_trials = 150, staircase = stair_list,
                                       conditional_on_data = FALSE, TC = list(UC = .35)))
  set.seed(22)
  dat_no <- suppressMessages(make_data(p, des, n_trials = 150, staircase = stair_list,
                                       conditional_on_data = FALSE))
  ssd_uc <- dat_uc$SSD[is.finite(dat_uc$SSD)]
  ssd_no <- dat_no$SSD[is.finite(dat_no$SSD)]
  expect_gt(median(ssd_uc), median(ssd_no))
})

# --- Unconditional branch: grouped make_ssd() staircases --------------------

grouped_ssd_fun <- make_ssd(factors = "S", p_stop = .5,
                            staircase = list(left  = list(SSD0 = .2),
                                             right = list(SSD0 = .4)))

expect_grouped_ladders_ok <- function(dat) {
  for (lvl in c("left", "right")) {
    ssd <- dat$SSD[dat$S == lvl & is.finite(dat$SSD)]
    expect_gt(length(ssd), 1)
    expect_equal(ssd[1], if (lvl == "left") .2 else .4)
    expect_true(all(round(abs(diff(ssd)), 3) %in% c(0, .05)))
  }
}

for (uv in c(FALSE, TRUE)) {
  test_that(paste0("grouped make_ssd staircase, unconditional",
                   if (uv) " vectorised" else ""), {
    des <- design(model = SSEXG,
      factors  = list(subjects = 1, S = c("left", "right")),
      Rlevels  = c("left", "right"),
      matchfun = function(d) as.character(d$S) == as.character(d$lR),
      functions = list(
        lI = function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)),
      formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                      muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
    p <- stair_pvector(des, "SSEXG")
    set.seed(31 + uv)
    dat <- suppressMessages(make_data(p, des, n_trials = 150,
                                      functions = list(SSD = grouped_ssd_fun),
                                      conditional_on_data = FALSE,
                                      use_vectorised = uv))
    expect_grouped_ladders_ok(dat)
  })
}

test_that("custom staircase functions error under conditional_on_data = FALSE", {
  des <- stair_design(SSRDEX, "SSRDEX")
  p   <- stair_pvector(des, "SSRDEX")
  sc  <- stair_list
  attr(sc, "staircase_function") <- function(dts, spec) NULL
  expect_error(
    suppressMessages(make_data(p, des, n_trials = 10, staircase = sc,
                               conditional_on_data = FALSE)),
    "Custom staircase functions")
})
