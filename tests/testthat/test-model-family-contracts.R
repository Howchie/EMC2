# Minimum-viable regression coverage for the numerical model families.
#
# The detailed independent quadrature, convergence, Monte Carlo, and solver
# validation lives in the files guarded by skip_model_validation().  These
# checks protect the package contracts that a routine refactor can break.

mvp_cdf_pdf_contract <- function(label, cdf, pdf, times, endpoint = NULL,
                                  defective = FALSE) {
  testthat::expect_true(all(is.finite(cdf)), label)
  testthat::expect_true(all(cdf >= -1e-12 & cdf <= 1 + 1e-12), label)
  testthat::expect_true(all(diff(cdf) >= -1e-10), label)
  testthat::expect_true(all(is.finite(pdf)) && all(pdf >= -1e-12), label)
  if (!is.null(endpoint)) {
    testthat::expect_equal(cdf[length(cdf)], endpoint, tolerance = 1e-10,
                           info = label)
  }
  if (defective) testthat::expect_lt(cdf[length(cdf)], 1, info = label)
}

test_that("model constructors expose one stable race-model contract", {
  models <- list(
    BAwF = BAwF(), BAwD = BAwD(), BAwR = BAwR(),
    BTAwL = BTAwLTransient(), FRQ = FRQ(),
    ROU = ROU(), ROUp = ROUp(), RLF = RLF(), RDMSWTN = RDMSWTN()
  )
  expect_true(all(vapply(models, function(m) identical(m$type, "RACE"), logical(1))))
  expect_true(all(vapply(models, function(m) is.character(m$c_name) && length(m$c_name) == 1L,
                         logical(1))))
  expect_true(all(vapply(models, function(m) "t0" %in% m$p_types_canonical,
                         logical(1))))
  threshold_models <- models[c("BAwF", "BAwD", "BAwR", "BTAwL", "ROU", "ROUp", "RLF", "RDMSWTN")]
  expect_true(all(vapply(threshold_models, function(m)
    all(c("B", "A") %in% m$p_types_canonical), logical(1))))
  expect_true(all(vapply(models, function(m) is.function(m$dfun) && is.function(m$pfun),
                         logical(1))))
  expect_match(BAwF()$c_name, "^BAwF")
  expect_match(BAwD()$c_name, "^BAwD")
  expect_match(BAwR()$c_name, "^BAwR")
  expect_match(BTAwLTransient()$c_name, "^BTAwL")

  # eta is an optional trailing parameter: its zero default preserves existing
  # designs, while omitting it from the canonical set keeps that default silent.
  ballistic <- list(
    LBA = LBA(),
    BAwL = BAwL(),
    BAwD = BAwD(),
    BAwDp = BAwDp(),
    BAwF = BAwF(),
    BAwR = BAwR(),
    BTAwL = BTAwL(),
    BTAwLTransient = BTAwLTransient(),
    BTAwLSustained = BTAwLSustained()
  )
  expect_true(all(vapply(ballistic, function(m) "eta" %in% names(m$p_types),
                         logical(1))))
  expect_true(all(vapply(ballistic, function(m) identical(unname(m$p_types[["eta"]]), 0),
                         logical(1))))
  expect_true(all(vapply(ballistic, function(m) !"eta" %in% m$p_types_canonical,
                         logical(1))))
})

test_that("representative model kernels preserve the CDF/PDF contract", {
  cases <- list(
    BAwF = list(
      times = EMC2:::bawf_tmax(.3, 1.2, 1.2, Inf) * c(.05, .3, .7, 1, 1.2),
      cdf = function(t) EMC2:::pbawf(t, .3, 1.2, 1, .5, 1.2, 1L),
      pdf = function(t) EMC2:::dbawf(t, .3, 1.2, 1, .5, 1.2, 1L)
    ),
    BAwD = list(
      times = EMC2:::bawd_tmax(.3, 1.2, .8, .4) * c(.05, .3, .7, 1, 1.2),
      cdf = function(t) EMC2:::pbawd(t, .3, 1.2, 1, .5, .8, .4, 1L),
      pdf = function(t) EMC2:::dbawd(t, .3, 1.2, 1, .5, .8, .4, 1L)
    ),
    BAwR = list(
      times = EMC2:::bawr_tmax(.3, 1.2, 1.2, 1) * c(.05, .3, .7, 1, 1.2),
      cdf = function(t) EMC2:::pbawr(t, .3, 1.2, 1, .5, 1.2, 1, 1L),
      pdf = function(t) EMC2:::dbawr(t, .3, 1.2, 1, .5, 1.2, 1, 1L)
    ),
    BTAwL = list(
      times = EMC2:::btawl_tmax_vec(.5, 1) * c(.05, .3, .7, 1, 1.2),
      cdf = function(t) EMC2:::pbtawl_transient(t, .3, 1.3, 1, .5, .5, 1, 0L),
      pdf = function(t) EMC2:::dbtawl_transient(t, .3, 1.3, 1, .5, .5, 1, 0L)
    ),
    FRQ = list(
      times = c(.01, .1, .4, 1, 4),
      cdf = function(t) EMC2:::pfrq(t, 2, 3, .8, .4),
      pdf = function(t) EMC2:::dfrq(t, 2, 3, .8, .4)
    )
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    cdf <- cs$cdf(cs$times)
    pdf <- cs$pdf(cs$times)
    mvp_cdf_pdf_contract(nm, cdf, pdf, cs$times)
    if (nm %in% c("BAwF", "BAwD", "BAwR", "BTAwL")) {
      expect_equal(cs$cdf(cs$times[4:5]), rep(cs$cdf(cs$times[4]), 2),
                   tolerance = 1e-10, info = paste(nm, "endpoint"))
      expect_equal(cs$pdf(cs$times[5]), 0, tolerance = 1e-12,
                   info = paste(nm, "endpoint density"))
    } else {
      expect_equal(cs$cdf(Inf), .8, tolerance = 1e-12,
                   info = "FRQ defective endpoint")
    }
  }
})

test_that("OU, Levy, and RDMSWTN wrappers preserve finite monotone outputs", {
  times <- c(.05, .2, .5, 1)
  rou <- cbind(v = rep(1.5, 4), k = rep(.5, 4), B = 1, A = .2,
               t0 = 0, s = 1)
  roup <- cbind(v_S = rep(1, 4), v_T = rep(.5, 4), tau_S = rep(2, 4),
                tau_T = rep(1, 4), k = rep(.5, 4), B = 1, A = .2,
                t0 = 0, s = 1)
  rlf <- cbind(v = 1.5, B = 1, A = .2, t0 = 0, s = 1, alpha = 1.7)
  rdmswtn <- cbind(v = 1.5, b = 1.2, A = .2, t0 = 0, s = 1, sv = .3,
                   lambda_g = 0, lambda_k = 0)
  for (x in list(
    ROU = list(p = pROU(times, rou), d = dROU(times, rou)),
    ROUp = list(p = pROUp(times, roup), d = dROUp(times, roup)),
    RLF = list(p = pRLF(times, rlf), d = dRLF(times, rlf)),
    RDMSWTN = list(p = pRDMSWTN(times, rdmswtn), d = dRDMSWTN(times, rdmswtn))
  )) {
    mvp_cdf_pdf_contract(names(x), x$p, x$d, times)
  }
  expect_equal(pROU(Inf, rou[1, , drop = FALSE]), 1, tolerance = 1e-8)
  expect_equal(pROUp(Inf, roup[1, , drop = FALSE]), 1, tolerance = 1e-8)
  expect_equal(pRDMSWTN(Inf, rdmswtn[1, , drop = FALSE]), 1, tolerance = 1e-8)
})

test_that("ballistic decay limits retain the LBA route", {
  times <- c(.1, .3, .7, 1.2)
  lba_p <- EMC2:::plba(times, .3, 1.2, 1, .5, posdrift = FALSE)
  lba_d <- EMC2:::dlba(times, .3, 1.2, 1, .5, posdrift = FALSE)

  expect_equal(EMC2:::pbawf(times, .3, 1.2, 1, .5, 0, 0L,
                             posdrift = FALSE), lba_p, tolerance = 1e-10)
  expect_equal(EMC2:::dbawf(times, .3, 1.2, 1, .5, 0, 0L,
                             posdrift = FALSE), lba_d, tolerance = 1e-10)
  expect_equal(EMC2:::pbawr(times, .3, 1.2, 1, .5, 0, 1, 0L,
                             posdrift = FALSE), lba_p, tolerance = 1e-10)
  expect_equal(EMC2:::dbawr(times, .3, 1.2, 1, .5, 0, 1, 0L,
                             posdrift = FALSE), lba_d, tolerance = 1e-10)

  # For BAwD, k = 0 leaves a constant clearance ell, so the same identity is
  # obtained after shifting the launch by ell on the untruncated route.
  shifted_p <- EMC2:::plba(times, .3, 1.2, 1 - .4, .5, posdrift = FALSE)
  shifted_d <- EMC2:::dlba(times, .3, 1.2, 1 - .4, .5, posdrift = FALSE)
  expect_equal(EMC2:::pbawd(times, .3, 1.2, 1, .5, 0, .4, 0L,
                             posdrift = FALSE), shifted_p, tolerance = 1e-8)
  expect_equal(EMC2:::dbawd(times, .3, 1.2, 1, .5, 0, .4, 0L,
                             posdrift = FALSE), shifted_d, tolerance = 1e-8)
})

test_that("the representative FPE solver keeps its closed-form anchor", {
  times <- c(.05, .2, .5, 1)
  wald_cdf <- function(t) {
    pnorm((t - 1) / sqrt(t)) + exp(2) * pnorm(-(1 + t) / sqrt(t))
  }
  out <- EMC2:::fpe_bm_fht_pdf_cdf_vec(
    times, 1, 1, 0, 1, 1, 1, 1, 128L, 256L
  )
  expect_true(all(out$pdf >= 0))
  expect_true(all(diff(out$cdf) >= -1e-12))
  expect_lt(max(abs(out$cdf - wald_cdf(times))), 2e-3)
})

test_that("one compiled ballistic simulator path remains wired", {
  lR <- factor(rep(c("left", "right"), 8), levels = c("left", "right"))
  pars <- cbind(mu = 1, sigma = .5, b = 1.2, A = .3, t0 = .1, k = 1)
  pars <- pars[rep(1, length(lR)), , drop = FALSE]
  out <- EMC2:::.rfun_BAwF(lR, pars, rep(TRUE, nrow(pars)), launch = 1L)
  expect_named(out, c("R", "rt"))
  expect_equal(nrow(out), 8L)
  expect_true(all(is.factor(out$R)))
  expect_true(all(is.infinite(out$rt) | out$rt >= .1))
})
