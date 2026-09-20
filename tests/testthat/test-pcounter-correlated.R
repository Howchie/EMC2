# Regression tests for the two-racer common-CV correlated-rate PCOUNTER.

pcorr_test_design <- function() {
  matchfun <- function(d) factor(as.character(d$S) == as.character(d$lR),
                                 levels = c(FALSE, TRUE))
  design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2,
    matchfun = matchfun,
    formula = list(nu ~ lM, c ~ 1, gamma ~ 1, k ~ 1,
                   omega ~ 1, t0 ~ 1, rho ~ 1),
    constants = c(gamma = log(0), omega = log(0)),
    model = PCOUNTERcorr
  )
}

pcorr_test_ll <- function(p, des, dat) {
  emc <- make_emc(dat, des, type = "single")
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  designs <- lapply(names(model$p_types), function(nm) {
    d <- attr(dadm, "designs")[[nm]]
    d[attr(d, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  list(
    ll = EMC2:::calc_ll_oo(
      matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
      constants = attr(dadm, "constants"), designs = designs,
      type = model$c_name, bounds = model$bound,
      transforms = model$transform, pretransforms = model$pre_transform,
      p_types = names(model$p_types), min_ll = log(1e-10),
      trend = model$trend
    ),
    dadm = dadm,
    pars = EMC2:::get_pars_matrix_oo(p, dadm, model)
  )
}

test_that("PCOUNTERcorr exposes the common-CV two-racer contract", {
  model <- PCOUNTERcorr()
  expect_equal(model$c_name, "PCOUNTERcorr")
  expect_true(isTRUE(model$correlated))
  expect_equal(model$p_types_canonical,
               c("nu", "c", "gamma", "k", "omega", "t0", "rho"))
  expect_equal(unname(model$transform$lower[["rho"]]), 0)
  expect_equal(unname(model$transform$upper[["rho"]]), 1)
  # K=2 is the default fixed counter threshold in the PM application.  The
  # exact raw constant k=0 therefore has to remain a valid model boundary,
  # just as it is for ordinary PCOUNTER.
  expect_true(0 %in% unname(model$bound$exception[names(model$bound$exception) == "k"]))
  expect_error(
    dPCOUNTERcorr(0.5, rbind(
      c(nu = 8, c = .3, gamma = 0, k = 1, omega = 0, t0 = 0, rho = .4),
      c(nu = 8, c = .2, gamma = 0, k = 1, omega = 0, t0 = 0, rho = .4)
    )),
    "common c/rho"
  )
  expect_error(
    dPCOUNTERcorr(0.5, rbind(
      c(nu = 8, c = .3, gamma = .1, k = 1, omega = 0, t0 = 0, rho = .4),
      c(nu = 8, c = .3, gamma = .1, k = 1, omega = 0, t0 = 0, rho = .4)
    )),
    "gamma=omega=0"
  )
})

test_that("rho zero nests the independent gamma-rate counter race", {
  pc <- PCOUNTER()
  pars_c <- rbind(
    c(nu = 8, c = .35, gamma = 0, k = 1, omega = 0, t0 = .1, rho = 0),
    c(nu = 8, c = .35, gamma = 0, k = 1, omega = 0, t0 = .1, rho = 0)
  )
  pars_i <- rbind(
    c(nu = 8, sv = 2.8, gamma = 0, k = 1, omega = 0, t0 = .1),
    c(nu = 8, sv = 2.8, gamma = 0, k = 1, omega = 0, t0 = .1)
  )
  rt <- c(.15, .2, .4, .8, 1.5)
  got <- dPCOUNTERcorr(rt, pars_c, winner = 1L)
  ref <- pc$dfun(rt, pars_i[rep(1, length(rt)), , drop = FALSE]) *
    (1 - pc$pfun(rt, pars_i[rep(2, length(rt)), , drop = FALSE]))
  expect_equal(got, ref, tolerance = 1e-11)

  # The point-rate boundary is the ordinary Erlang race.
  pars_c[, "c"] <- 0
  pars_c[, "rho"] <- 0
  fixed <- stats::dgamma(rt - .1, shape = 3, rate = 8) *
    (1 - stats::pgamma(rt - .1, shape = 3, rate = 8))
  expect_equal(dPCOUNTERcorr(rt, pars_c, winner = 1L), fixed,
               tolerance = 1e-11)
})

test_that("compiled correlated PCOUNTER likelihood matches the analytic R density",
          {
            des <- pcorr_test_design()
            p <- c(nu = log(8), nu_lMTRUE = log(1.4), c = log(.35),
                   k = log(1), t0 = log(.1), rho = qnorm(.4))
            set.seed(22)
            dat <- make_data(p, des, n_trials = 20, rt_resolution = NULL)
            out <- pcorr_test_ll(p, des, dat)
            dadm <- out$dadm
            pm <- out$pars
            direct <- vapply(seq(1, nrow(pm), by = 2), function(j) {
              winner <- which(dadm$winner[j:(j + 1)])
              log(dPCOUNTERcorr(dadm$rt[j], pm[j:(j + 1), , drop = FALSE],
                                winner = winner))
            }, numeric(1))
            expect_true(all(is.finite(direct)))
            expected <- sum(direct[attr(dadm, "expand")])
            expect_equal(unname(out$ll), expected, tolerance = 1e-10)

            # Pointwise likelihood is the same kernel and must preserve the
            # trial grouping used by WAIC.
            emc <- make_emc(dat, des, type = "single")
            dadm2 <- emc[[1]]$data[[1]]
            model <- emc[[1]]$model()
            designs <- lapply(names(model$p_types), function(nm) {
              d <- attr(dadm2, "designs")[[nm]]
              d[attr(d, "expand"), , drop = FALSE]
            })
            names(designs) <- names(model$p_types)
            pw <- EMC2:::calc_ll_oo_pw(
              matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm2,
              constants = attr(dadm2, "constants"), designs = designs,
              type = model$c_name, bounds = model$bound,
              transforms = model$transform, pretransforms = model$pre_transform,
              p_types = names(model$p_types), min_ll = log(1e-10),
              trend = model$trend
            )
            expect_equal(sum(pw), expected, tolerance = 1e-10)
          })
