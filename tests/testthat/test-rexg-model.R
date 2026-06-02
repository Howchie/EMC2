test_that("REXG dfun/pfun match ex-Gaussian references", {
  rt <- c(-0.2, 0.05, 0.2, 0.5, 1.0, 2.0)
  mu <- c(-0.1, 0.0, 0.1, 0.2, 0.3, 0.4)
  sigma <- c(0.15, 0.2, 0.12, 0.18, 0.25, 0.3)
  tau <- c(0.2, 0.25, 0.3, 0.18, 0.22, 0.35)
  pars <- cbind(mu = mu, sigma = sigma, tau = tau)
  model <- REXG()

  expect_equal(model$dfun(rt, pars), EMC2:::dexGaussian(rt, pars), tolerance = 1e-12)
  expect_equal(model$pfun(rt, pars), EMC2:::pexGaussian(rt, pars), tolerance = 1e-12)
})

test_that("REXG evaluates censored race likelihood in C++ path", {
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  design_rexg <- design(
    factors = list(subjects = 1, S = 1),
    Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1),
    functions = list(lI = lI_fun),
    model = REXG
  )

  p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18))
  dat <- make_data(p, design_rexg, n_trials = 80, rt_resolution = 1e-3)
  datc <- make_missing(dat, LC = 0.2, UC = 0.9, verbose = FALSE, rt_resolution = 1e-3)
  emc <- make_emc(datc, design_rexg, type = "single")
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }

  ll <- EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = p_types,
    min_ll = log(1e-10),
    trend = model$trend
  )

  expect_true(is.finite(ll))
})
