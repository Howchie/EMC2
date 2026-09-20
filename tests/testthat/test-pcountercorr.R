library(EMC2)

test_that("PCOUNTERcorr exposes the shared-gamma fixed-K reference", {
  p <- rbind(c(nu = 8, c = .3, gamma = 0, k = 2, omega = 0, t0 = .1, rho = .45),
             c(nu = 5, c = .3, gamma = 0, k = 3, omega = 0, t0 = .1, rho = .45))
  expect_gt(dPCOUNTERcorr(.5, p, 1), 0)
  expect_equal(dPCOUNTERcorr(.5, p, 1), dPCOUNTERcorr(.5, p[2:1, ], 2), tolerance = 1e-12)
  # The two winner densities exhaust the response probability.
  total <- integrate(function(x) dPCOUNTERcorr(x, p, 1), .1, Inf)$value +
    integrate(function(x) dPCOUNTERcorr(x, p, 2), .1, Inf)$value
  expect_equal(total, 1, tolerance = 2e-7)
})

test_that("PCOUNTERcorr compiled likelihood agrees with the R reference", {
  matchfun <- function(d) factor(as.character(d$S) == as.character(d$lR),
                                 levels = c(FALSE, TRUE))
  des <- design(factors = list(subjects = 1, S = 1:2), Rlevels = 1:2,
    matchfun = matchfun,
    formula = list(nu ~ lM, c ~ 1, gamma ~ 1, k ~ 1, omega ~ 1, t0 ~ 1, rho ~ 1),
    constants = c(gamma = log(0), omega = log(0)), model = PCOUNTERcorr)
  theta <- c(nu = log(8), nu_lMTRUE = log(1.5), c = log(.3), k = log(2),
             t0 = log(.1), rho = qnorm(.4))
  set.seed(71)
  dat <- make_data(theta, des, n_trials = 30, rt_resolution = NULL)
  emc <- make_emc(dat, des, type = "single")
  dadm <- emc[[1]]$data[[1]]; mod <- emc[[1]]$model(); pt <- names(mod$p_types)
  designs <- lapply(pt, function(nm) { x <- attr(dadm, "designs")[[nm]]; x[attr(x, "expand"), , drop = FALSE] })
  names(designs) <- pt
  got <- EMC2:::calc_ll_oo(matrix(theta, 1, dimnames = list(NULL, names(theta))), dadm,
    attr(dadm, "constants"), designs, mod$c_name, mod$bound, mod$transform,
    mod$pre_transform, pt, log(1e-10), mod$trend)
  pm <- EMC2:::get_pars_c_wrapper_oo(matrix(theta, 1, dimnames = list(NULL, names(theta))),
    dadm, attr(dadm, "constants"), designs, mod$bound, mod$transform,
    mod$pre_transform, mod$trend)
  nacc <- length(levels(dadm$lR)); nuniq <- nrow(dadm) / nacc
  ref_unique <- vapply(seq_len(nuniq), function(j) {
    rows <- ((j - 1L) * nacc + 1L):(j * nacc)
    pp <- pm[rows, c("nu", "c", "gamma", "k", "omega", "t0", "rho"), drop = FALSE]
    win <- which(dadm$R[rows] == dadm$lR[rows])
    log(dPCOUNTERcorr(dadm$rt[rows[1]], pp, win))
  }, numeric(1))
  ref <- sum(ref_unique[attr(dadm, "expand")])
  expect_equal(as.numeric(got), ref, tolerance = 1e-8)
})
