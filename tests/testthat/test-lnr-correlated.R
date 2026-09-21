test_that("LNRcorr exposes the Reynolds Gaussian-copula contract", {
  model <- LNRcorr()
  expect_equal(model$c_name, "LNR_CORR")
  expect_identical(model$correlation_type, "lnr_gaussian_copula")
  expect_true(isTRUE(model$correlated))
  expect_equal(unname(model$transform$lower[["rho"]]), -1)
  expect_equal(unname(model$transform$upper[["rho"]]), 1)
  expect_error(
    model$Ttransform(data.frame(m = 0, s = 1, t0 = 0, rho = 1.1), NULL),
    "rho"
  )
})

test_that("LNRcorr likelihood agrees with the bivariate normal conditional form", {
  matchfun <- function(d) factor(as.character(d$S) == as.character(d$lR),
                                  levels = c(FALSE, TRUE))
  des <- design(
    factors = list(subjects = 1, S = c("a", "b")),
    Rlevels = c("a", "b"), matchfun = matchfun,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1,
                   rho = as.formula("rho ~ 1")),
    constants = c(t0 = log(0)), model = LNRcorr
  )
  dat <- data.frame(
    subjects = 1,
    S = factor(c("a", "b"), levels = c("a", "b")),
    R = factor(c("a", "b"), levels = c("a", "b")),
    rt = c(0.4, 0.5)
  )
  emc <- make_emc(dat, des, type = "single", compress = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[["rho"]] <- qnorm(0.7)
  designs <- lapply(names(model$p_types), function(nm) {
    d <- attr(dadm, "designs")[[nm]]
    d[attr(d, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  got <- EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound,
    transforms = model$transform, pretransforms = model$pre_transform,
    p_types = names(model$p_types), min_ll = log(1e-10), trend = model$trend
  )
  rho <- 0.4
  expected <- sum(vapply(seq_along(dat$rt), function(i) {
    z <- log(dat$rt[i])
    winner <- i
    loser <- 3L - winner
    cm <- rho * z
    cs <- sqrt(1 - rho^2)
    dlnorm(dat$rt[i], 0, 1, log = TRUE) +
      pnorm((z - cm) / cs, lower.tail = FALSE, log.p = TRUE)
  }, numeric(1)))
  expect_equal(as.numeric(got), expected, tolerance = 1e-10)
})

test_that("the compiled LNRcorr simulator follows the R reference", {
  skip_on_cran()
  n <- 10000L
  lR <- factor(rep(c("a", "b"), n), levels = c("a", "b"))
  pars <- cbind(m = rep(c(-0.6, -0.2), n), s = rep(c(0.5, 0.55), n),
                t0 = 0, rho = 0.4)
  ok <- rep(TRUE, nrow(pars))
  set.seed(1)
  cpp <- EMC2:::rlnr_corr_cpp(pars, levels(lR), ok)
  options_old <- options(emc2.cpp_rfun = FALSE)
  on.exit(options(options_old), add = TRUE)
  set.seed(2)
  ref <- EMC2:::.rLNRcorr(lR, as.data.frame(pars), ok)
  expect_lt(abs(mean(cpp$R == 1) - mean(ref$R == "a")), 0.03)
  expect_lt(abs(median(cpp$rt) - median(ref$rt)), 0.03)
})
