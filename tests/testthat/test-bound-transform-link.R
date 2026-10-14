mf_bt <- function(d) d$S == d$lR
dat_bt <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ]
dat_bt$subjects <- droplevels(dat_bt$subjects)

lba_bt <- function(..., data = dat_bt) {
  suppressMessages(design(data = data, model = LBA, matchfun = mf_bt,
                          report_p_vector = FALSE, ...))
}

ll_bt <- function(d, p, data = dat_bt) {
  emc <- suppressMessages(make_emc(data, d, type = "single", compress = FALSE,
                                   rt_resolution = NULL))
  EMC2:::calc_ll_manager(matrix(p, 1, dimnames = list(NULL, names(p))),
                         emc[[1]]$data[[1]], emc[[1]]$model)
}

test_that("the transform takes the bound's lower side", {
  d <- lba_bt(formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1))
  tf <- d$model()$transform
  expect_equal(unname(tf$lower[c("t0", "A", "sv")]), c(0.05, 1e-4, 1e-4))
  p <- sampled_pars(d)
  p[] <- 0
  for (x in c(-30, -5, 0, 3)) {
    p["t0"] <- x
    mp <- mapped_pars(d, p, digits = 15)
    expect_true(all(mp$t0 > 0.05))
    expect_equal(unique(mp$t0), 0.05 + exp(x))
  }
})

test_that("parameters with a constant keep the plain link and its exceptions", {
  d <- lba_bt(formula = list(v ~ lM, sv ~ lM, B ~ 1, A ~ 1),
              constants = c(sv = log(1), A = log(0)))
  tf <- d$model()$transform
  expect_equal(unname(tf$lower[c("sv", "A", "t0", "pContaminant")]),
               c(0, 0, 0, 0))
  p <- sampled_pars(d)
  p[] <- 0
  mp <- mapped_pars(d, p)
  expect_true(all(mp$t0 == 0))
  expect_true(all(mp$A == 0))
  expect_true(all(mp$sv[mp$lM == "FALSE"] == 1))
  expect_true(is.finite(ll_bt(d, p)))
  expect_gt(ll_bt(d, p), length(dat_bt$rt) * log(1e-10))
})

test_that("bound overrides move the link and explicit transforms win", {
  d <- lba_bt(formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
              bound = list(minmax = cbind(t0 = c(0.1, Inf))))
  expect_equal(unname(d$model()$transform$lower["t0"]), 0.1)
  d <- lba_bt(formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
              transform = list(lower = c(t0 = 0)))
  expect_equal(unname(d$model()$transform$lower["t0"]), 0)
  expect_error(lba_bt(formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                      transform = list(func = c(t0 = "pnorm"),
                                       lower = c(t0 = 0), upper = c(t0 = 0.04))),
               "lies outside")
})

test_that("explicit model sides are kept (FRQ shapes, RLF alpha)", {
  tf <- EMC2:::fill_transform(NULL, FRQ, bound = EMC2:::fill_bound(NULL, FRQ))
  expect_equal(unname(tf$lower[c("alpha", "beta", "t0")]), c(1, 1, 0.05))
  tf <- EMC2:::fill_transform(NULL, RLF, bound = EMC2:::fill_bound(NULL, RLF))
  expect_equal(unname(c(tf$lower["alpha"], tf$upper["alpha"])), c(1, 2))
})

test_that("the linked and plain designs agree at the same natural values", {
  f <- list(v ~ lM, sv ~ lM, B ~ E + lR, A ~ 1, t0 ~ 1)
  d_plain <- lba_bt(formula = f, constants = c(sv = log(1)),
                    transform = list(lower = c(A = 0, t0 = 0)))
  d_link <- lba_bt(formula = f, constants = c(sv = log(1)))
  p <- sampled_pars(d_plain)
  p[] <- 0
  p[c("v", "v_lMTRUE", "sv_lMTRUE", "A", "t0")] <-
    c(log(2), log(1.5), log(0.8), log(0.5), log(0.25))
  q <- p
  q["A"] <- log(0.5 - 1e-4)
  q["t0"] <- log(0.25 - 0.05)
  expect_equal(ll_bt(d_link, q), ll_bt(d_plain, p))
})

test_that("a posttransform trend that crosses the bound is still floored", {
  dat <- dat_bt
  dat$cov <- seq(0, 1, length.out = nrow(dat))
  tr <- make_trend(par_names = "t0", cov_names = "cov", kernels = "lin_incr",
                   bases = "lin", phase = "posttransform")
  d <- suppressMessages(design(data = dat, model = LBA, matchfun = mf_bt,
                               formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                               trend = tr, report_p_vector = FALSE))
  expect_equal(unname(d$model()$transform$lower["t0"]), 0.05)
  p <- sampled_pars(d)
  p[] <- 0
  p["t0"] <- log(0.25)
  p["t0.w"] <- -1
  bad <- ll_bt(d, p, data = dat)
  # Floored trials contribute log(1e-10) each; no other trial's density
  # exceeds 10 at these values.
  n_out <- sum(0.3 - dat$cov <= 0.05)
  expect_lt(bad, n_out * log(1e-10) + (nrow(dat) - n_out) * log(10))
})

test_that("saved fits keep the transform they were made with", {
  tf <- samples_LNR[[1]]$model()$transform
  expect_equal(unname(tf$lower["t0"]), 0)
})

test_that("a sampled -Inf reaches the plain endpoint and its exception", {
  d <- lba_bt(formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1))
  tf <- d$model()$transform
  expect_equal(unname(tf$lower_inf[c("A", "t0")]), c(0, 0))
  p <- sampled_pars(d)
  p[] <- 0
  p["t0"] <- log(0.2)
  p["A"] <- -Inf
  expect_true(all(mapped_pars(d, p, digits = 15)$A == 0))
  expect_gt(ll_bt(d, p), length(dat_bt$rt) * log(1e-10))
  # A finite value that rounds onto the linked side is not the exception.
  p["A"] <- -50
  expect_true(all(mapped_pars(d, p, digits = 15)$A == 1e-4))
  expect_equal(ll_bt(d, p), nrow(dat_bt) * log(1e-10))
})

test_that("split transforms use the same link and endpoint as the plain path", {
  f <- list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1)
  d <- lba_bt(formula = f)
  d_split <- lba_bt(formula = f, pre_transform_terms = list(t0 = "t0", A = "A"))
  dadm <- EMC2:::design_model(dat_bt, d, compress = FALSE, verbose = FALSE)
  dadm_split <- EMC2:::design_model(dat_bt, d_split, compress = FALSE, verbose = FALSE)
  p <- sampled_pars(d)
  p[] <- 0
  for (x in c(-3, -Inf)) {
    p["t0"] <- x
    p["A"] <- x
    m <- EMC2:::get_pars_matrix_oo(p, dadm, d$model())
    m_split <- EMC2:::get_pars_matrix_oo(p, dadm_split, d_split$model())
    expect_equal(m_split[, c("t0", "A")], m[, c("t0", "A")])
  }
  expect_true(all(m_split[, c("t0", "A")] == 0))
})
