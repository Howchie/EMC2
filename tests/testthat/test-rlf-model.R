skip_model_validation()

rlf_pars <- function(n, v = 1.5, B = 1, A = 0, t0 = 0, s = 1,
                     alpha = 1.7) {
  cbind(
    v = rep(v, length.out = n),
    B = rep(B, length.out = n),
    A = rep(A, length.out = n),
    t0 = rep(t0, length.out = n),
    s = rep(s, length.out = n),
    alpha = rep(alpha, length.out = n)
  )
}

test_that("RLF cache keys exclude t0 and scale out s", {
  n <- 40L
  rt <- seq(0.2, 1.5, length.out = n)
  scale <- rep(c(1, 2), length.out = n)
  out <- EMC2:::rlf_pdf_cdf_vec(
    rt, 1.5 * scale, 1 * scale, 0.2 * scale,
    rep(c(0, 0.05), length.out = n), scale, rep(1.7, n),
    50L, 0.02, FALSE, TRUE, TRUE
  )
  expect_identical(out$n_keys, 1L)
  expect_identical(out$n_solves, 1L)
  expect_true(all(is.finite(out$pdf)))
  expect_true(all(out$cdf >= 0 & out$cdf <= 1))

  alpha <- rep(c(1.5, 1.7), length.out = n)
  two <- EMC2:::rlf_pdf_cdf_vec(
    rt, 1.5 * scale, 1 * scale, 0.2 * scale,
    rep(0, n), scale, alpha,
    50L, 0.02, FALSE, TRUE, TRUE
  )
  expect_identical(two$n_keys, 2L)
  expect_identical(two$n_solves, 2L)
})

test_that("RLF production grid defaults stay synchronized", {
  withr::local_options(list(
    emc2.rlf_nx = NULL,
    emc2.rlf_dt = NULL,
    emc2.rlf_richardson = NULL,
    emc2.rlf_richardson_ratio = NULL
  ))
  grid <- EMC2:::.rlf_grid()
  expect_identical(grid$nx, 70L)
  expect_identical(grid$dt_target, 1.6e-2)
  expect_true(grid$richardson)
  expect_identical(grid$richardson_ratio, 1.4)

  direct <- formals(EMC2:::rlf_pdf_cdf_vec)
  expect_identical(direct$nx, 70L)
  expect_identical(direct$dt_target, 1.6e-2)
  expect_true(eval(direct$sparse_output))
  expect_true(eval(direct$richardson))
  expect_false(eval(direct$adaptive))
  expect_identical(direct$richardson_ratio, 1.4)
})

test_that("RLF production grid controls low-alpha discretisation error", {
  rt <- seq(0.05, 4, length.out = 31)
  n <- length(rt)
  solve <- function(nx, dt, ratio) {
    EMC2:::rlf_pdf_cdf_vec(
      rt, rep(1.5, n), rep(1.2, n), rep(0.2, n), rep(0, n),
      rep(1, n), rep(1.1, n),
      nx, dt, FALSE, TRUE,
      horizon_split = TRUE, richardson = TRUE,
      richardson_ratio = ratio
    )
  }
  production <- solve(160L, 0.016, 1.25)
  former <- solve(128L, 0.008, 1.5)
  reference <- solve(320L, 0.004, 1.5)
  central <- reference$cdf >= 0.02 & reference$cdf <= 0.98 &
    reference$pdf > 1e-10
  log_pdf_error <- function(x) {
    mean(abs(log(x$pdf[central]) - log(reference$pdf[central])))
  }

  expect_lt(log_pdf_error(production), 0.01)
  expect_lt(max(abs(production$cdf - reference$cdf)), 0.006)
  expect_lt(log_pdf_error(production), log_pdf_error(former))
  expect_lt(
    max(abs(production$cdf - reference$cdf)),
    max(abs(former$cdf - reference$cdf))
  )
})

test_that("RLF sparse and complete-grid output agree", {
  rt <- c(0.08, 0.17, 0.29, 0.44, 0.63, 0.86, 1.13, 1.5)
  v <- seq(1, 2.4, length.out = 8)
  B <- rep(c(0.9, 1.1), 4)
  A <- rep(c(0, 0.25), 4)
  t0 <- rep(c(0, 0.03), 4)
  s <- rep(1, 8)
  alpha <- seq(1.3, 2, length.out = 8)
  call <- function(sparse) {
    EMC2:::rlf_pdf_cdf_vec(
      rt, v, B, A, t0, s, alpha,
      60L, 0.02, FALSE, sparse,
      horizon_split = TRUE, richardson = FALSE
    )
  }

  # Sparse output evaluates the modal curves exactly at the requested times;
  # the complete grid samples them on a uniform grid and interpolates, so the
  # two agree only to the resolution of that grid.
  sparse <- call(TRUE)
  full <- call(FALSE)
  expect_true(all(is.finite(sparse$pdf)))
  expect_true(all(sparse$pdf >= 0))
  expect_true(all(sparse$cdf >= 0 & sparse$cdf <= 1))
  expect_equal(sparse$cdf, full$cdf, tolerance = 5e-3)
  expect_equal(sparse$pdf, full$pdf, tolerance = 5e-2)
})

test_that("RLF sparse output does not depend on the output-grid spacing", {
  # dt_target sets the spacing of the cached curve, not a time step: with the
  # query times supplied it is not consulted at all, and the answer must be
  # bit-for-bit independent of it.
  rt <- seq(0.1, 1.5, length.out = 16)
  call <- function(dt) {
    EMC2:::rlf_pdf_cdf_vec(
      rt, rep(1.5, 16), rep(1, 16), rep(0.2, 16),
      rep(0, 16), rep(1, 16), rep(1.7, 16),
      60L, dt, FALSE, TRUE, TRUE
    )
  }
  coarse <- call(0.08)
  fine <- call(0.005)
  expect_identical(coarse$n_solves, 1L)
  expect_equal(coarse$pdf, fine$pdf, tolerance = 1e-14)
  expect_equal(coarse$cdf, fine$cdf, tolerance = 1e-14)
})

test_that("RLF alpha = 2 approaches the independent RDM oracle", {
  rt <- seq(0.1, 2, by = 0.1)
  pars <- rlf_pars(length(rt), alpha = 2)
  levy <- EMC2:::rlf_pdf_cdf_vec(
    rt, pars[, "v"], pars[, "B"], pars[, "A"], pars[, "t0"],
    pars[, "s"], pars[, "alpha"],
    200L, 0.005, FALSE, TRUE, TRUE
  )
  rdm <- pars[, c("v", "B", "A", "t0", "s"), drop = FALSE]
  expect_lt(max(abs(levy$pdf - EMC2:::dRDM(rt, rdm))), 3e-3)
  expect_lt(max(abs(levy$cdf - EMC2:::pRDM(rt, rdm))), 3e-4)
})

test_that("RLF model transforms alpha and reaches the C++ race likelihood", {
  withr::local_options(list(
    emc2.rlf_nx = 40L,
    emc2.rlf_dt = 0.02
  ))
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[
    forstmann$subjects %in% unique(forstmann$subjects)[1],
  ][1:12, ]
  dat$subjects <- droplevels(dat$subjects)
  design_rlf <- suppressMessages(design(
    data = dat, model = RLF, matchfun = matchfun,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
    constants = c(s = log(1), alpha = qnorm(0.7))
  ))
  emc <- suppressMessages(make_emc(
    dat, design_rlf, type = "single", n_chains = 1, compress = TRUE
  ))
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- lapply(p_types, function(name) {
    matrix <- attr(dadm, "designs")[[name]]
    matrix[attr(matrix, "expand"), , drop = FALSE]
  })
  names(designs) <- p_types

  p <- c(
    v = log(1.6), v_lMTRUE = 0.5, B = log(1),
    A = log(0.2), t0 = log(0.15)
  )[names(sampled_pars(design_rlf))]
  particle <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  mapped <- EMC2:::get_pars_c_wrapper_oo(
    particle, dadm, attr(dadm, "constants"), designs,
    model$bound, model$transform, model$pre_transform,
    return_all_pars = TRUE
  )
  expect_equal(unique(mapped[, "alpha"]), 1.7)

  ll <- EMC2:::calc_ll_oo(
    particle, dadm, attr(dadm, "constants"), designs,
    model$c_name, model$bound, model$transform, model$pre_transform,
    p_types, log(1e-10), model$trend
  )
  expect_true(is.finite(ll))
  expect_gt(ll, nrow(dat) * log(1e-10))
})
