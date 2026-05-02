test_that("RGAMMA dfun/pfun match gamma references with shift accounted for", {
  rt <- c(0.05, 0.2, 0.5, 1.0, 2.0, 4.0)
  lambda <- c(0.7, 1.3, 2.1, 0.9, 1.7, 3.0)
  shift <- c(0.01, 0.05, 0.1, 0.2, 0.15, 0.3)
  model <- RGAMMA()

  for (shape in c(1, 2)) {
    pars <- cbind(lambda = lambda, shape = rep(shape, length(rt)), shift = shift)
    d_mod <- model$dfun(rt, pars)
    p_mod <- model$pfun(rt, pars)

    tt <- rt - shift
    d_ref <- ifelse(tt > 0, stats::dgamma(tt, shape = shape, rate = lambda), 0)
    p_ref <- ifelse(tt > 0, stats::pgamma(tt, shape = shape, rate = lambda), 0)

    expect_equal(d_mod, d_ref, tolerance = 1e-12)
    expect_equal(p_mod, p_ref, tolerance = 1e-12)

    # Shift accounting check:
    # f_shift(rt) = dgamma(rt - shift), F_shift(rt) = pgamma(rt - shift), 0 below shift.
    shift2 <- 0.3
    pars0 <- cbind(lambda = lambda, shape = rep(shape, length(rt)), shift = rep(0, length(rt)))
    rt_shifted <- rt + shift2
    d_shift_mod <- model$dfun(rt_shifted, cbind(pars0[, c("lambda", "shape"), drop = FALSE], shift = rep(shift2, length(rt))))
    p_shift_mod <- model$pfun(rt_shifted, cbind(pars0[, c("lambda", "shape"), drop = FALSE], shift = rep(shift2, length(rt))))
    d_shift_ref <- stats::dgamma(rt_shifted - shift2, shape = shape, rate = lambda)
    p_shift_ref <- stats::pgamma(rt_shifted - shift2, shape = shape, rate = lambda)

    expect_equal(d_shift_mod, d_shift_ref, tolerance = 1e-12)
    expect_equal(p_shift_mod, p_shift_ref, tolerance = 1e-12)

    rt_below_shift <- c(0.0, shift2 / 2)
    pars_short <- pars[seq_along(rt_below_shift), , drop = FALSE]
    pars_short[, "shift"] <- shift2
    expect_equal(model$dfun(rt_below_shift, pars_short), c(0, 0))
    expect_equal(model$pfun(rt_below_shift, pars_short), c(0, 0))
  }
})
