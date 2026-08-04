test_that("Gompertz is exposed with the physical parameter contract", {
  m <- EMC2::GOM()
  expect_identical(m$c_name, "GOM")
  expect_identical(names(m$p_types)[1:6],
                   c("alpha", "beta", "K", "B", "A", "t0"))
  expect_identical(m$p_types_canonical, names(m$p_types)[1:6])
})

test_that("Gompertz uses physical start and moving-boundary scales", {
  withr::local_options(emc2.fpe_nx = 64L, emc2.fpe_dt = 0.02,
                       emc2.fpe_grade = 4, emc2.fpe_tgrade = 8)
  rt <- c(0.2, 0.5, 1.0, Inf)
  pars <- cbind(alpha = rep(1, length(rt)), beta = rep(0.5, length(rt)),
                K = rep(5, length(rt)), B = rep(1, length(rt)),
                A = rep(0.5, length(rt)), t0 = rep(0.1, length(rt)))
  fixed <- EMC2:::pGOM(rt, pars)
  expect_true(all(diff(fixed[is.finite(fixed)]) >= 0))
  expect_equal(fixed[length(fixed)], 1)

  collapsing <- cbind(pars, Binf = rep(1, length(rt)), tau = rep(1, length(rt)))
  moving <- EMC2:::pGOM(rt, collapsing, kind = "exponential")
  expect_gt(moving[3], fixed[3])
})

test_that(".rGOM_R reads collapsing boundary parameters from matrix colnames", {
  lR <- factor(c("left", "right"), levels = c("left", "right"))
  pars <- cbind(alpha = c(1, 1), beta = c(0.5, 0.5), K = c(5, 5),
                B = c(1, 1), A = c(0, 0), t0 = c(0.1, 0.1),
                Binf = c(0.2, 0.2), tau = c(0.1, 0.1), pw = c(1, 1))
  set.seed(42)
  sim_fixed <- EMC2:::.rGOM_R(lR, pars, kind = "fixed")
  set.seed(42)
  sim_coll <- EMC2:::.rGOM_R(lR, pars, kind = "exponential")
  # Rapid collapse to lower Binf should shorten RTs compared to fixed bound
  expect_lt(sim_coll$rt[1], sim_fixed$rt[1])
})
