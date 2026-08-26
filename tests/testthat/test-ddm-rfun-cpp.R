test_that("DDM constructor uses the compiled bounded single-process simulator", {
  old_opt <- getOption("emc2.cpp_rfun")
  on.exit(options(emc2.cpp_rfun = old_opt), add = TRUE)
  options(emc2.cpp_rfun = TRUE)

  n <- 3000L
  raw <- matrix(
    rep(c(.8, 1.3, .45, .2, .1, 1.1, .35, .24), n),
    nrow = n, byrow = TRUE,
    dimnames = list(NULL, c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"))
  )
  pars <- DDM()$Ttransform(raw, NULL)
  attr(pars, "ok") <- rep(TRUE, n)
  R <- factor(rep(c("left", "right"), length.out = n),
              levels = c("left", "right"))

  set.seed(41)
  direct <- EMC2:::rddm_cpp(pars, levels(R), attr(pars, "ok"))
  set.seed(41)
  through_model <- DDM()$rfun(list(R = R), pars)
  expect_equal(as.integer(through_model$R), direct$R)
  expect_equal(through_model$rt, direct$rt)
  expect_true(all(is.finite(direct$rt)))

  options(emc2.cpp_rfun = FALSE)
  set.seed(42)
  fallback <- DDM()$rfun(list(R = R), pars)
  set.seed(42)
  reference <- EMC2:::rDDM(R, pars, attr(pars, "ok"))
  expect_equal(fallback, reference)

  options(emc2.cpp_rfun = TRUE)
  set.seed(43)
  compiled <- EMC2:::rddm_cpp(pars, levels(R), attr(pars, "ok"))
  expect_lt(abs(mean(reference$R == "left") - mean(compiled$R == 1)), .05)
  expect_lt(
    suppressWarnings(ks.test(reference$rt[reference$R == "left"],
                              compiled$rt[compiled$R == 1])$statistic),
    .08
  )
  expect_lt(
    suppressWarnings(ks.test(reference$rt[reference$R == "right"],
                              compiled$rt[compiled$R == 2])$statistic),
    .08
  )
})

test_that("DDM endpoint probabilities come from one bounded diffusion", {
  n <- 5000L
  raw <- matrix(
    rep(c(1, 1, 0, 0, 0, 1, .5, 0), n),
    nrow = n, byrow = TRUE,
    dimnames = list(NULL, c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"))
  )
  pars <- DDM()$Ttransform(raw, NULL)
  set.seed(44)
  simulated <- EMC2:::rddm_cpp(
    pars, c("left", "right"), rep(TRUE, n)
  )
  expected_lower <- EMC2:::pDDM(
    Inf, factor("lower", levels = c("lower", "upper")), pars[1, , drop = FALSE]
  )[1]
  expect_lt(abs(mean(simulated$R == 1) - expected_lower), .03)
  expect_true(all(is.finite(simulated$rt)))
})

test_that("DDMGNG reuses the compiled bounded DDM path", {
  old_opt <- getOption("emc2.cpp_rfun")
  on.exit(options(emc2.cpp_rfun = old_opt), add = TRUE)
  options(emc2.cpp_rfun = TRUE)

  n <- 200L
  raw <- matrix(
    rep(c(1.2, 1, 0, .1, 0, 1, .5, 0), n),
    nrow = n, byrow = TRUE,
    dimnames = list(NULL, c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"))
  )
  R <- factor(rep(c("no", "yes"), length.out = n), levels = c("no", "yes"))
  pars <- DDMGNG()$Ttransform(
    raw,
    data.frame(TIMEOUT = rep(.4, n), Rnogo = rep(1, n))
  )
  attr(pars, "ok") <- rep(TRUE, n)

  set.seed(45)
  out <- DDMGNG()$rfun(list(R = R), pars)
  expect_equal(nrow(out), n)
  expect_true(all(is.na(out$R) == is.infinite(out$rt)))
  expect_true(all(out$rt[is.finite(out$rt)] <= .4))
})
