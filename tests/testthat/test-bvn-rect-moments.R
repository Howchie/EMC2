library(EMC2)

rect_moment_reference <- function(mu1, sd1, mu2, sd2, rho,
                                  lo1, hi1, lo2, hi2) {
  sd_cond <- sd2 * sqrt(1 - rho^2)
  zl1 <- if (is.infinite(lo1)) -Inf else (lo1 - mu1) / sd1
  zh1 <- if (is.infinite(hi1)) Inf else (hi1 - mu1) / sd1

  integrand <- function(x, what) {
    v1 <- mu1 + sd1 * x
    cm <- mu2 + rho * sd2 * x
    zl2 <- (lo2 - cm) / sd_cond
    zh2 <- (hi2 - cm) / sd_cond
    p2 <- pnorm(zh2) - pnorm(zl2)
    m2 <- cm * p2 + sd_cond * (dnorm(zl2) - dnorm(zh2))
    w <- dnorm(x)
    switch(what,
           p = w * p2,
           m1 = w * v1 * p2,
           m2 = w * m2,
           m12 = w * v1 * m2)
  }

  c(p = integrate(integrand, zl1, zh1, what = "p",
                  subdivisions = 300L, rel.tol = 1e-9)$value,
    m1 = integrate(integrand, zl1, zh1, what = "m1",
                   subdivisions = 300L, rel.tol = 1e-9)$value,
    m2 = integrate(integrand, zl1, zh1, what = "m2",
                   subdivisions = 300L, rel.tol = 1e-9)$value,
    m12 = integrate(integrand, zl1, zh1, what = "m12",
                    subdivisions = 300L, rel.tol = 1e-9)$value)
}

test_that("BVN rectangle probabilities and moments match conditional integration", {
  skip_on_cran()
  rectangles <- list(
    c(-Inf, .5, -Inf, .7),
    c(-.3, .4, .2, 1.2),
    c(.1, .11, .6, .8),
    c(.5, Inf, -.2, Inf)
  )
  for (rho in c(-.8, 0, .8)) {
    for (mu in list(c(0, 0), c(.4, -.3))) {
      for (bounds in rectangles) {
        got <- EMC2:::bawl_corr_bvn_rect_probe(
          mu[1], .9, mu[2], 1.1, rho,
          bounds[1], bounds[2], bounds[3], bounds[4]
        )
        expect_equal(got$status, 0L)
        ref <- rect_moment_reference(
          mu[1], .9, mu[2], 1.1, rho,
          bounds[1], bounds[2], bounds[3], bounds[4]
        )
        expect_equal(unname(unlist(got[c("p", "m1", "m2", "m12")])),
                     unname(ref), tolerance = 1e-4)
      }
    }
  }
})

test_that("exact pair values agree with the independent numeric-pair route", {
  skip_on_cran()
  for (pos in c(FALSE, TRUE)) {
    for (rho in c(-.8, -.3, .3, .8)) {
      for (ks in list(c(0, 0), c(1e-6, .3), c(.3, 1e-6))) {
        exact <- EMC2:::bawl_corr_pair_probe(
          .7, .1, .2, 1.4, ks[1], 1, .8,
          .13, .35, 1.1, ks[2], .4, 1.1,
          rho, pos, FALSE
        )
        numeric <- EMC2:::bawl_corr_pair_probe(
          .7, .1, .2, 1.4, ks[1], 1, .8,
          .13, .35, 1.1, ks[2], .4, 1.1,
          rho, pos, TRUE
        )
        expect_true(numeric$survival_status %in% c(0L, 1L))
        expect_true(numeric$cause1_status %in% c(0L, 1L))
        expect_true(numeric$cause2_status %in% c(0L, 1L))
        if (exact$survival_status == 0L && numeric$survival_status == 0L)
          expect_equal(exact$survival, numeric$survival, tolerance = 2e-5)
        if (exact$cause1_status == 0L && numeric$cause1_status == 0L)
          expect_equal(exact$cause1, numeric$cause1, tolerance = 2e-5)
        if (exact$cause2_status == 0L && numeric$cause2_status == 0L)
          expect_equal(exact$cause2, numeric$cause2, tolerance = 2e-5)
      }
    }
  }
})

test_that("positive pair normalizer is the BVN orthant probability", {
  skip_on_cran()
  cases <- list(
    c(1.0, .8, .4, 1.1, -.8),
    c(.2, .9, -.3, 1.2, 0),
    c(1.5, .7, .1, .6, .8)
  )
  for (z in cases) {
    got <- EMC2:::bawl_corr_pair_probe(
      .7, .1, .2, 1.4, 0, z[1], z[2],
      .13, .35, 1.1, 0, z[3], z[4],
      z[5], TRUE, FALSE
    )$normalizer
    ref <- rect_moment_reference(
      z[1], z[2], z[3], z[4], z[5], 0, Inf, 0, Inf
    )["p"]
    expect_equal(got, unname(ref), tolerance = 1e-6)
  }
})

test_that("point-start pair branches remain continuous and finite", {
  skip_on_cran()
  for (starts in list(c(0, .2), c(.2, 0), c(0, 0))) {
    for (pos in c(FALSE, TRUE)) {
      got <- EMC2:::bawl_corr_pair_probe(
        .7, .1, starts[1], 1.4, .2, 1, .8,
        .13, starts[2], 1.1, .1, .4, 1.1,
        .6, pos, FALSE
      )
      ref <- EMC2:::bawl_corr_pair_probe(
        .7, .1, starts[1], 1.4, .2, 1, .8,
        .13, starts[2], 1.1, .1, .4, 1.1,
        .6, pos, TRUE
      )
      expect_true(got$survival_status %in% c(0L, 1L))
      expect_true(got$cause1_status %in% c(0L, 1L))
      expect_true(got$cause2_status %in% c(0L, 1L))
      expect_equal(got$survival, ref$survival, tolerance = 3e-5)
      expect_equal(got$cause1, ref$cause1, tolerance = 3e-5)
      expect_equal(got$cause2, ref$cause2, tolerance = 3e-5)
    }
  }
})

make_exact_context <- function(dat, posdrift = TRUE, rho = .6) {
  model <- BAwLcorr(posdrift = posdrift)
  formula <- list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1,
                  k ~ 1, mG ~ 1, mK ~ 1, rho ~ 1)
  des <- design(
    data = dat, Rlevels = levels(dat$R),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    formula = formula, model = model, report_p_vector = FALSE
  )
  emc <- make_emc(dat, des, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE, rt_resolution = NULL)
  dadm <- emc[[1]]$data[[1]]
  model_fit <- emc[[1]]$model()
  designs <- lapply(names(model_fit$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model_fit$p_types)
  p <- sampled_pars(des, doMap = FALSE)
  p["v"] <- 1
  p[c("sv", "B", "A", "t0", "k", "mG", "mK")] <-
    log(c(.8, 1.4, .2, .1, 0, 1, 1))
  p["rho"] <- qnorm((rho + 1) / 2)
  list(
    args = list(
      particle_matrix = matrix(p, nrow = 1,
                               dimnames = list(NULL, names(p))),
      data = dadm, constants = attr(dadm, "constants"), designs = designs,
      type = model_fit$c_name, bounds = model_fit$bound,
      transforms = model_fit$transform, pretransforms = model_fit$pre_transform,
      p_types = names(model_fit$p_types), min_ll = log(1e-10),
      trend = model_fit$trend
    ),
    data = dadm
  )
}

test_that("exact pair routing covers finite, censored, missing, and pointwise paths", {
  skip_on_cran()
  cases <- list(
    data.frame(subjects = factor(1), S = factor("correct", levels = c("correct", "error")),
               R = factor("correct", levels = c("correct", "error")), rt = .7),
    data.frame(subjects = factor(1), S = factor("correct", levels = c("correct", "error")),
               R = factor(NA, levels = c("correct", "error")), rt = Inf,
               LC = 0, UC = .4),
    data.frame(subjects = factor(1), S = factor("correct", levels = c("correct", "error")),
               R = factor(NA, levels = c("correct", "error")), rt = -Inf,
               LC = .4),
    data.frame(subjects = factor(1), S = factor("correct", levels = c("correct", "error")),
               R = factor(NA, levels = c("correct", "error")), rt = NA_real_,
               LT = .2, LC = .6, UC = .8, UT = 1.5),
    data.frame(subjects = factor(1), S = factor("correct", levels = c("correct", "error")),
               R = factor("correct", levels = c("correct", "error")), rt = NA_real_,
               LT = .2, LC = .6, UC = .8, UT = 1.5)
  )
  for (dat in cases) {
    ctx <- make_exact_context(dat)
    EMC2:::bawl_corr_counters_reset()
    Sys.setenv(EMC2_BAWLCORR_COUNTERS = "1")
    on.exit(Sys.unsetenv("EMC2_BAWLCORR_COUNTERS"), add = TRUE)
    total <- as.numeric(do.call(EMC2:::calc_ll_oo, ctx$args))
    pointwise <- as.numeric(do.call(EMC2:::calc_ll_oo_pw, ctx$args))
    counters <- unlist(EMC2:::bawl_corr_counter_values())
    expect_true(is.finite(total))
    expect_true(all(is.finite(pointwise)))
    expect_equal(sum(pointwise), total, tolerance = 1e-10)
    expect_gt(unname(counters[["exact_pair_trials"]] +
                     counters[["numeric_pair_trials"]]), 0)
  }
})
