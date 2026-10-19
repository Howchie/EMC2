# Gamma rate mixtures (Math/gamma_mixtures.tex) for OPR (ratio chart, and the
# trial gain of the background chart) and RSUB (both evidence laws).  References integrate the fixed-rate kernels
# against the Gamma(1/cv^2, mean v) law of the trial rate by adaptive
# quadrature between its quantiles.

gm_ref <- function(fn, v, cv) {
  beta <- 1 / cv^2
  q <- qgamma(c(1e-18, seq(.02, .98, .02), 1 - 1e-16), beta, rate = beta / v)
  sum(vapply(seq_len(length(q) - 1), function(k)
    integrate(function(V) dgamma(V, beta, rate = beta / v) * fn(V), q[k], q[k + 1],
              rel.tol = 1e-13, subdivisions = 2000)$value, 0))
}

test_that("OPR mixture density, CDF and survivor match quadrature over the rate", {
  skip_on_cran()
  cases <- list(c(3, 10, 1, .3), c(1, 5, .3, .8), c(2.5, 8, 2, .2),
                c(7.3, 40, .5, .05), c(1.5, 3, 1, 1.5))
  for (cs in cases) {
    K <- cs[1]; v <- cs[2]; rho <- cs[3]; cv <- cs[4]
    for (x in K / v * c(.1, .5, 1, 2, 6)) {
      info <- paste(c(cs, x), collapse = "/")
      expect_equal(EMC2:::dopr_gamma(x, K, v, rho, cv),
                   gm_ref(function(V) EMC2:::dopr(x, K, (1 + rho) * V, rho * V), v, cv),
                   tolerance = 1e-9, info = info)
      expect_equal(EMC2:::popr_gamma(x, K, v, rho, cv),
                   gm_ref(function(V) EMC2:::popr(x, K, (1 + rho) * V, rho * V), v, cv),
                   tolerance = 1e-9, info = info)
      expect_equal(EMC2:::popr_gamma(x, K, v, rho, cv) +
                     EMC2:::popr_gamma(x, K, v, rho, cv, lower_tail = FALSE), 1,
                   tolerance = 1e-13, info = info)
    }
  }
})

test_that("the OPR mixture density is the hypergeometric closed form", {
  hyp <- function(x, K, v, rho, cv) {
    beta <- 1 / cv^2; rv <- beta / v
    zz <- 4 * rho * (1 + rho) * x^2 / (rv + (1 + 2 * rho) * x)^2
    a <- (beta + K) / 2; b <- a + .5; c <- K + 1
    t <- 1; s <- 1
    for (n in 0:1e5) {
      t <- t * (a + n) * (b + n) / ((c + n) * (n + 1)) * zz; s <- s + t
      if (t < 1e-17 * s) break
    }
    lgamma(beta + K) - lgamma(beta) - lgamma(K) + beta * log(rv) + K * log1p(rho) +
      (K - 1) * log(x) - (beta + K) * log(rv + (1 + 2 * rho) * x) + log(s)
  }
  for (cs in list(c(3, 10, 1, .3), c(12, 30, .2, .1), c(2, 6, 3, .7)))
    for (x in cs[1] / cs[2] * c(.1, 1, 5))
      expect_equal(EMC2:::dopr_gamma(x, cs[1], cs[2], cs[3], cs[4], log_out = TRUE),
                   hyp(x, cs[1], cs[2], cs[3], cs[4]), tolerance = 1e-11)
})

test_that("RSUB mixtures match quadrature over the rate", {
  skip_on_cran()
  cases <- list(c(1.5, 2, 0, .7, .4), c(1, .5, 0, .7, .7), c(1.3, .9, .5, 1, .3),
                c(1.5, 2, 0, .7, .1), c(1, 0, .8, 1, .5), c(1, 1, .3, .5, 1.3))
  for (packets in c(FALSE, TRUE)) for (cs in cases) {
    v <- cs[1]; B <- cs[2]; A <- cs[3]; s <- cs[4]; cv <- cs[5]
    for (x in (B + A / 2) / v * c(.05, .5, 1, 2, 6)) {
      info <- paste(c(packets, cs, x), collapse = "/")
      expect_equal(EMC2:::drsub(x, v, B, A, s, packets, cv = cv),
                   gm_ref(function(V) EMC2:::drsub(rep(x, length(V)), V, B, A, s, packets),
                          v, cv), tolerance = 1e-8, info = info)
      expect_equal(EMC2:::prsub(x, v, B, A, s, packets, cv = cv),
                   gm_ref(function(V) EMC2:::prsub(rep(x, length(V)), V, B, A, s, packets),
                          v, cv), tolerance = 1e-8, info = info)
      expect_equal(EMC2:::prsub(x, v, B, A, s, packets, cv = cv) +
                     EMC2:::prsub(x, v, B, A, s, packets, lower_tail = FALSE, cv = cv), 1,
                   tolerance = 1e-12, info = info)
    }
  }
})

test_that("the background-chart gain matches quadrature over the gain", {
  skip_on_cran()
  cases <- list(c(3, 10, 4, .3), c(1.4, 4, 7, .8), c(2.6, 8, 0, .5))
  for (cs in cases) {
    K <- cs[1]; v <- cs[2]; lambda0 <- cs[3]; cv <- cs[4]
    pars <- cbind(K = K, v = v, lambda0 = lambda0, t0 = 0, cv = cv)
    for (x in K / v * c(.1, .5, 1, 2, 6)) {
      info <- paste(c(cs, x), collapse = "/")
      expect_equal(EMC2:::dOPR(x, pars),
                   gm_ref(function(G) EMC2:::dopr(x, K, G * (v + lambda0), G * lambda0), 1, cv),
                   tolerance = 1e-9, info = info)
      expect_equal(EMC2:::pOPR(x, pars),
                   gm_ref(function(G) EMC2:::popr(x, K, G * (v + lambda0), G * lambda0), 1, cv),
                   tolerance = 1e-9, info = info)
    }
  }
})

test_that("small cv approaches the fixed rate", {
  u <- c(.1, .3, 1)
  expect_equal(EMC2:::popr_gamma(u, 3, 10, 1, 1e-3), EMC2:::popr(u, 3, 20, 10),
               tolerance = 1e-5)
  for (packets in c(FALSE, TRUE))
    expect_equal(EMC2:::prsub(u, 1.3, .9, .5, .8, packets, cv = 1e-3),
                 EMC2:::prsub(u, 1.3, .9, .5, .8, packets), tolerance = 1e-5)
})

test_that("the simulators reproduce the mixture CDFs", {
  skip_on_cran()
  local_test_rng(kind = "Mersenne-Twister")
  n <- 20000
  probe <- c(.3, .5, .8, 1.2, 2)
  lR <- factor(rep("r", n))
  po <- cbind(K = 2.6, v = 8, rho = .7, t0 = .1, cv = .5)[rep(1, n), ]
  th <- EMC2:::popr_gamma(probe - .1, 2.6, 8, .7, .5)
  for (cpp in c(TRUE, FALSE)) {
    m <- if (cpp) n else 4000
    sim <- withr::with_options(list(emc2.cpp_rfun = cpp),
                               EMC2:::.rfun_OPR(lR[seq_len(m)], po[seq_len(m), ], 0L))
    emp <- vapply(probe, function(x) mean(sim$rt <= x), 0)
    expect_lt(max(abs(emp - th)), if (cpp) .015 else .03)
  }
  pb <- cbind(K = 2.6, v = 8, lambda0 = 5.6, t0 = .1, cv = .5)[rep(1, n), ]
  for (cpp in c(TRUE, FALSE)) {
    m <- if (cpp) n else 4000
    sim <- withr::with_options(list(emc2.cpp_rfun = cpp),
                               EMC2:::.rfun_OPR(lR[seq_len(m)], pb[seq_len(m), ], 2L))
    emp <- vapply(probe, function(x) mean(sim$rt <= x), 0)
    expect_lt(max(abs(emp - th)), if (cpp) .015 else .03)
  }
  pr <- cbind(v = 1.3, B = .9, A = 1.1, t0 = .2, s = .8, cv = .6)[rep(1, n), ]
  for (packets in c(FALSE, TRUE)) for (cpp in c(TRUE, FALSE)) {
    sim <- withr::with_options(list(emc2.cpp_rfun = cpp),
                               EMC2:::.rfun_RSUB(lR, pr, packets))
    th <- EMC2:::prsub(probe - .2, 1.3, .9, 1.1, .8, packets, cv = .6)
    emp <- vapply(probe, function(x) mean(sim$rt <= x), 0)
    expect_lt(max(abs(emp - th)), .015)
  }
})

gm_fixture <- function(model, formula, p) {
  set.seed(20261017)
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  invisible(utils::capture.output(
    des <- suppressMessages(design(data = dat, model = model, formula = formula,
                                   matchfun = function(d) d$S == d$lR))))
  invisible(utils::capture.output(p0 <- sampled_pars(des)))
  expect_setequal(names(p0), names(p))
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  p <- p[names(p0)]
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  designs <- lapply(setNames(nm = names(model$p_types)), function(nm)
    attr(dadm, "designs")[[nm]][attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE])
  ll <- EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                          constants = attr(dadm, "constants"), designs = designs,
                          type = model$c_name, bounds = model$bound,
                          transforms = model$transform,
                          pretransforms = model$pre_transform,
                          p_types = names(model$p_types), min_ll = log(1e-10),
                          trend = model$trend)
  list(dadm = dadm, pars = EMC2:::get_pars_matrix_oo(p, dadm, model), ll = ll)
}

# Race likelihood from the kernel functions: winner density times loser survivors.
gm_ref_ll <- function(dadm, pars, dens, surv) {
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  u <- dadm$rt - pars[, "t0"]
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    w <- ix[win[ix]]
    log(dens(u[w], pars[w, ])) +
      sum(vapply(ix[!win[ix]], function(i) log(surv(u[i], pars[i, ])), 0))
  }, 0)
  sum(pmax(ll, log(1e-10)))
}

test_that("the compiled race likelihood uses the mixture kernels", {
  skip_on_cran()
  fx <- gm_fixture(function() OPR(parameterization = "ratio"),
                   list(K ~ 1, v ~ lM, rho ~ 1, t0 ~ 1, cv ~ lM),
                   c(K = log(1.2), v = log(5), v_lMTRUE = .8, rho = log(.6),
                     t0 = log(.15), cv = log(.4), cv_lMTRUE = -.3))
  expect_true(all(fx$pars[, "cv"] > 0))
  ref <- gm_ref_ll(fx$dadm, fx$pars,
    function(u, p) EMC2:::dopr_gamma(u, p["K"], p["v"], p["rho"], p["cv"]),
    function(u, p) EMC2:::popr_gamma(u, p["K"], p["v"], p["rho"], p["cv"],
                                     lower_tail = FALSE))
  expect_equal(fx$ll, ref, tolerance = 1e-8)
  fx <- gm_fixture(function() OPR(parameterization = "background"),
                   list(K ~ 1, v ~ lM, lambda0 ~ 1, t0 ~ 1, cv ~ lM),
                   c(K = log(1.4), v = log(4), v_lMTRUE = .9, lambda0 = log(7),
                     t0 = log(.15), cv = log(.4), cv_lMTRUE = -.3))
  expect_true(all(fx$pars[, "cv"] > 0))
  ref <- gm_ref_ll(fx$dadm, fx$pars,
    function(u, p) EMC2:::dopr_gamma(u, p["K"], p["v"], p["lambda0"] / p["v"], p["cv"]),
    function(u, p) EMC2:::popr_gamma(u, p["K"], p["v"], p["lambda0"] / p["v"], p["cv"],
                                     lower_tail = FALSE))
  expect_equal(fx$ll, ref, tolerance = 1e-8)
  for (packets in c(FALSE, TRUE)) {
    fx <- gm_fixture(if (packets) function() RSUB("packets") else RSUB,
                     list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, cv ~ 1),
                     c(v = log(1.5), v_lMTRUE = .9, B = log(.8), A = log(.3),
                       t0 = log(.15), s = log(.7), cv = log(.5)))
    ref <- gm_ref_ll(fx$dadm, fx$pars,
      function(u, p) EMC2:::drsub(u, p["v"], p["B"], p["A"], p["s"], packets, cv = p["cv"]),
      function(u, p) EMC2:::prsub(u, p["v"], p["B"], p["A"], p["s"], packets,
                                  lower_tail = FALSE, cv = p["cv"]))
    expect_equal(fx$ll, ref, tolerance = 1e-8, info = packets)
  }
})
