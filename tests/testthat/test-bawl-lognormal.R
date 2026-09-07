skip_model_validation()

# BAwL with a lognormal launch strength (drift_distribution = "lognormal").
#
# The crossing condition V >= W(a) = k (b - a e^{-kt}) / (1 - e^{-kt}) is affine
# in the start point a ~ U(0, A), so substituting w = W(a) turns both integrals
# over (V, a) into integrals over w alone and the lognormal launch stays closed
# form.  Every reference below integrates over the START POINT instead, so it
# shares no algebra with the kernels: it tests the substitution, not just its
# arithmetic.

# ---------------------------------------------------------------------------
# References: quadrature over the start point
# ---------------------------------------------------------------------------

bawl_W <- function(a, t, b, k) {
  if (k <= 1e-12) return((b - a) / t)
  E <- exp(-k * t)
  k * (b - a * E) / (1 - E)
}

# -dW/dt, the Jacobian that turns P(V >= W) into a density.
bawl_dW <- function(a, t, b, k) {
  if (k <= 1e-12) return((b - a) / t^2)
  E <- exp(-k * t)
  k^2 * E * (b - a) / (1 - E)^2
}

ref_F_logn <- function(t, A, b, mu, sigma, k) {
  if (!isTRUE(t > 0)) return(0)
  if (is.infinite(t))
    return(if (k <= 1e-12) 1 else
      plnorm(k * b, mu, sigma, lower.tail = FALSE))
  if (A <= 0) return(plnorm(bawl_W(0, t, b, k), mu, sigma, lower.tail = FALSE))
  integrate(function(a) plnorm(bawl_W(a, t, b, k), mu, sigma, lower.tail = FALSE),
            0, A, rel.tol = .Machine$double.eps^0.75)$value / A
}

ref_f_logn <- function(t, A, b, mu, sigma, k) {
  if (!isTRUE(t > 0) || is.infinite(t)) return(0)
  if (A <= 0)
    return(dlnorm(bawl_W(0, t, b, k), mu, sigma) * bawl_dW(0, t, b, k))
  integrate(function(a) dlnorm(bawl_W(a, t, b, k), mu, sigma) *
              vapply(a, bawl_dW, numeric(1), t = t, b = b, k = k),
            0, A, rel.tol = .Machine$double.eps^0.75)$value / A
}

# launch = 1L is sampled on the natural scale as (mean, cv); the reference
# integrals below stay in (mu, sigma), so convert at the call boundary.
ln_meancv <- function(mu, sigma)
  c(mean = unname(exp(mu + sigma^2 / 2)), cv = unname(sqrt(expm1(sigma^2))))
ln_musigma <- function(mean, cv) {
  s2 <- log1p(unname(cv)^2)
  c(mu = unname(log(mean) - s2 / 2), sigma = unname(sqrt(s2)))
}
cpp_p_logn <- function(t, A, b, mu, sigma, k) {
  mc <- ln_meancv(mu, sigma)
  EMC2:::pleakyba(t, A, b, mc[["mean"]], mc[["cv"]], k, TRUE, 1L)
}
cpp_d_logn <- function(t, A, b, mu, sigma, k) {
  mc <- ln_meancv(mu, sigma)
  EMC2:::dleakyba(t, A, b, mc[["mean"]], mc[["cv"]], k, TRUE, 1L)
}

# ---------------------------------------------------------------------------
# 1. Kernels
# ---------------------------------------------------------------------------

test_that("the lognormal BAwL CDF and PDF match start-point quadrature", {
  grid <- data.frame(
    A = c(0, 0.3, 1),
    b = c(1, 2.5, 2.5),
    mu = c(-0.5, 0.2, 1),
    sigma = c(0.3, 0.9, 0.5),
    k = c(0, 0.5, 3)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    for (t in c(0.2, 1)) {
      want_F <- ref_F_logn(t, g$A, g$b, g$mu, g$sigma, g$k)
      want_f <- ref_f_logn(t, g$A, g$b, g$mu, g$sigma, g$k)
      expect_equal(cpp_p_logn(t, g$A, g$b, g$mu, g$sigma, g$k), want_F,
                   tolerance = 1e-7)
      expect_equal(cpp_d_logn(t, g$A, g$b, g$mu, g$sigma, g$k), want_f,
                   tolerance = 1e-6)
    }
  }
})

test_that("the lognormal BAwL density integrates to its CDF", {
  for (g in list(c(A = 0.5, b = 2, mu = 0.3, sigma = 0.5, k = 0),
                 c(A = 0.5, b = 2, mu = 0.3, sigma = 0.5, k = 1.5),
                 c(A = 0, b = 1.5, mu = 0.8, sigma = 0.6, k = 0.7),
                 c(A = 1.2, b = 1.2, mu = -0.2, sigma = 1.1, k = 2.5))) {
    got <- integrate(function(u)
      cpp_d_logn(u, g["A"], g["b"], g["mu"], g["sigma"], g["k"]),
      0, 4, rel.tol = 1e-10)$value
    expect_equal(got, cpp_p_logn(4, g["A"], g["b"], g["mu"], g["sigma"], g["k"]),
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("the lognormal BAwL upper tail is defective exactly at P(V > k b)", {
  # With leak, only launch strengths above k b ever reach threshold, whatever
  # the start point; k = 0 is proper because the launch strength is positive.
  for (g in list(c(A = 0.5, b = 2, mu = 0.3, sigma = 0.5, k = 1.5),
                 c(A = 0, b = 1, mu = 1, sigma = 0.4, k = 2))) {
    expect_equal(cpp_p_logn(Inf, g["A"], g["b"], g["mu"], g["sigma"], g["k"]),
                 plnorm(g["k"] * g["b"], g["mu"], g["sigma"], lower.tail = FALSE),
                 tolerance = 1e-12, ignore_attr = TRUE)
  }
  expect_equal(cpp_p_logn(Inf, 0.5, 2, 0.3, 0.5, 0), 1)
  expect_equal(cpp_d_logn(Inf, 0.5, 2, 0.3, 0.5, 1.5), 0)
})

test_that("the lognormal BAwL kernels stay finite in the far tails", {
  tt <- c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100)
  d <- cpp_d_logn(tt, 0.5, 2, 0.3, 0.5, 1.2)
  p <- cpp_p_logn(tt, 0.5, 2, 0.3, 0.5, 1.2)
  expect_true(all(is.finite(d)) && all(d >= 0))
  expect_true(all(is.finite(p)) && all(p >= 0 & p <= 1))
  expect_false(is.unsorted(p))                       # CDF is monotone
  expect_true(all(is.finite(cpp_d_logn(tt, 0, 2, 0.3, 0.5, 1.2))))
  # A = b is the extreme start-point range (the accumulator can start at
  # threshold); it must not produce a negative density or NaN.
  expect_true(all(is.finite(cpp_d_logn(tt, 2, 2, 0.3, 0.5, 1.2))))
})

test_that("the normal-launch BAwL kernels are untouched", {
  # Regression guard: adding the launch switch must not perturb the Gaussian
  # branch, which is what the existing BAwL/LBA fits use.
  expect_equal(EMC2:::pleakyba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8, 0.4, TRUE, 0L),
               EMC2:::pleakyba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8, 0.4, TRUE))
  expect_equal(EMC2:::dleakyba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8, 0.4, TRUE, 0L),
               EMC2:::dleakyba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8, 0.4, TRUE))
  expect_equal(EMC2:::plba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8),
               EMC2:::pleakyba(c(0.3, 1, 2), 0.5, 2, 1.2, 0.8, 0, TRUE))
})

# ---------------------------------------------------------------------------
# 2. Model construction
# ---------------------------------------------------------------------------

test_that("BAwL(drift_distribution) selects the launch pair and c_name", {
  ln <- BAwL(drift_distribution = "lognormal")
  no <- BAwL()
  expect_equal(ln$c_name, "BAwL_LOGN")
  expect_equal(no$c_name, "BAwL")
  expect_equal(no$drift_distribution, "normal")     # default is unchanged
  expect_equal(names(ln$p_types)[1:2], c("mean", "cv"))
  expect_equal(names(no$p_types)[1:2], c("v", "sv"))
  expect_equal(ln$p_types_canonical, c("mean", "cv", "B", "A", "t0", "k"))
  expect_equal(unname(ln$transform$func[c("mean", "cv")]),
               c("exp", "exp"))
  # The clock suffixes still compose, and _LOGN never meets IO.
  expect_equal(BAwL(drift_distribution = "lognormal",
                    erlang_type = "local_kill")$c_name,
               "BAwL_LOGN_LOCAL_KILL")
  expect_equal(BAwL(drift_distribution = "lognormal", erlang_shape = 2L)$c_name,
               "BAwL_LOGN_E2")
})

test_that("lognormal BAwL refuses the combinations it does not implement", {
  # posdrift is meaningless (V > 0 by construction) and the correlated path is
  # a one-factor decomposition of the *Gaussian* drift vector.
  expect_error(BAwL(drift_distribution = "lognormal", posdrift = FALSE),
               "positive by construction")
  expect_error(BAwL(drift_distribution = "lognormal", correlated = TRUE),
               "only implemented for")
  expect_error(BAwLcorr(), NA)
})

test_that("the lognormal Ttransform emits the kernel column order", {
  # p_types prefix IS the kernel column order, so mean/cv must lead.
  m <- BAwL(drift_distribution = "lognormal")
  pars <- cbind(mean = 1.4, cv = 0.5, B = 1, A = 0.3, t0 = 0.2, k = 0.8,
                mG = 1, mK = 1)
  out <- m$Ttransform(pars, NULL)
  expect_equal(colnames(out)[1:8],
               c("mean", "cv", "B", "A", "t0", "k", "lambda_g", "lambda_k"))
  expect_equal(unname(out[, "b"]), 1.3)
})

test_that("the model's dfun and pfun call the lognormal kernels", {
  m <- BAwL(drift_distribution = "lognormal")
  mc <- ln_meancv(0.3, 0.5)
  pars <- cbind(mean = mc[["mean"]], cv = mc[["cv"]], B = 1, A = 0.3, t0 = 0.2,
                k = 0.8, lambda_g = 0, lambda_k = 0, b = 1.3)
  rt <- 0.9
  expect_equal(m$dfun(rt, pars), cpp_d_logn(rt - 0.2, 0.3, 1.3, 0.3, 0.5, 0.8),
               ignore_attr = TRUE)
  expect_equal(m$pfun(rt, pars), cpp_p_logn(rt - 0.2, 0.3, 1.3, 0.3, 0.5, 0.8),
               ignore_attr = TRUE)
  # A lognormal design must not silently accept v/sv columns.
  bad <- pars
  colnames(bad)[1:2] <- c("v", "sv")
  expect_error(m$dfun(rt, bad), "requires parameter columns")
})

# ---------------------------------------------------------------------------
# 3. Simulation
# ---------------------------------------------------------------------------

test_that("the lognormal BAwL simulators produce valid draws", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), 100), levels = c("left", "right"))
  mc <- ln_meancv(0.4, 0.5)
  pars <- cbind(mean = mc[["mean"]], cv = mc[["cv"]], b = 1.3, A = 0.3,
                t0 = 0.2, k = 0.6, lambda_g = 0, lambda_k = 0)
  pars <- pars[rep(1, length(lR)), , drop = FALSE]
  for (cpp in c(TRUE, FALSE)) {
    withr::local_options(emc2.cpp_rfun = cpp)
    set.seed(11)
    sim <- EMC2:::.rfun_BAwL(lR, pars, posdrift = TRUE, launch = 1L)
    expect_equal(nrow(sim), 100)
    expect_true(all(sim$R %in% c("left", "right", NA)))
    expect_true(all(sim$rt[!is.na(sim$rt)] >= 0.2))
  }
})

# ---------------------------------------------------------------------------
# 4. The compiled likelihood
# ---------------------------------------------------------------------------

# An independent R race likelihood built from the start-point references above:
# nothing is shared with the compiled path except the parameter values, so this
# covers design, Ttransform, column order, adapter dispatch and race assembly.
bawl_logn_ref_ll <- function(dadm, pars, min_ll = log(1e-10)) {
  Fu <- function(i, u) {
    q <- ln_musigma(pars[i, "mean"], pars[i, "cv"])
    ref_F_logn(u, pars[i, "A"], pars[i, "b"], q[["mu"]], q[["sigma"]], pars[i, "k"])
  }
  fu <- function(i, u) {
    q <- ln_musigma(pars[i, "mean"], pars[i, "cv"])
    ref_f_logn(u, pars[i, "A"], pars[i, "b"], q[["mu"]], q[["sigma"]], pars[i, "k"])
  }
  n_lR <- length(levels(dadm$lR))
  lt <- numeric(nrow(dadm) / n_lR)
  for (j in seq_along(lt)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    if (is.infinite(rt) && rt > 0) {
      term <- prod(vapply(idx, function(i) 1 - Fu(i, Inf), numeric(1)))
    } else {
      w <- idx[which(dadm$winner[idx])]
      term <- fu(w, rt - pars[w, "t0"])
      for (i in setdiff(idx, w)) term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    }
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

bawl_logn_fixture <- function(form, consts, model = NULL) {
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  if (is.null(model))
    model <- function() BAwL(drift_distribution = "lognormal")
  des <- design(data = dat, model = model, matchfun = matchfun,
                formula = form, constants = consts)
  list(dat = dat, des = des,
       emc = suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                       compress = FALSE, rt_resolution = NULL)))
}

bawl_logn_ll <- function(fx, p) {
  model <- fx$emc[[1]]$model()
  dadm <- fx$emc[[1]]$data[[1]]
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform,
                    p_types = names(model$p_types), min_ll = log(1e-10),
                    trend = model$trend)
}

test_that("the compiled lognormal BAwL likelihood matches the R reference", {
  skip_on_cran()
  fx <- suppressMessages(bawl_logn_fixture(
    list(mean ~ lM, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    c(mG = log(1), mK = log(1))))
  p <- c(mean = 0.9, mean_lMTRUE = 0.4, cv = log(0.6), B = log(0.8),
         A = log(0.3), t0 = log(0.15), k = log(0.8))
  p <- p[names(sampled_pars(fx$des))]
  got <- bawl_logn_ll(fx, p)
  expect_true(is.finite(got))
  pars <- EMC2:::get_pars_matrix_oo(p, fx$emc[[1]]$data[[1]],
                                    fx$emc[[1]]$model())
  expect_equal(got, bawl_logn_ref_ll(fx$emc[[1]]$data[[1]], pars),
               tolerance = 1e-6)
})

test_that("the lognormal likelihood is sensitive to every free parameter", {
  # Guards against a column-order or dispatch slip that silently ignores one of
  # the launch parameters (the failure mode a plain finite-value check misses).
  skip_on_cran()
  fx <- suppressMessages(bawl_logn_fixture(
    list(mean ~ 1, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    c(mG = log(1), mK = log(1))))
  p <- c(mean = 0.9, cv = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(0.8))
  base <- bawl_logn_ll(fx, p)
  for (nm in names(p)) {
    q <- p; q[nm] <- q[nm] + 0.1
    expect_false(isTRUE(all.equal(base, bawl_logn_ll(fx, q))))
  }
})

test_that("data simulated from the lognormal model round-trips through the likelihood", {
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = function() BAwL(drift_distribution = "lognormal"),
                matchfun = matchfun,
                formula = list(mean ~ lM, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
                constants = c(mean = 0))
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- c(mean_lMTRUE = 0.8, cv = log(0.5), B = log(0.9), A = log(0.3),
           t0 = log(0.2), k = log(0.5))[names(p)]
  set.seed(7)
  sim <- suppressMessages(make_data(p, design = des, n_trials = 100))
  expect_true(all(is.finite(sim$rt) | sim$rt == Inf))
  emc <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  ll <- EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
                          dadm, constants = attr(dadm, "constants"),
                          designs = designs, type = model$c_name,
                          bounds = model$bound, transforms = model$transform,
                          pretransforms = model$pre_transform,
                          p_types = names(model$p_types), min_ll = log(1e-10),
                          trend = model$trend)
  expect_true(is.finite(ll))
  # The data-generating values must beat a clearly wrong launch location.
  q <- p; q["mean_lMTRUE"] <- 0
  ll_wrong <- EMC2:::calc_ll_oo(matrix(q, nrow = 1, dimnames = list(NULL, names(q))),
                                dadm, constants = attr(dadm, "constants"),
                                designs = designs, type = model$c_name,
                                bounds = model$bound, transforms = model$transform,
                                pretransforms = model$pre_transform,
                                p_types = names(model$p_types),
                                min_ll = log(1e-10), trend = model$trend)
  expect_gt(ll, ll_wrong)
})

test_that("the lognormal evidence scale is fixed by mu, B or A but not sigma", {
  skip_on_cran()
  # (V, b, A) -> (cV, cb, cA) leaves every crossing time unchanged.  On the
  # sampled scale that is the single direction mu += log c, B += log c,
  # A += log c; sigma is dimensionless and does NOT identify it, which is why
  # the normal model's constants = c(sv = log(1)) has no lognormal analogue.
  fx <- suppressMessages(bawl_logn_fixture(
    list(mean ~ lM, cv ~ 1, B ~ E, A ~ 1, t0 ~ 1, k ~ 1),
    c(mG = log(1), mK = log(1))))
  p <- c(mean = 0.9, mean_lMTRUE = 0.4, cv = log(0.6), B = log(0.8),
         B_Eneutral = 0.1, B_Eaccuracy = 0.3, A = log(0.3), t0 = log(0.15),
         k = log(0.8))
  p <- p[names(sampled_pars(fx$des))]
  base <- bawl_logn_ll(fx, p)
  for (lc in c(-0.7, 0.35, 1.2)) {
    q <- p
    q[c("mean", "B", "A")] <- q[c("mean", "B", "A")] + lc
    expect_equal(bawl_logn_ll(fx, q), base, tolerance = 1e-10)
  }
  # Any single leg of that direction on its own must move the likelihood,
  # otherwise the constraint would not be removing anything.
  for (nm in c("mean", "B", "A")) {
    q <- p; q[nm] <- q[nm] + 0.35
    expect_false(isTRUE(all.equal(bawl_logn_ll(fx, q), base)))
  }
  # The rescaling shifts intercepts only, so contrast coefficients are
  # unaffected by it and stay identified whichever intercept is pinned.
  for (nm in c("mean_lMTRUE", "B_Eneutral", "B_Eaccuracy", "cv")) {
    q <- p; q[nm] <- q[nm] + 0.15
    expect_false(isTRUE(all.equal(bawl_logn_ll(fx, q), base)))
  }
})

test_that("the lognormal launch composes with the kill and guess clocks", {
  skip_on_cran()
  # The clocks wrap the base kernels, so this checks the wiring rather than new
  # mathematics: a kill clock must remove response mass, a guess clock add it.
  m0 <- BAwL(drift_distribution = "lognormal")
  mk <- BAwL(drift_distribution = "lognormal", erlang_type = "local_kill")
  mg <- BAwL(drift_distribution = "lognormal", erlang_type = "local_guess")
  mc <- ln_meancv(0.4, 0.5)
  pars <- cbind(mean = mc[["mean"]], cv = mc[["cv"]], B = 1, A = 0.3, t0 = 0.2, k = 0.6,
                lambda_g = 0, lambda_k = 0, b = 1.3)
  pars_k <- pars; pars_k[, "lambda_k"] <- 2
  pars_g <- pars; pars_g[, "lambda_g"] <- 2
  rt <- 0.8
  expect_lt(mk$pfun(rt, pars_k), m0$pfun(rt, pars))
  expect_gt(mg$pfun(rt, pars_g), m0$pfun(rt, pars))
  expect_equal(mk$pfun(rt, pars), m0$pfun(rt, pars))   # clocks off = base model
})
