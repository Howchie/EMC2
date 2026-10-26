skip_model_validation()

# BAwL with a gamma launch strength (drift_distribution = "gamma"), sampled as
# the arithmetic mean and SD of V: shape (mean/sd)^2, scale sd^2/mean.  As in
# test-bawl-lognormal.R, every reference integrates over the START POINT, so
# it shares no algebra with the w-space kernels.

bawl_W <- function(a, t, b, k) {
  if (k <= 1e-12) return((b - a) / t)
  E <- exp(-k * t)
  k * (b - a * E) / (1 - E)
}

bawl_dW <- function(a, t, b, k) {
  if (k <= 1e-12) return((b - a) / t^2)
  E <- exp(-k * t)
  k^2 * E * (b - a) / (1 - E)^2
}

gam_shape <- function(mean, sd) (mean / sd)^2
gam_scale <- function(mean, sd) sd^2 / mean

ref_F_gam <- function(t, A, b, mean, sd, k) {
  a <- gam_shape(mean, sd); s <- gam_scale(mean, sd)
  if (!isTRUE(t > 0)) return(0)
  if (is.infinite(t))
    return(if (k <= 1e-12) 1 else pgamma(k * b, a, scale = s, lower.tail = FALSE))
  if (A <= 0) return(pgamma(bawl_W(0, t, b, k), a, scale = s, lower.tail = FALSE))
  integrate(function(x) pgamma(bawl_W(x, t, b, k), a, scale = s, lower.tail = FALSE),
            0, A, rel.tol = .Machine$double.eps^0.75)$value / A
}

ref_f_gam <- function(t, A, b, mean, sd, k) {
  a <- gam_shape(mean, sd); s <- gam_scale(mean, sd)
  if (!isTRUE(t > 0) || is.infinite(t)) return(0)
  if (A <= 0)
    return(dgamma(bawl_W(0, t, b, k), a, scale = s) * bawl_dW(0, t, b, k))
  integrate(function(x) dgamma(bawl_W(x, t, b, k), a, scale = s) *
              vapply(x, bawl_dW, numeric(1), t = t, b = b, k = k),
            0, A, rel.tol = .Machine$double.eps^0.75)$value / A
}

cpp_p_gam <- function(t, A, b, mean, sd, k, log_out = FALSE)
  EMC2:::pleakyba_norm(t, A, b, mean, sd, k, TRUE, log_out, 4L)
cpp_d_gam <- function(t, A, b, mean, sd, k)
  EMC2:::dleakyba_norm(t, A, b, mean, sd, k, TRUE, FALSE, 4L)

# Shapes 0.25 to 2500: the last lies outside the incomplete-gamma band and
# exercises the R::pgamma fallbacks.
gam_grid <- data.frame(
  A    = c(0,   0.3, 1,   0.5, 0.4, 0.8),
  b    = c(1,   2.5, 2.5, 1.5, 2,   1.8),
  mean = c(1.5, 2,   3,   1.2, 2.5, 2),
  sd   = c(0.8, 4,   1,   0.3, 2.5, 0.04),
  k    = c(0,   0.5, 3,   0,   1.2, 0.4)
)

# ---------------------------------------------------------------------------
# 1. Kernels
# ---------------------------------------------------------------------------

test_that("the gamma BAwL CDF and PDF match start-point quadrature", {
  for (i in seq_len(nrow(gam_grid))) {
    g <- gam_grid[i, ]
    for (t in c(0.2, 0.7, 1.5)) {
      want_F <- ref_F_gam(t, g$A, g$b, g$mean, g$sd, g$k)
      want_f <- ref_f_gam(t, g$A, g$b, g$mean, g$sd, g$k)
      expect_equal(cpp_p_gam(t, g$A, g$b, g$mean, g$sd, g$k), want_F,
                   tolerance = 1e-7)
      expect_equal(cpp_d_gam(t, g$A, g$b, g$mean, g$sd, g$k), want_f,
                   tolerance = 1e-6)
    }
  }
})

test_that("the gamma BAwL density integrates to its CDF", {
  for (g in list(c(A = 0.5, b = 2, mean = 2, sd = 1, k = 0),
                 c(A = 0.5, b = 2, mean = 2, sd = 3, k = 1.5),
                 c(A = 0, b = 1.5, mean = 2.5, sd = 1.2, k = 0.7),
                 c(A = 1.2, b = 1.2, mean = 1.5, sd = 0.5, k = 0.3))) {
    got <- integrate(function(u)
      vapply(u, cpp_d_gam, numeric(1), A = g[["A"]], b = g[["b"]],
             mean = g[["mean"]], sd = g[["sd"]], k = g[["k"]]),
      0, 4, rel.tol = 1e-10)$value
    expect_equal(got, cpp_p_gam(4, g[["A"]], g[["b"]], g[["mean"]], g[["sd"]], g[["k"]]),
                 tolerance = 1e-7)
  }
})

test_that("the gamma BAwL upper tail is defective exactly at P(V > k b)", {
  for (g in list(c(A = 0.5, b = 2, mean = 2, sd = 1, k = 1.5),
                 c(A = 0, b = 1, mean = 3, sd = 2, k = 2))) {
    expect_equal(cpp_p_gam(Inf, g[["A"]], g[["b"]], g[["mean"]], g[["sd"]], g[["k"]]),
                 pgamma(g[["k"]] * g[["b"]], gam_shape(g[["mean"]], g[["sd"]]),
                        scale = gam_scale(g[["mean"]], g[["sd"]]),
                        lower.tail = FALSE),
                 tolerance = 1e-12)
  }
  expect_equal(cpp_p_gam(Inf, 0.5, 2, 2, 1, 0), 1)
  expect_equal(cpp_d_gam(Inf, 0.5, 2, 2, 1, 1.5), 0)
})

test_that("the gamma BAwL log CDF keeps its digits deep in the early tail", {
  # Very short times need very large launch strengths, so log F is far below
  # anything the natural-scale path can represent.
  for (g in list(c(A = 0.5, b = 2, mean = 2, sd = 1, k = 0.5),
                 c(A = 0.3, b = 1.5, mean = 1, sd = 2, k = 0),
                 c(A = 0.4, b = 2, mean = 2, sd = 0.04, k = 0.4))) {
    for (t in c(0.05, 0.2)) {
      a <- gam_shape(g[["mean"]], g[["sd"]]); s <- gam_scale(g[["mean"]], g[["sd"]])
      l0 <- pgamma(bawl_W(g[["A"]], t, g[["b"]], g[["k"]]), a, scale = s,
                   lower.tail = FALSE, log.p = TRUE)
      wantF <- log(integrate(function(x)
        exp(pgamma(bawl_W(x, t, g[["b"]], g[["k"]]), a, scale = s,
                   lower.tail = FALSE, log.p = TRUE) - l0),
        0, g[["A"]], rel.tol = 1e-12)$value / g[["A"]]) + l0
      got <- cpp_p_gam(t, g[["A"]], g[["b"]], g[["mean"]], g[["sd"]], g[["k"]],
                       log_out = TRUE)
      expect_true(is.finite(got))
      expect_equal(got, wantF, tolerance = 1e-6)
    }
  }
})

test_that("the gamma BAwL kernels stay finite in the far tails", {
  tt <- c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100)
  for (sd in c(0.05, 1, 6)) {
    d <- vapply(tt, cpp_d_gam, numeric(1), A = 0.5, b = 2, mean = 2, sd = sd, k = 1.2)
    p <- vapply(tt, cpp_p_gam, numeric(1), A = 0.5, b = 2, mean = 2, sd = sd, k = 1.2)
    expect_true(all(is.finite(d)) && all(d >= 0))
    expect_true(all(is.finite(p)) && all(p >= 0 & p <= 1))
    expect_false(is.unsorted(p))
    expect_true(all(is.finite(vapply(tt, cpp_d_gam, numeric(1), A = 0, b = 2,
                                     mean = 2, sd = sd, k = 1.2))))
    expect_true(all(is.finite(vapply(tt, cpp_d_gam, numeric(1), A = 2, b = 2,
                                     mean = 2, sd = sd, k = 1.2))))
  }
})

test_that("a high-shape gamma launch approaches the positive-normal launch", {
  # With mean and SD held fixed the gamma tends to the normal as the shape
  # grows, so the two launches must agree to O(CV).
  for (t in c(0.6, 1, 1.6)) {
    pg <- cpp_p_gam(t, 0.4, 1.5, 2, 0.1, 0.3)
    pn <- EMC2:::pleakyba_norm(t, 0.4, 1.5, 2, 0.1, 0.3, TRUE, FALSE, 0L)
    expect_lt(abs(pg - pn), 0.02)
  }
})

# ---------------------------------------------------------------------------
# 2. Model construction
# ---------------------------------------------------------------------------

test_that("BAwL(drift_distribution = \"gamma\") selects the launch pair and c_name", {
  m <- BAwL(drift_distribution = "gamma")
  expect_equal(m$c_name, "BAwL_GAMMA")
  expect_equal(m$drift_distribution, "gamma")
  expect_equal(m$p_types_canonical, c("mean", "sd", "B", "A", "t0", "k"))
  expect_equal(unname(m$transform$func[c("mean", "sd")]), c("exp", "exp"))
  expect_equal(BAwL(drift_distribution = "gamma", erlang_type = "local_kill")$c_name,
               "BAwL_GAMMA_LOCAL_KILL")
  expect_equal(BAwL(drift_distribution = "gamma", erlang_shape = 2L)$c_name,
               "BAwL_GAMMA_E2")
  expect_equal(BAwL(drift_distribution = "gamma", parameterization = "ratio")$c_name,
               "BAwLRAT_GAMMA")
  expect_error(BAwL(drift_distribution = "gamma", posdrift = FALSE),
               "positive by construction")
  expect_error(BAwL(drift_distribution = "gamma", correlated = TRUE),
               "only implemented for")
})

test_that("the gamma Ttransform emits the kernel column order", {
  m <- BAwL(drift_distribution = "gamma")
  pars <- cbind(mean = 1.4, sd = 0.5, B = 1, A = 0.3, t0 = 0.2, k = 0.8,
                mG = 1, mK = 1)
  out <- m$Ttransform(pars, NULL)
  expect_equal(colnames(out)[1:8],
               c("mean", "sd", "B", "A", "t0", "k", "lambda_g", "lambda_k"))
  r <- BAwL(drift_distribution = "gamma", parameterization = "ratio")
  pars_r <- pars; colnames(pars_r)[6] <- "r"
  expect_equal(unname(r$Ttransform(pars_r, NULL)[, "k"]), 0.8 * 1.4)
})

test_that("the model's dfun and pfun call the gamma kernels", {
  m <- BAwL(drift_distribution = "gamma")
  pars <- cbind(mean = 2, sd = 1.3, B = 1, A = 0.3, t0 = 0.2, k = 0.8,
                lambda_g = 0, lambda_k = 0, b = 1.3)
  rt <- 0.9
  expect_equal(m$dfun(rt, pars), ref_f_gam(rt - 0.2, 0.3, 1.3, 2, 1.3, 0.8),
               tolerance = 1e-6)
  expect_equal(m$pfun(rt, pars), ref_F_gam(rt - 0.2, 0.3, 1.3, 2, 1.3, 0.8),
               tolerance = 1e-7)
  bad <- pars
  colnames(bad)[1:2] <- c("v", "sv")
  expect_error(m$dfun(rt, bad), "requires parameter columns")
})

# ---------------------------------------------------------------------------
# 3. Simulation
# ---------------------------------------------------------------------------

test_that("the gamma launch draws have the sampled mean and SD", {
  set.seed(20261026)
  x <- EMC2:::.rgamma_meansd(200000L, 2.5, 0.8)
  expect_equal(mean(x), 2.5, tolerance = 0.01)
  expect_equal(sd(x), 0.8, tolerance = 0.02)
})

test_that("both gamma BAwL simulators match the CDF", {
  skip_on_cran()
  # One accumulator per trial, so the simulated hit time is the racer itself.
  n <- 20000L
  lR <- factor(rep("a", n), levels = "a")
  pars <- cbind(mean = 2, sd = 1.2, b = 1.3, A = 0.4, t0 = 0.1, k = 0.5,
                lambda_g = 0, lambda_k = 0)
  pars <- pars[rep(1, n), , drop = FALSE]
  for (cpp in c(TRUE, FALSE)) {
    withr::local_options(emc2.cpp_rfun = cpp)
    set.seed(31)
    sim <- EMC2:::.rfun_BAwL(lR, pars, posdrift = TRUE, launch = 4L)
    rt <- sim$rt
    for (q in c(0.4, 0.8, 1.5)) {
      emp <- mean(!is.na(rt) & is.finite(rt) & rt <= q)
      want <- ref_F_gam(q - 0.1, 0.4, 1.3, 2, 1.2, 0.5)
      expect_lt(abs(emp - want), 4 * sqrt(want * (1 - want) / n))
    }
  }
})

# ---------------------------------------------------------------------------
# 4. The compiled likelihood
# ---------------------------------------------------------------------------

bawl_gam_ref_ll <- function(dadm, pars, min_ll = log(1e-10)) {
  Fu <- function(i, u) ref_F_gam(u, pars[i, "A"], pars[i, "b"], pars[i, "mean"],
                                 pars[i, "sd"], pars[i, "k"])
  fu <- function(i, u) ref_f_gam(u, pars[i, "A"], pars[i, "b"], pars[i, "mean"],
                                 pars[i, "sd"], pars[i, "k"])
  n_lR <- length(levels(dadm$lR))
  lt <- numeric(nrow(dadm) / n_lR)
  for (j in seq_along(lt)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    w <- idx[which(dadm$winner[idx])]
    term <- fu(w, rt - pars[w, "t0"])
    for (i in setdiff(idx, w)) term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

bawl_gam_fixture <- function(form, consts) {
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(60), ]
  des <- design(data = dat, model = function() BAwL(drift_distribution = "gamma"),
                matchfun = matchfun, formula = form, constants = consts)
  list(dat = dat, des = des,
       emc = suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                       compress = FALSE, rt_resolution = NULL)))
}

bawl_gam_ll <- function(fx, p) {
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

test_that("the compiled gamma BAwL likelihood matches the R reference", {
  skip_on_cran()
  fx <- suppressMessages(bawl_gam_fixture(
    list(mean ~ lM, sd ~ lM, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    c(mG = log(1), mK = log(1))))
  p <- c(mean = log(1.5), mean_lMTRUE = 0.5, sd = log(0.9), sd_lMTRUE = 0.2,
         B = log(0.8), A = log(0.3), t0 = log(0.15), k = log(0.8))
  p <- p[names(sampled_pars(fx$des))]
  got <- bawl_gam_ll(fx, p)
  expect_true(is.finite(got))
  pars <- EMC2:::get_pars_matrix_oo(p, fx$emc[[1]]$data[[1]],
                                    fx$emc[[1]]$model())
  expect_equal(got, bawl_gam_ref_ll(fx$emc[[1]]$data[[1]], pars),
               tolerance = 1e-6)
})

test_that("the gamma launch scales with (V, b, A) and every parameter matters", {
  skip_on_cran()
  # (V, b, A) -> (cV, cb, cA) leaves every crossing time unchanged, so sd
  # carries the evidence scale as sv does for the normal launch.  Checked on
  # the natural scale: the sampled scale is not an exact rescaling once the
  # bound link puts the lower bound into the transform.
  for (cc in c(0.4, 2.5)) {
    for (t in c(0.3, 0.8, 2)) {
      expect_equal(cpp_d_gam(t, cc * 0.3, cc * 1.3, cc * 2, cc * 1.1, 0.6),
                   cpp_d_gam(t, 0.3, 1.3, 2, 1.1, 0.6), tolerance = 1e-10)
      expect_equal(cpp_p_gam(t, cc * 0.3, cc * 1.3, cc * 2, cc * 1.1, 0.6),
                   cpp_p_gam(t, 0.3, 1.3, 2, 1.1, 0.6), tolerance = 1e-10)
    }
  }
  fx <- suppressMessages(bawl_gam_fixture(
    list(mean ~ lM, sd ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    c(mG = log(1), mK = log(1))))
  p <- c(mean = log(1.5), mean_lMTRUE = 0.5, sd = log(0.9), B = log(0.8),
         A = log(0.3), t0 = log(0.15), k = log(0.8))
  p <- p[names(sampled_pars(fx$des))]
  base <- bawl_gam_ll(fx, p)
  for (nm in names(p)) {
    q <- p; q[nm] <- q[nm] + 0.1
    expect_false(isTRUE(all.equal(base, bawl_gam_ll(fx, q))))
  }
})

test_that("the gamma launch composes with the kill and guess clocks", {
  skip_on_cran()
  mk <- BAwL(drift_distribution = "gamma", erlang_type = "local_kill")
  mg <- BAwL(drift_distribution = "gamma", erlang_type = "local_guess")
  base <- cbind(mean = 2, sd = 1, B = 1, A = 0.3, t0 = 0.2, k = 0.6,
                mG = 1, mK = 1)
  pk <- mk$Ttransform(base, NULL)
  pg <- mg$Ttransform(base, NULL)
  p0 <- pk; p0[, "lambda_k"] <- 0
  expect_lt(mk$pfun(Inf, pk), mk$pfun(Inf, p0))
  expect_gt(mg$pfun(0.5, pg), mg$pfun(0.5, p0))
})
