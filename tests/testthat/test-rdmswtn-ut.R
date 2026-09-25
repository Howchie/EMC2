skip_model_validation()
library(EMC2)

# RDMSWTN_UT: the ordinary RDMSWTN accumulator on an accelerating clock,
# q(x) = x + u x^2/2 (linear) or expm1(u x)/u (exponential), x = t - t0.
# The compiled likelihood is the RDMSWTN kernels behind the generic
# operational-time warp, so most checks compare against RDMSWTN at q(x).

clocks <- c("linear", "exponential")
ut_q <- function(x, u, clock)
  if (clock == "linear") x + u * x^2 / 2 else expm1(u * x) / u
ut_theta <- function(x, u, clock)
  if (clock == "linear") 1 + u * x else exp(u * x)
ut_code <- c(linear = 1L, exponential = 2L)

ut_pars <- function(u = .8, posdrift = TRUE, v = 1.2, A = .2, sv = .3) {
  RDMSWTN_UT(posdrift = posdrift)$Ttransform(
    cbind(v = v, B = 1, A = A, t0 = .1, s = .9, sv = sv, u = u), NULL
  )
}
as_ordinary <- function(pars, x = NULL) {
  out <- if (is.null(x)) pars else pars[rep(1L, length(x)), , drop = FALSE]
  out[, "t0"] <- 0
  cbind(out, lambda_g = 0, lambda_k = 0)
}

# Compiled-likelihood context shared by the design-level tests.
ut_formula <- function(correlated = FALSE)
  c(list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, u ~ 1),
    if (correlated) list(rho ~ 1))
ut_context <- function(dat, model, formula = ut_formula(isTRUE(model$correlated)),
                       constants = c(s = 0), matchfun = NULL) {
  des <- design(data = dat, Rlevels = levels(dat$R), model = model,
                formula = formula, constants = constants, matchfun = matchfun,
                report_p_vector = FALSE)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = NULL)[[1]]
  list(design = des, emc = emc)
}
ut_values <- function(p, v = c(1.3, .9), u = .8, sv = .3, rho = NULL,
                      natural_v = FALSE) {
  p[grep("^v_", names(p))] <- if (natural_v) v else log(v)
  p["B"] <- 0; p["A"] <- log(.2); p["t0"] <- log(.1)
  # sv and u are sampled on the log scale: 0 is log(0), never a plain 0.
  p["sv"] <- log(sv)
  if ("u" %in% names(p)) p["u"] <- log(u)
  if (!is.null(rho)) p["rho"] <- qnorm((rho + 1) / 2)
  p
}
ut_ll <- function(ctx, p, pointwise = FALSE) {
  dadm <- ctx$emc$data[[1]]
  m <- ctx$emc$model()
  designs <- lapply(names(m$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(m$p_types)
  args <- list(
    particle_matrix = matrix(p, 1, dimnames = list(NULL, names(p))),
    data = dadm, constants = attr(dadm, "constants"), designs = designs,
    type = m$c_name, bounds = m$bound, transforms = m$transform,
    pretransforms = m$pre_transform, p_types = names(m$p_types),
    min_ll = log(1e-10), trend = m$trend)
  if (pointwise) as.numeric(do.call(EMC2:::calc_ll_oo_pw, args))
  else as.numeric(do.call(EMC2:::calc_ll_oo, c(args, list(marginalise = NULL))))
}
ut_mapped <- function(ctx, p)
  EMC2:::get_pars_matrix_oo(p, ctx$emc$data[[1]], ctx$emc$model())

test_that("constructors expose the urgency-clock contract", {
  for (clock in clocks) {
    model <- RDMSWTN_UT(clock = clock)
    expect_identical(model$p_types_canonical,
                     c("v", "B", "A", "t0", "s", "sv", "u"))
    expect_identical(names(model$p_types)[1:7], model$p_types_canonical)
    expect_identical(model$transform$func[["u"]], "exp")
    expect_identical(unname(model$p_types[["u"]]), log(0))
    expect_equal(unname(model$bound$minmax[, "u"]), c(0, Inf))
    expect_equal(unname(model$bound$exception[["u"]]), 0)
    expect_identical(model$clock, clock)
  }
  sfx <- c(linear = "", exponential = "_EXP")
  for (clock in clocks) {
    expect_identical(RDMSWTN_UT(clock = clock)$c_name,
                     paste0("RDMSWTN_UT", sfx[[clock]]))
    expect_identical(RDMSWTN_UT(posdrift = FALSE, clock = clock)$c_name,
                     paste0("RDMSWTN_UT_IO", sfx[[clock]]))
    expect_identical(RDMSWTN_UTcorr(clock = clock)$c_name,
                     paste0("RDMSWTN_UT", sfx[[clock]], "_CORR"))
    expect_identical(RDMSWTN_UTcorr(clock = clock, correlate = "drifts")$c_name,
                     paste0("RDMSWTN_UT", sfx[[clock]], "_CORRD"))
  }
  # The shared constructor leaves RDMSWTN_TT's contract untouched.
  tt <- RDMSWTN_TT()
  expect_identical(tt$c_name, "RDMSWTN_TT")
  expect_identical(tt$p_types_canonical, c("v", "B", "A", "t0", "s", "sv", "tau"))
  expect_equal(unname(tt$bound$minmax[, "tau"]), c(1e-4, Inf))
  expect_false("tau" %in% names(tt$bound$exception))
  expect_null(tt$clock)
})

test_that("both clocks are exact operational-time transforms of RDMSWTN", {
  x <- c(.05, .3, .8, 1.5, 3)
  for (clock in clocks) for (pd in c(TRUE, FALSE))
    for (cfg in list(list(), list(A = 0), list(sv = 0), list(A = 0, sv = 0))) {
      pars <- do.call(ut_pars, c(list(posdrift = pd), cfg))
      u <- pars[1, "u"]
      t <- x + pars[1, "t0"]
      ordinary <- as_ordinary(pars, x)
      q <- ut_q(x, u, clock)
      expect_equal(EMC2:::dRDMSWTN_UT(t, pars, posdrift = pd, clock = clock),
                   EMC2:::dRDMSWTN(q, ordinary, posdrift = pd) *
                     ut_theta(x, u, clock), tolerance = 1e-10)
      expect_equal(EMC2:::pRDMSWTN_UT(t, pars, posdrift = pd, clock = clock),
                   EMC2:::pRDMSWTN(q, ordinary, posdrift = pd),
                   tolerance = 1e-10)
      # The clock is unbounded, so the eventual mass is RDMSWTN's own:
      # 1 under posdrift, the intrinsic defect otherwise.
      F_inf <- EMC2:::pRDMSWTN_UT(Inf, pars, posdrift = pd, clock = clock)
      expect_equal(F_inf, EMC2:::pRDMSWTN(Inf, as_ordinary(pars),
                                          posdrift = pd), tolerance = 1e-12)
      if (pd) expect_equal(F_inf, 1, tolerance = 1e-12)
    }
})

test_that("the density is the derivative of the CDF and integrates to F(Inf)", {
  for (clock in clocks) for (pd in c(TRUE, FALSE)) {
    pars <- ut_pars(u = 1.5, posdrift = pd, v = if (pd) 1.2 else .3)
    t <- c(.2, .4, .7, 1.1)
    h <- 1e-5
    dF <- (EMC2:::pRDMSWTN_UT(t + h, pars, pd, clock = clock) -
             EMC2:::pRDMSWTN_UT(t - h, pars, pd, clock = clock)) / (2 * h)
    expect_equal(dF, EMC2:::dRDMSWTN_UT(t, pars, pd, clock = clock),
                 tolerance = 1e-6)
    mass <- integrate(function(z) EMC2:::dRDMSWTN_UT(z, pars, pd, clock = clock),
                      pars[1, "t0"], Inf, rel.tol = 1e-10,
                      subdivisions = 2000L)$value
    expect_equal(mass, EMC2:::pRDMSWTN_UT(Inf, pars, pd, clock = clock),
                 tolerance = 1e-7)
    if (!pd) expect_lt(mass, 1 - 1e-3)       # the IO defect is material
  }
})

test_that("u = 0 is exactly RDMSWTN, analytically and in the compiled race", {
  t <- c(.2, .5, 1.2, Inf)
  for (clock in clocks) for (pd in c(TRUE, FALSE)) {
    zero <- ut_pars(u = 0, posdrift = pd)
    ordinary <- cbind(zero, lambda_g = 0, lambda_k = 0)
    expect_identical(EMC2:::dRDMSWTN_UT(t, zero, pd, clock = clock),
                     EMC2:::dRDMSWTN(t, ordinary, posdrift = pd))
    expect_equal(EMC2:::pRDMSWTN_UT(t, zero, pd, clock = clock),
                 EMC2:::pRDMSWTN(t, ordinary, posdrift = pd), tolerance = 1e-13)
  }
  set.seed(401)
  n <- 40
  dat <- data.frame(subjects = factor(1),
                    R = factor(sample(c("a", "b"), n, TRUE), levels = c("a", "b")),
                    rt = .2 + rgamma(n, 3, 6))
  dat$R[1:3] <- NA; dat$rt[1:3] <- Inf
  trunc <- dat; trunc$LT <- .25; trunc$UT <- 1.5
  trunc <- trunc[is.infinite(trunc$rt) | (trunc$rt > .25 & trunc$rt < 1.5), ]
  base_f <- list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1)
  for (d in list(dat, trunc)) for (pd in c(TRUE, FALSE)) for (sv in c(0, .3)) {
    ref_ctx <- ut_context(d, RDMSWTN(posdrift = pd), base_f)
    ref <- ut_ll(ref_ctx, ut_values(sampled_pars(ref_ctx$design, doMap = FALSE),
                                   v = if (pd) c(1.3, .9) else c(1.3, -.2),
                                   sv = sv, natural_v = !pd))
    for (clock in clocks) {
      ctx <- ut_context(d, RDMSWTN_UT(posdrift = pd, clock = clock))
      p <- ut_values(sampled_pars(ctx$design, doMap = FALSE), u = 0, sv = sv,
                     v = if (pd) c(1.3, .9) else c(1.3, -.2), natural_v = !pd)
      expect_equal(ut_ll(ctx, p), ref, tolerance = 1e-12)
    }
  }
})

test_that("clock inversion is stable across scales and at overflow", {
  x <- 10^seq(-6, 2, by = .25)
  for (clock in clocks) for (u in 10^seq(-12, 3)) {
    y <- ut_q(x, u, clock)
    keep <- is.finite(y)
    back <- EMC2:::rdmswtn_clock_qinv(y[keep], u, ut_code[[clock]])
    expect_equal(back, x[keep], tolerance = 1e-13)
  }
  expect_equal(EMC2:::rdmswtn_clock_qinv(1e308, 10, 1L), sqrt(2 / 10) * sqrt(1e308),
               tolerance = 1e-12)
  expect_equal(EMC2:::rdmswtn_clock_qinv(1e308, 10, 2L),
               (log(10) + log(1e308)) / 10, tolerance = 1e-14)
  expect_identical(EMC2:::rdmswtn_clock_qinv(c(0, Inf), .5, 1L), c(0, Inf))
  expect_true(is.na(EMC2:::rdmswtn_clock_qinv(1, -1, 1L)))
})

test_that("large operational times stay finite and consistent", {
  # Slow trials under strong urgency push q far into the ordinary tail.
  pars <- ut_pars(u = 20)
  t <- c(5, 10, 30)
  x <- t - pars[1, "t0"]
  q <- ut_q(x, 20, "linear")
  logf <- EMC2:::dRDMSWTN_UT(t, pars, log = TRUE)
  expect_true(all(is.finite(logf)))
  ordinary_logf <- vapply(q, function(qi) EMC2:::drdmswtn(
    qi, pars[1, "v"], pars[1, "b"], pars[1, "A"], pars[1, "s"], 0,
    pars[1, "sv"], log_out = TRUE), numeric(1))
  expect_equal(logf, ordinary_logf + log1p(20 * x), tolerance = 1e-9)
  expect_true(all(diff(logf) < 0))
  # Exponential operational time overflows: zero density, full CDF.
  expect_identical(EMC2:::dRDMSWTN_UT(200, pars, clock = "exponential"), 0)
  expect_equal(EMC2:::pRDMSWTN_UT(200, pars, clock = "exponential"), 1)
  # Invalid clocks are zero-density rows, never NaN.
  bad <- pars; bad[, "u"] <- -1
  expect_identical(EMC2:::dRDMSWTN_UT(.5, bad), 0)
  expect_identical(EMC2:::pRDMSWTN_UT(.5, bad), 0)
})

test_that("urgency shortens every quantile, the slow ones most", {
  quantile_at <- function(p, pars, clock) uniroot(
    function(t) EMC2:::pRDMSWTN_UT(t, pars, clock = clock) - p,
    c(pars[1, "t0"] + 1e-9, 50), tol = 1e-12)$root
  probs <- c(.1, .5, .9)
  for (clock in clocks) {
    qs <- sapply(c(0, .5, 1, 2), function(u)
      sapply(probs, quantile_at, pars = ut_pars(u = u), clock = clock))
    expect_true(all(diff(t(qs)) < 0))          # every quantile falls with u
    shrink <- qs[, 1] - qs[, 4]
    expect_true(all(diff(shrink) > 0))         # most at the slow quantiles
  }
})

test_that("a clock shared by every accumulator leaves choice probabilities unchanged", {
  a <- ut_pars(v = 1.3); b <- ut_pars(v = .9)
  choice <- function(u, clock) {
    a[, "u"] <- u; b[, "u"] <- u
    integrate(function(t) EMC2:::dRDMSWTN_UT(t, a, clock = clock) *
                (1 - EMC2:::pRDMSWTN_UT(t, b, clock = clock)),
              .1, Inf, rel.tol = 1e-10)$value
  }
  for (clock in clocks)
    expect_equal(choice(2, clock), choice(0, clock), tolerance = 1e-7)
})

test_that("compiled and R simulators follow the analytic CDF", {
  n <- 4000L
  lR <- factor(rep("a", n), levels = "a")
  for (clock in clocks) for (pd in c(TRUE, FALSE)) {
    row <- ut_pars(posdrift = pd, v = if (pd) 1.2 else .4)
    pars <- row[rep(1L, n), , drop = FALSE]
    F_inf <- EMC2:::pRDMSWTN_UT(Inf, row, pd, clock = clock)
    cdf <- function(t) EMC2:::pRDMSWTN_UT(t, row, pd, clock = clock) / F_inf
    set.seed(203)
    cpp <- EMC2:::rrdmswtn_clock_cpp(pars, "a", rep(TRUE, n), pd, ut_code[[clock]])$rt
    set.seed(204)
    ref <- EMC2:::rRDMSWTN_UT(lR, pars, posdrift = pd, clock = clock)$rt
    for (rt in list(cpp, ref)) {
      expect_equal(mean(is.finite(rt)), F_inf, tolerance = 4 * sqrt(.25 / n))
      expect_gt(suppressWarnings(
        ks.test(rt[is.finite(rt)], function(t) vapply(t, cdf, 0))$p.value), 1e-3)
    }
  }
  # make_data takes the compiled route.
  des <- design(factors = list(S = "target", subjects = 1), Rlevels = c("a", "b"),
                formula = ut_formula(), constants = c(s = 0),
                model = RDMSWTN_UT(clock = "exponential"), report_p_vector = FALSE)
  p <- ut_values(sampled_pars(des, doMap = FALSE))
  set.seed(205)
  sim <- make_data(p, design = des, n_trials = 40, rt_resolution = NULL)
  expect_true(all(is.finite(sim$rt)))
})

test_that("the compiled and R copula simulators draw the same finishing times", {
  n <- 200L
  lR <- factor(rep(c("a", "b"), n), levels = c("a", "b"))
  for (clock in clocks) {
    base <- rbind(c(ut_pars(v = 1.3)[1, ], rho = .6),
                  c(ut_pars(v = .9)[1, ], rho = .6))
    pars <- base[rep(1:2, n), ]
    set.seed(206)
    cpp <- EMC2:::rrdmswtn_clock_corr_cpp(pars, levels(lR), rep(TRUE, 2L * n),
                                          TRUE, ut_code[[clock]])
    set.seed(206)
    ref <- EMC2:::rRDMSWTN_UT_corr(lR, pars, clock = clock)
    expect_equal(cpp$rt, ref$rt, tolerance = 1e-7)
    expect_identical(levels(lR)[cpp$R], as.character(ref$R))
    # rho = 0 falls back to the independent simulator exactly.
    pars[, "rho"] <- 0
    set.seed(207)
    a <- EMC2:::rrdmswtn_clock_corr_cpp(pars, levels(lR), rep(TRUE, 2L * n),
                                        TRUE, ut_code[[clock]])
    set.seed(207)
    b <- EMC2:::rrdmswtn_clock_cpp(pars, levels(lR), rep(TRUE, 2L * n),
                                   TRUE, ut_code[[clock]])
    expect_identical(a$rt, b$rt)
  }
  bad <- pars; bad[, "u"] <- -1
  expect_error(EMC2:::rrdmswtn_clock_cpp(bad, levels(lR), rep(TRUE, 2L * n),
                                         TRUE, 1L),
               "u must be finite and nonnegative")
})

test_that("the compiled likelihood matches a direct race computation", {
  dat <- data.frame(
    subjects = factor(1), R = factor(c("a", "b", NA), levels = c("a", "b")),
    rt = c(.7, 1.3, Inf), LT = c(.2, 0, 0), UT = c(1.1, Inf, Inf)
  )
  for (clock in clocks) for (pd in c(TRUE, FALSE)) {
    ctx <- ut_context(dat, RDMSWTN_UT(posdrift = pd, clock = clock))
    p <- ut_values(sampled_pars(ctx$design, doMap = FALSE),
                   v = if (pd) c(1.3, .9) else c(.2, -.1), natural_v = !pd)
    pars <- ut_mapped(ctx, p)
    d <- function(t, r) EMC2:::dRDMSWTN_UT(t, pars[r, , drop = FALSE], pd, clock = clock)
    S <- function(t, rows) prod(1 - EMC2:::pRDMSWTN_UT(
      t, pars[rows, , drop = FALSE], pd, clock = clock))
    direct <- c(
      log(d(.7, 1)) + log(S(.7, 2)) -
        log(S(.2, 1:2) - S(1.1, 1:2) + S(Inf, 1:2)),
      log(d(1.3, 4)) + log(S(1.3, 3)),
      log(max(S(Inf, 5:6), 1e-10))
    )
    expect_equal(ut_ll(ctx, p, pointwise = TRUE), direct, tolerance = 1e-8)
  }
})

test_that("correlated likelihoods nest rho = 0 and match the copula density", {
  dat <- data.frame(subjects = factor(1),
                    R = factor(c("a", "b", "a"), levels = c("a", "b")),
                    rt = c(.3, .6, 1.2))
  for (clock in clocks) {
    ind_ctx <- ut_context(dat, RDMSWTN_UT(clock = clock))
    ind <- ut_ll(ind_ctx, ut_values(sampled_pars(ind_ctx$design, doMap = FALSE)))
    ctx <- ut_context(dat, RDMSWTN_UTcorr(clock = clock))
    p0 <- ut_values(sampled_pars(ctx$design, doMap = FALSE), rho = 0)
    expect_equal(ut_ll(ctx, p0), ind, tolerance = 1e-10)

    p <- ut_values(p0, rho = .6)
    pars <- ut_mapped(ctx, p)
    rho <- pars[1, "rho"]
    direct <- vapply(1:3, function(i) {
      rows <- 2 * i - 1:0
      w <- rows[dat$R[i] == c("a", "b")]; l <- setdiff(rows, w)
      Fw <- EMC2:::pRDMSWTN_UT(dat$rt[i], pars[w, , drop = FALSE], clock = clock)
      Fl <- EMC2:::pRDMSWTN_UT(dat$rt[i], pars[l, , drop = FALSE], clock = clock)
      log(EMC2:::dRDMSWTN_UT(dat$rt[i], pars[w, , drop = FALSE], clock = clock)) +
        pnorm((qnorm(Fl) - rho * qnorm(Fw)) / sqrt(1 - rho^2),
              lower.tail = FALSE, log.p = TRUE)
    }, numeric(1))
    expect_equal(ut_ll(ctx, p, pointwise = TRUE), direct, tolerance = 1e-10)

    # Correlated drift draws: rho = 0 is the independent race.
    dctx <- ut_context(dat, RDMSWTN_UTcorr(clock = clock, correlate = "drifts"),
                       matchfun = function(d) d$lR == "a")
    expect_equal(ut_ll(dctx, ut_values(sampled_pars(dctx$design, doMap = FALSE),
                                       rho = 0)),
                 ind, tolerance = 1e-8)
  }
})

test_that("malformed correlated designs are rejected", {
  dat3 <- data.frame(subjects = factor(1),
                     R = factor("a", levels = c("a", "b", "c")), rt = .7)
  expect_error(ut_context(dat3, RDMSWTN_UTcorr()), "at most two accumulator rows")
  model <- RDMSWTN_UTcorr(posdrift = FALSE)
  pars <- matrix(c(1, 1, .2, .1, 1, .3, .5, 0, 0, .5), nrow = 1,
                 dimnames = list(NULL, names(model$p_types)))
  expect_error(model$Ttransform(pars, NULL),
               "requires sv = 0 on the correlated rows")
  pars[, "sv"] <- 0
  expect_silent(model$Ttransform(pars, NULL))
})
