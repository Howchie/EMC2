skip_model_validation()

# BAwR, the ballistic accumulator with a ramping clearance rate --
# src/model_BAwR.h, R/model_BAwR.R.
#
# dX/du = V - kappa u^p, so X(u) = z + V u - kappa u^(p+1)/(p+1),
# b = B + A, z ~ U(0, A).
#
# Oracles first.  The kappa = 0 member is exactly the LBA, which pins the
# shared geometry; past that there is no closed-form oracle, so the kernels are
# checked against references written from the model definition (a numerical
# minimisation of the required-launch curve and quadrature over the start
# point), against Monte Carlo, and against internal consistency
# (integrate(f) == F).  The references below share no code with the compiled
# kernels.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

# Required launch strength to be exactly at b at time u from start point z.
ref_vstar_p <- function(u, z, b, kap, pw) (b - z) / u + kap * u^pw / (pw + 1)

# Survivor and density of the launch strength.
# For launch = 1L the public pair is (mean, cv); log V ~ N(mu, sigma^2) with
# sigma^2 = log1p(cv^2) and mu = log(mean) - sigma^2/2.
ln_musigma <- function(m, cv) {
  s2 <- log1p(unname(cv)^2)
  c(mu = log(unname(m)) - s2 / 2, sigma = sqrt(s2))
}
ref_gbar_p <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) {
    q <- ln_musigma(p1, p2)
    return(plnorm(w, q[["mu"]], q[["sigma"]], lower.tail = FALSE))
  }
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  pnorm(w, p1, p2, lower.tail = FALSE) / den
}
ref_g_p <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) {
    q <- ln_musigma(p1, p2)
    return(dlnorm(w, q[["mu"]], q[["sigma"]]))
  }
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  dnorm(w, p1, p2) / den
}

# F(u) = (1/A) int_0^A P(V >= min_{s <= u} V*(s, z)) dz.  The inner minimum is
# found numerically from the crossing law, with no appeal to the tangency
# formula the kernel uses.
ref_F_bawr <- function(u, p1, p2, b, A, kap, pw, launch = 1L,
                       posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  vmin <- function(z) {
    if (kap <= 0) return(ref_vstar_p(u, z, b, kap, pw))
    hi <- min(u, ((pw + 1) * (b - z) / (kap * pw))^(1 / (pw + 1)))
    if (!isTRUE(hi > 1e-12)) return(ref_vstar_p(u, z, b, kap, pw))
    o <- stats::optimize(function(s) ref_vstar_p(s, z, b, kap, pw),
                         c(1e-12, hi), tol = 1e-12)$objective
    min(o, ref_vstar_p(hi, z, b, kap, pw))
  }
  if (A <= 0) return(ref_gbar_p(vmin(0), p1, p2, launch, posdrift))
  stats::integrate(function(zz) {
    vapply(zz, function(z) ref_gbar_p(vmin(z), p1, p2, launch, posdrift),
           numeric(1))
  }, 0, A, rel.tol = 1e-11)$value / A
}

# f(u) = (1/A) int_0^A g(V*(u, z)) [-dV*/du]_+ dz, with
# -dV*/du = (V* - kappa u^p) / u from elementary calculus on V*.
ref_f_bawr <- function(u, p1, p2, b, A, kap, pw, launch = 1L,
                       posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  cc <- kap * u^pw
  dens <- function(z) {
    w <- ref_vstar_p(u, z, b, kap, pw)
    jac <- (w - cc) / u
    if (!isTRUE(jac > 0)) return(0)
    ref_g_p(w, p1, p2, launch, posdrift) * jac
  }
  if (A <= 0) return(dens(0))
  stats::integrate(function(zz) vapply(zz, dens, numeric(1)), 0, A,
                   rel.tol = 1e-11)$value / A
}

# First passage by root finding on the trajectory itself.
ref_fp_bawr <- function(V, z, b, kap, pw) {
  if (!isTRUE(b > z)) return(0)
  if (!isTRUE(V > 0)) return(Inf)
  if (kap <= 1e-12) return((b - z) / V)
  up <- (V / kap)^(1 / pw)                       # trajectory peak
  Gx <- function(u) V * u - kap * u^(pw + 1) / (pw + 1) - (b - z)
  if (Gx(up) < 0) return(Inf)
  stats::uniroot(Gx, c(1e-14, up), tol = 1e-13)$root
}

# The plain lognormal launch (launch = 1L) is sampled as the natural-scale
# mean and CV; the reference integrals in this file work in (mu, sigma), so
# these wrappers convert at the kernel boundary.  The split launch (2L) is not
# reparameterized and passes straight through.
d_bawr <- function(t, A, b, p1, p2, kap, pw, launch = 1L, posdrift = TRUE,
                   log_out = FALSE) {
  EMC2:::dbawr(t, A, b, p1, p2, kap, pw, launch = as.integer(launch),
               posdrift = posdrift, log_out = log_out)
}
p_bawr <- function(t, A, b, p1, p2, kap, pw, launch = 1L, posdrift = TRUE,
                   log_out = FALSE) {
  EMC2:::pbawr(t, A, b, p1, p2, kap, pw, launch = as.integer(launch),
               posdrift = posdrift, log_out = log_out)
}
t_max_p <- function(A, b, kap, pw) EMC2:::bawr_tmax(A, b, kap, pw)

# ---------------------------------------------------------------------------
# 1. Geometry: the endpoint, and its dependence on b
# ---------------------------------------------------------------------------

test_that("T_max is the argmin of the required-launch curve", {
  for (pw in c(0.5, 1, 2.3)) {
    for (kap in c(0.3, 1.1, 3)) {
      for (b in c(0.6, 1.3)) {
        got <- t_max_p(0, b, kap, pw)
        num <- stats::optimize(function(u) ref_vstar_p(u, 0, b, kap, pw),
                               c(1e-8, 60), tol = 1e-14)$minimum
        expect_equal(got, num, tolerance = 1e-5)
        # Closed form from the model definition, independently written.
        expect_equal(got, ((pw + 1) * b / (kap * pw))^(1 / (pw + 1)),
                     tolerance = 1e-12)
      }
    }
  }
  expect_true(is.infinite(t_max_p(0.2, 1, 0, 1)))   # LBA limit
  expect_true(is.na(t_max_p(0.2, 1, 1, 0)))        # p = 0 is not in the family
  expect_true(is.na(t_max_p(0.2, 1, 1, -1)))
})

test_that("T_max depends on b but not on A or the launch distribution", {
  # The substantive contrast with BAwF, whose endpoint is free of b: here the
  # decay is in physical time and knows nothing about the threshold, so a more
  # cautious accumulator simply gets longer before the drive expires.
  # T_max scales as b^(1/(p+1)), so at p = 1 doubling b multiplies it by
  # sqrt(2).
  tm <- EMC2:::bawr_tmax_vec(0.3, c(0.6, 1.2, 2.4, 4.8), 1.1, 1)
  expect_length(tm, 4L)
  expect_equal(tm[-1] / tm[-4], rep(sqrt(2), 3), tolerance = 1e-12)
  tm3 <- EMC2:::bawr_tmax_vec(0.3, c(1, 8), 1.1, 2)
  expect_equal(tm3[2] / tm3[1], 2, tolerance = 1e-12)   # b^(1/3)
  # A is absent from the endpoint (it sets where the EARLIEST saturation is,
  # not the latest), and so is every launch parameter.
  expect_equal(t_max_p(0, 1.4, 1.1, 1), t_max_p(0.9, 1.4, 1.1, 1),
               tolerance = 1e-12)
})

test_that("V_c(0) = kappa T_max^p is the omission boundary", {
  A <- c(0, 0.3); b <- c(1, 2); kap <- c(1.5, 0.7); pw <- c(1, 2.3)
  vc <- EMC2:::bawr_vcrit_vec(A, b, kap, pw)
  tm <- EMC2:::bawr_tmax_vec(A, b, kap, pw)
  expect_equal(vc, kap * tm^pw, tolerance = 1e-12)
  # A launch exactly at V_c(0) from z = 0 is tangent: it reaches b at T_max and
  # nowhere else, so F(Inf) for a point start is P(V >= V_c(0)).
  vc0 <- EMC2:::bawr_vcrit_vec(0, 1.2, 1.3, 1)
  q <- ln_musigma(1.0, 0.5)          # public (mean, cv) -> (mu, sigma)
  expect_equal(p_bawr(Inf, 0, 1.2, 1.0, 0.5, 1.3, 1, launch = 1L),
               plnorm(vc0, q[["mu"]], q[["sigma"]], lower.tail = FALSE),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. CDF against an independent start-point quadrature
# ---------------------------------------------------------------------------

bawr_cases <- list(
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, kap = 1.2, pw = 1,   launch = 1L),
  list(A = 0.0, b = 1.0, p1 = 1.2, p2 = 0.4, kap = 0.8, pw = 1,   launch = 1L),
  list(A = 0.6, b = 1.5, p1 = 4.0, p2 = 1.5, kap = 1.0, pw = 1,   launch = 0L),
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, kap = 1.2, pw = 0.5, launch = 1L),
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, kap = 1.2, pw = 2.3, launch = 1L),
  list(A = 0.5, b = 1.5, p1 = 3.0, p2 = 1.2, kap = 2.0, pw = 2.3, launch = 0L),
  list(A = 0.3, b = 1.0, p1 = 0.8, p2 = 0.9, kap = 3.0, pw = 1,   launch = 1L))

test_that("the CDF matches quadrature over the start point", {
  for (cs in bawr_cases) {
    Tm <- t_max_p(cs$A, cs$b, cs$kap, cs$pw)
    us <- Tm * c(0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
    got <- p_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                  launch = cs$launch)
    ref <- vapply(us, ref_F_bawr, numeric(1), cs$p1, cs$p2, cs$b, cs$A,
                  cs$kap, cs$pw, cs$launch)
    expect_equal(got, ref, tolerance = 1e-6)
    expect_true(all(diff(got) >= -1e-12))          # monotone
  }
})

test_that("the PDF matches quadrature over the start point", {
  for (cs in bawr_cases) {
    Tm <- t_max_p(cs$A, cs$b, cs$kap, cs$pw)
    us <- Tm * c(0.1, 0.3, 0.5, 0.7, 0.9)
    got <- d_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                  launch = cs$launch)
    ref <- vapply(us, ref_f_bawr, numeric(1), cs$p1, cs$p2, cs$b, cs$A,
                  cs$kap, cs$pw, cs$launch)
    expect_equal(got, ref, tolerance = 1e-6)
    expect_true(all(got >= 0))
  }
})

test_that("the log arguments agree with the natural ones", {
  for (cs in bawr_cases[1:4]) {
    Tm <- t_max_p(cs$A, cs$b, cs$kap, cs$pw)
    us <- Tm * c(0.2, 0.6, 0.9)
    expect_equal(d_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                        launch = cs$launch, log_out = TRUE),
                 log(d_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                            launch = cs$launch)), tolerance = 1e-10)
    expect_equal(p_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                        launch = cs$launch, log_out = TRUE),
                 log(p_bawr(us, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                            launch = cs$launch)), tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# 3. Internal consistency
# ---------------------------------------------------------------------------

test_that("the PDF integrates to the CDF", {
  for (cs in bawr_cases) {
    Tm <- t_max_p(cs$A, cs$b, cs$kap, cs$pw)
    u1 <- Tm * 0.85
    I <- stats::integrate(function(x)
      d_bawr(x, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw, launch = cs$launch),
      1e-9, u1, subdivisions = 2000, rel.tol = 1e-10)$value
    expect_equal(I, p_bawr(u1, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                           launch = cs$launch), tolerance = 1e-8)
  }
})

test_that("F(Inf) is the finishing probability and the tail is defective", {
  for (cs in bawr_cases) {
    Finf <- p_bawr(Inf, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                   launch = cs$launch)
    ref <- ref_F_bawr(t_max_p(cs$A, cs$b, cs$kap, cs$pw), cs$p1, cs$p2, cs$b,
                      cs$A, cs$kap, cs$pw, cs$launch)
    expect_equal(Finf, ref, tolerance = 1e-6)
    expect_lt(Finf, 1)                 # always defective: weak launches omit
    # Past T_max the CDF is flat: no response can occur after the endpoint.
    Tm <- t_max_p(cs$A, cs$b, cs$kap, cs$pw)
    expect_equal(p_bawr(Tm * 3, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                        launch = cs$launch), Finf, tolerance = 1e-12)
    expect_equal(d_bawr(Tm * 3, cs$A, cs$b, cs$p1, cs$p2, cs$kap, cs$pw,
                        launch = cs$launch), 0)
  }
})

test_that("the endpoint shutdown is quadratic for A > 0 and linear for A = 0", {
  # Tangency kills the temporal Jacobian linearly; with A > 0 the live
  # start-point interval collapses too, supplying the second factor.  A point
  # start has no such interval, so it keeps the linear order.
  for (pw in c(1, 2.3)) {
    kap <- 1.3; b <- 1.2
    eps <- c(1e-2, 5e-3, 2.5e-3)
    for (A in c(0, 0.4)) {
      Tm <- t_max_p(A, b, kap, pw)
      fv <- d_bawr(Tm - eps, A, b, 1.0, 0.5, kap, pw, launch = 1L)
      expect_true(all(fv > 0))
      ord <- log(fv[-1] / fv[-length(fv)]) / log(eps[-1] / eps[-length(eps)])
      expect_equal(as.numeric(ord), rep(if (A > 0) 2 else 1, 2),
                   tolerance = 0.05)
    }
  }
})

# ---------------------------------------------------------------------------
# 4. The kappa = 0 oracle
# ---------------------------------------------------------------------------

test_that("kappa = 0 is exactly the LBA, for every p", {
  ts <- c(0.2, 0.5, 1, 2)
  lba_d <- EMC2:::dlba(ts, 0.3, 1.0, 2.0, 0.8, posdrift = TRUE)
  lba_p <- EMC2:::plba(ts, 0.3, 1.0, 2.0, 0.8, posdrift = TRUE)
  for (pw in c(0.5, 1, 2.3)) {
    expect_equal(d_bawr(ts, 0.3, 1.0, 2.0, 0.8, 0, pw, launch = 0L), lba_d,
                 tolerance = 1e-12)
    expect_equal(p_bawr(ts, 0.3, 1.0, 2.0, 0.8, 0, pw, launch = 0L), lba_p,
                 tolerance = 1e-12)
    # Continuity into the limit, not just the exact branch.
    expect_equal(d_bawr(ts, 0.3, 1.0, 2.0, 0.8, 1e-9, pw, launch = 0L), lba_d,
                 tolerance = 1e-7)
  }
})

# ---------------------------------------------------------------------------
# 5. Monte Carlo
# ---------------------------------------------------------------------------

test_that("the CDF matches a brute-force simulation of the trajectory", {
  skip_on_cran()
  set.seed(11)
  A <- 0.4; b <- 1.2; mu <- 1.0; sg <- 0.5; kap <- 1.2; pw <- 1  # (mean, cv)
  N <- 1000
  qms <- ln_musigma(mu, sg)
  V <- rlnorm(N, qms[["mu"]], qms[["sigma"]])
  z <- runif(N, 0, A)
  hit <- vapply(seq_len(N), function(i) ref_fp_bawr(V[i], z[i], b, kap, pw),
                numeric(1))
  Tm <- t_max_p(A, b, kap, pw)
  us <- Tm * c(0.2, 0.5, 0.8, 0.95)
  emp <- vapply(us, function(u) mean(hit <= u), numeric(1))
  got <- p_bawr(us, A, b, mu, sg, kap, pw, launch = 1L)
  # Absolute tolerance: expect_equal's is relative, which is far too loose in
  # the left tail and far too tight nowhere.
  expect_true(all(abs(got - emp) < 4 * sqrt(0.25 / N) + 1e-3))
  expect_lt(abs(mean(is.finite(hit)) -
                p_bawr(Inf, A, b, mu, sg, kap, pw, launch = 1L)),
            4 * sqrt(0.25 / N) + 1e-3)
  expect_true(all(hit[is.finite(hit)] <= Tm + 1e-9))   # the hard endpoint
})

# ---------------------------------------------------------------------------
# 6. The model object
# ---------------------------------------------------------------------------

test_that("BAwR() exposes the documented parameters and c_name", {
  m <- BAwR()
  expect_equal(m$c_name, "BAwR_LOGN")
  expect_equal(m$p_types_canonical,
               c("mean", "cv", "B", "A", "t0", "kappa", "p"))
  expect_equal(BAwR("normal")$c_name, "BAwR")
  expect_equal(BAwR("normal", posdrift = FALSE)$c_name, "BAwRIO")
  expect_equal(BAwR("normal")$p_types_canonical,
               c("v", "sv", "B", "A", "t0", "kappa", "p"))
  # kappa = 0 must stay reachable (it is the LBA limit); p = 0 must not.
  expect_true("kappa" %in% names(m$bound$exception))
  expect_false("p" %in% names(m$bound$exception))
  expect_error(BAwR("lognormal", posdrift = FALSE), "posdrift only applies")
})

test_that("Ttransform reports b, Tmax, rt_max and Vcrit", {
  m <- BAwR()
  pars <- cbind(mean = c(1, 1), cv = c(0.5, 0.5), B = c(0.9, 3.6),
                A = c(0.3, 0.3), t0 = c(0.15, 0.15), kappa = c(1.1, 1.1),
                p = c(1, 1))
  out <- m$Ttransform(pars, NULL)
  expect_equal(out[, "b"], pars[, "B"] + pars[, "A"])
  expect_equal(out[, "Tmax"],
               EMC2:::bawr_tmax_vec(pars[, "A"], out[, "b"], pars[, "kappa"],
                                    pars[, "p"]))
  expect_equal(out[, "rt_max"], pars[, "t0"] + out[, "Tmax"])
  expect_equal(out[, "Vcrit"],
               EMC2:::bawr_vcrit_vec(pars[, "A"], out[, "b"], pars[, "kappa"],
                                     pars[, "p"]))
  # Quadrupling B multiplies b by ~4 and hence Tmax by 2 at p = 1: caution
  # buys time here, unlike in BAwF.
  expect_gt(out[2, "Tmax"], out[1, "Tmax"])
})

# ---------------------------------------------------------------------------
# 7. The compiled race likelihood
# ---------------------------------------------------------------------------

ref_race_ll_bawr <- function(dadm, pars, launch, min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mean", "cv") else c("v", "sv")
  lp <- function(i) unname(c(pars[i, nm[1]], pars[i, nm[2]]))
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    if (is.infinite(u)) u <- t_max_p(pars[i, "A"], pars[i, "b"],
                                     pars[i, "kappa"], pars[i, "p"])
    q <- lp(i)
    ref_F_bawr(u, q[1], q[2], pars[i, "b"], pars[i, "A"],
               pars[i, "kappa"], pars[i, "p"], launch)
  }
  fu <- function(i, u) {
    q <- lp(i)
    ref_f_bawr(u, q[1], q[2], pars[i, "b"], pars[i, "A"],
               pars[i, "kappa"], pars[i, "p"], launch)
  }
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    if (is.infinite(rt) && rt > 0) {
      term <- prod(vapply(idx, function(i) 1 - Fu(i, Inf), numeric(1)))
    } else {
      w <- idx[which(dadm$winner[idx])]
      term <- fu(w, rt - pars[w, "t0"])
      for (i in setdiff(idx, w))
        term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    }
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

bawr_mk <- function(model, form, consts, dat) {
  des <- design(data = dat, model = model, matchfun = function(d) d$S == d$lR,
                formula = form, constants = consts)
  list(des = des, emc = suppressMessages(
    make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
             rt_resolution = NULL)))
}

bawr_ll <- function(fx, p, p_types_override = NULL) {
  model <- fx$emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- fx$emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  if (!is.null(p_types_override)) p_types <- p_types_override
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform, p_types = p_types,
                    min_ll = log(1e-10), trend = model$trend)
}

bawr_dat <- function(n = 60) {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat[seq_len(n), ]
}

bawr_form <- list(mean ~ 1, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, kappa ~ 1, p ~ 1)
bawr_p <- c(mean = log(2.5), cv = log(0.6), B = log(0.8), A = log(0.3),
            t0 = log(0.15), kappa = log(0.8), p = log(1))

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  dat <- bawr_dat()
  p <- c(bawr_p, v = 3, sv = log(1))
  fx <- list(
    ln = suppressMessages(bawr_mk(BAwR, bawr_form, NULL, dat)),
    no = suppressMessages(bawr_mk(function() BAwR("normal"),
      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, kappa ~ 1, p ~ 1), c(sv = log(1)),
      dat)))
  for (which in c("ln", "no")) {
    launch <- if (which == "ln") 1L else 0L
    pp <- p[names(sampled_pars(fx[[which]]$des))]
    got <- bawr_ll(fx[[which]], pp)
    expect_true(is.finite(got))
    dadm <- fx[[which]]$emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(pp, dadm, fx[[which]]$emc[[1]]$model())
    expect_equal(got, ref_race_ll_bawr(dadm, pars, launch), tolerance = 1e-6)
  }
})

test_that("a non-unit p survives the whole compiled path", {
  # p = 1 is the member with the tidiest algebra; a fit with p != 1 exercises
  # the power-stop-loss frozen branch through the sampled likelihood.
  skip_on_cran()
  dat <- bawr_dat()
  p <- bawr_p; p["p"] <- log(2.3)
  fx <- suppressMessages(bawr_mk(BAwR, bawr_form, NULL, dat))
  got <- bawr_ll(fx, p[names(sampled_pars(fx$des))])
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p[names(sampled_pars(fx$des))], dadm,
                                    fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll_bawr(dadm, pars, 1L), tolerance = 1e-6)
})

test_that("the compiled BAwR likelihood reproduces the LBA at kappa = 0", {
  # A second, independent route to the kappa = 0 oracle: two models, two
  # adapter branches, two column layouts, one answer.  This is what catches a
  # wiring error the kernel tests cannot see (wrong column order, wrong
  # t0_index, a stray posdrift default).
  skip_on_cran()
  dat <- bawr_dat()
  fl <- suppressMessages(bawr_mk(LBA, list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                                 c(sv = log(1)), dat))
  fp <- suppressMessages(bawr_mk(function() BAwR("normal"),
                                 list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                                 c(sv = log(1), kappa = log(0), p = log(1)),
                                 dat))
  pv <- c(v = 3, B = log(0.9), A = log(0.3), t0 = log(0.15))
  l_lba <- bawr_ll(fl, pv[names(sampled_pars(fl$des))])
  l_bawr <- bawr_ll(fp, pv[names(sampled_pars(fp$des))])
  expect_true(is.finite(l_lba))
  expect_equal(l_bawr, l_lba, tolerance = 1e-8)
})

test_that("a p_types reordering is caught by the column contract", {
  skip_on_cran()
  dat <- bawr_dat()
  fx <- suppressMessages(bawr_mk(BAwR, bawr_form, NULL, dat))
  swapped <- names(BAwR()$p_types)
  swapped[1:2] <- swapped[2:1]
  expect_error(bawr_ll(fx, bawr_p[names(sampled_pars(fx$des))],
                       p_types_override = swapped),
               "BAwR_LOGN kernels expect parameter column")
  # kappa and p must not be interchangeable either: they are adjacent, share a
  # log transform, and swapping them would otherwise fit silently.
  swapped2 <- names(BAwR()$p_types)
  i <- match(c("kappa", "p"), swapped2)
  swapped2[i] <- swapped2[rev(i)]
  expect_error(bawr_ll(fx, bawr_p[names(sampled_pars(fx$des))],
                       p_types_override = swapped2),
               "BAwR_LOGN kernels expect parameter column")
})

test_that("omissions past t0 + T_max stay well posed", {
  skip_on_cran()
  dat <- bawr_dat(80)
  dat$rt[1:8] <- Inf                 # intrinsic omissions
  dat$R[1:8] <- NA                   # with an unknown response
  fx <- suppressMessages(bawr_mk(BAwR, bawr_form, NULL, dat))
  got <- bawr_ll(fx, bawr_p[names(sampled_pars(fx$des))])
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(bawr_p[names(sampled_pars(fx$des))], dadm,
                                    fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll_bawr(dadm, pars, 1L), tolerance = 1e-6)
  # A finite response beyond t0 + T_max has exactly zero density and must
  # floor rather than return -Inf or NaN.
  dat2 <- dat
  dat2$rt[1:8] <- 0.15 + t_max_p(0.3, 1.1, 0.8, 1) + 5
  fx2 <- suppressMessages(bawr_mk(BAwR, bawr_form, NULL, dat2))
  expect_true(is.finite(bawr_ll(fx2, bawr_p[names(sampled_pars(fx2$des))])))
})

# ---------------------------------------------------------------------------
# 8. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree with each other and with the CDF", {
  skip_on_cran()
  pars <- cbind(mean = 1.0, cv = 0.5, b = 1.2, A = 0.4, t0 = 0.15,
                kappa = 1.2, p = 1)
  n <- 1000
  pm <- pars[rep(1, 2 * n), , drop = FALSE]
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  set.seed(7)
  cpp <- EMC2:::.rfun_BAwR(lR, pm, rep(TRUE, 2 * n), launch = 1L)
  set.seed(7)
  withr::with_options(list(emc2.cpp_rfun = FALSE), {
    rr <- EMC2:::.rfun_BAwR(lR, pm, rep(TRUE, 2 * n), launch = 1L)
  })
  # The two paths consume the RNG stream differently (rlnorm vs
  # exp(mu + sigma * norm_rand)), so only the distributions may be compared.
  fin <- is.finite(cpp$rt)
  expect_equal(mean(fin), mean(is.finite(rr$rt)), tolerance = 0.03)
  # The winner's RT is the min of two iid first passages, so its CDF is
  # 1 - (1 - F)^2 with F the single-accumulator CDF.
  qs <- stats::quantile(cpp$rt[fin], c(0.25, 0.5, 0.75))
  emp <- vapply(qs, function(q) mean(cpp$rt <= q, na.rm = TRUE), numeric(1))
  Fq <- p_bawr(as.numeric(qs) - 0.15, 0.4, 1.2, 1.0, 0.5, 1.2, 1, launch = 1L)
  expect_lt(max(abs(as.numeric(emp) - (1 - (1 - Fq)^2))), 0.06)
  expect_true(all(cpp$rt[fin] <= 0.15 + t_max_p(0.4, 1.2, 1.2, 1) + 1e-8))
})

test_that("the simulator's hit times match the reference root solve", {
  # The C++ and R solvers share a bracket argument that is easy to get wrong:
  # unlike BAwF, this one takes the DISTANCE b - z, not b and z separately.
  V <- c(0.8, 1.5, 2.0, 4.0, 0.2)
  z <- c(0.0, 0.1, 0.3, 0.2, 0.0)
  for (pw in c(0.5, 1, 2.3)) {
    got <- mapply(EMC2:::.bawr_hit_time, V, 1.2 - z, 1.1, pw)
    ref <- mapply(ref_fp_bawr, V, z, 1.2, 1.1, pw)
    expect_equal(got, ref, tolerance = 1e-9)
  }
})

test_that("make_data produces omissions the design can be fit back through", {
  skip_on_cran()
  dat <- bawr_dat(40)
  des <- suppressMessages(design(
    data = dat, model = BAwR, matchfun = function(d) d$S == d$lR,
    formula = bawr_form))
  set.seed(3)
  sim <- make_data(bawr_p, design = des, n_trials = 60)
  expect_true(any(is.infinite(sim$rt)))          # intrinsic omissions
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] <= 0.15 + t_max_p(0.3, 1.1, 0.8, 1) + 1e-8))
  emc <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(bawr_ll(list(emc = emc, des = des), bawr_p)))
})
