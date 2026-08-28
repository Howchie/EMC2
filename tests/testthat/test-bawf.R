skip_model_validation()

# BAwF, the ballistic accumulator with global fading -- src/model_BAwF.h,
# R/model_BAwF.R.
#
# X(u) = h_rho(u) [z + V u], b = B + A, z ~ U(0, A).
#
# Oracles first.  The k = 0 member is exactly the LBA, which pins the shared
# geometry; past that there is no closed-form oracle, so the kernels are
# checked against references written from the model definition (a numerical
# minimisation of the required-launch curve and quadrature over the start
# point), against Monte Carlo, and against internal consistency
# (integrate(f) == F).  The references below share no code with the compiled
# kernels.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

ref_H <- function(x, rho) if (is.finite(rho)) exp(rho * log1p(x / rho)) else exp(x)
ref_Hp <- function(x, rho) {
  if (is.finite(rho)) exp((rho - 1) * log1p(x / rho)) else exp(x)
}
ref_xmax <- function(rho) if (is.finite(rho)) rho / (rho - 1) else 1

# Required launch strength to be exactly at b at time u from start point z.
ref_vstar <- function(u, z, b, k, rho) (b * ref_H(k * u, rho) - z) / u

# Survivor and density of the launch strength.
ref_gbar <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) return(plnorm(w, p1, p2, lower.tail = FALSE))
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  pnorm(w, p1, p2, lower.tail = FALSE) / den
}
ref_g <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) return(dlnorm(w, p1, p2))
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  dnorm(w, p1, p2) / den
}

# F(u) = (1/A) int_0^A P(V >= min_{s <= u} V*(s, z)) dz.  The inner minimum is
# found numerically from the crossing law, with no appeal to the tangency
# formula the kernel uses.
ref_F_bawf <- function(u, p1, p2, b, A, k, rho = Inf, launch = 1L,
                       posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  vmin <- function(z) {
    if (k <= 0) return(ref_vstar(u, z, b, k, rho))
    hi <- min(u, ref_xmax(rho) / k)
    o <- stats::optimize(function(s) ref_vstar(s, z, b, k, rho),
                         c(1e-12, hi), tol = 1e-12)$objective
    min(o, ref_vstar(hi, z, b, k, rho))
  }
  if (A <= 0) return(ref_gbar(vmin(0), p1, p2, launch, posdrift))
  stats::integrate(function(zz) {
    vapply(zz, function(z) ref_gbar(vmin(z), p1, p2, launch, posdrift),
           numeric(1))
  }, 0, A, rel.tol = 1e-11)$value / A
}

# f(u) = (1/A) int_0^A g(V*(u, z)) [-dV*/du]_+ dz, with
# -dV*/du = (V* - k b H'(k u)) / u from elementary calculus on V*.
ref_f_bawf <- function(u, p1, p2, b, A, k, rho = Inf, launch = 1L,
                       posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  cc <- k * b * ref_Hp(k * u, rho)
  dens <- function(z) {
    w <- ref_vstar(u, z, b, k, rho)
    jac <- (w - cc) / u
    if (!isTRUE(jac > 0)) return(0)
    ref_g(w, p1, p2, launch, posdrift) * jac
  }
  if (A <= 0) return(dens(0))
  stats::integrate(function(zz) vapply(zz, dens, numeric(1)), 0, A,
                   rel.tol = 1e-11)$value / A
}

# First passage by root finding on the trajectory itself.
ref_fp_bawf <- function(V, z, b, k, rho = Inf) {
  if (!isTRUE(b > z)) return(0)
  if (!isTRUE(V > 0)) return(Inf)
  if (k <= 1e-12) return((b - z) / V)
  xm <- ref_xmax(rho)
  Fx <- function(x) b * ref_H(x, rho) - z - V * x / k
  if (Fx(xm) > 0) return(Inf)
  stats::uniroot(Fx, c(0, xm), tol = 1e-13)$root / k
}

d_bawf <- function(t, A, b, p1, p2, k, launch = 1L, posdrift = TRUE,
                   rho = Inf, log_out = FALSE) {
  EMC2:::dbawf(t, A, b, p1, p2, k, launch = as.integer(launch),
               posdrift = posdrift, log_out = log_out,
               rho = if (is.finite(rho)) rho else 0)
}
p_bawf <- function(t, A, b, p1, p2, k, launch = 1L, posdrift = TRUE,
                   rho = Inf, log_out = FALSE) {
  EMC2:::pbawf(t, A, b, p1, p2, k, launch = as.integer(launch),
               posdrift = posdrift, log_out = log_out,
               rho = if (is.finite(rho)) rho else 0)
}
t_max <- function(A, b, k, rho = Inf) {
  EMC2:::bawf_tmax(A, b, k, if (is.finite(rho)) rho else 0)
}

# ---------------------------------------------------------------------------
# 1. Geometry: the endpoint belongs to k alone
# ---------------------------------------------------------------------------

test_that("T_max is x_max/k and no other parameter can move it", {
  # This is the whole reason the model exists: B, A and the launch
  # distribution are absent from the endpoint.
  for (rho in c(Inf, 2, 4)) {
    xm <- ref_xmax(rho)
    for (k in c(0.4, 1, 3.7)) {
      vals <- c(t_max(0, 1, k, rho), t_max(0.5, 1, k, rho),
                t_max(0.5, 9, k, rho), t_max(0, 0.2, k, rho))
      expect_equal(vals, rep(xm / k, 4), tolerance = 1e-12)
    }
  }
  expect_true(is.infinite(t_max(0.2, 1, 0)))          # LBA limit
  expect_true(is.na(t_max(0.2, 1, 1, 1)))             # rho = 1 has no endpoint
})

test_that("V_c(0) is k b H'(x_max) and scales with b, not with k alone", {
  vc <- EMC2:::bawf_vcrit_vec(c(0, 0.3), c(1, 2), c(1.5, 0.7), 0)
  expect_equal(vc, exp(1) * c(1.5, 0.7) * c(1, 2), tolerance = 1e-12)
  expect_equal(EMC2:::bawf_vcrit_vec(0.3, 2, 0.7, 4), 0.7 * 2 * (4 / 3)^3,
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 2. The k = 0 oracle
# ---------------------------------------------------------------------------

test_that("k = 0 is exactly the LBA", {
  ts <- c(0.05, 0.2, 0.5, 1, 3)
  for (A in c(0, 0.3)) {
    for (pd in c(TRUE, FALSE)) {
      expect_equal(d_bawf(ts, A, 1.0, 2.0, 0.8, 0, launch = 0L, posdrift = pd),
                   EMC2:::dlba(ts, A, 1.0, 2.0, 0.8, posdrift = pd),
                   tolerance = 1e-12)
      expect_equal(p_bawf(ts, A, 1.0, 2.0, 0.8, 0, launch = 0L, posdrift = pd),
                   EMC2:::plba(ts, A, 1.0, 2.0, 0.8, posdrift = pd),
                   tolerance = 1e-12)
    }
  }
  # and the limit is approached, not just hit at exactly zero
  expect_equal(d_bawf(ts, 0.3, 1.0, 2.0, 0.8, 1e-9, launch = 0L),
               EMC2:::dlba(ts, 0.3, 1.0, 2.0, 0.8), tolerance = 1e-7)
})

# ---------------------------------------------------------------------------
# 3. Both launches against the independent references
# ---------------------------------------------------------------------------

bawf_cases <- list(
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, k = 1.2, launch = 1L, rho = Inf),
  list(A = 0.0, b = 1.0, p1 = 1.2, p2 = 0.4, k = 0.8, launch = 1L, rho = Inf),
  list(A = 0.6, b = 1.5, p1 = 4.0, p2 = 1.5, k = 1.0, launch = 0L, rho = Inf),
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, k = 1.2, launch = 1L, rho = 2),
  list(A = 0.4, b = 1.2, p1 = 1.0, p2 = 0.5, k = 1.2, launch = 1L, rho = 4),
  list(A = 0.5, b = 1.1, p1 = 3.0, p2 = 1.0, k = 2.0, launch = 0L, rho = 2),
  list(A = 0.3, b = 1.0, p1 = 0.8, p2 = 0.9, k = 3.0, launch = 1L, rho = Inf)
)

test_that("the CDF matches an independent start-point quadrature", {
  skip_on_cran()
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    us <- Tm * c(0.05, 0.2, 0.5, 0.75, 0.9, 0.99)
    got <- p_bawf(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                  rho = cs$rho)
    ref <- vapply(us, ref_F_bawf, numeric(1), cs$p1, cs$p2, cs$b, cs$A, cs$k,
                  cs$rho, cs$launch)
    expect_equal(got, ref, tolerance = 1e-7)
  }
})

test_that("the density matches an independent start-point quadrature", {
  skip_on_cran()
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    us <- Tm * c(0.05, 0.2, 0.5, 0.75, 0.9, 0.99)
    got <- d_bawf(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                  rho = cs$rho)
    ref <- vapply(us, ref_f_bawf, numeric(1), cs$p1, cs$p2, cs$b, cs$A, cs$k,
                  cs$rho, cs$launch)
    expect_equal(got, ref, tolerance = 1e-7)
  }
})

test_that("integrate(f) reproduces F through the live/frozen seam", {
  skip_on_cran()
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    for (frac in c(0.5, 0.95)) {
      u1 <- Tm * frac
      I <- stats::integrate(function(x) {
        d_bawf(x, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
               rho = cs$rho)
      }, 1e-10, u1, subdivisions = 2000, rel.tol = 1e-11)$value
      expect_equal(I, p_bawf(u1, cs$A, cs$b, cs$p1, cs$p2, cs$k,
                             launch = cs$launch, rho = cs$rho),
                   tolerance = 1e-8)
    }
  }
})

# ---------------------------------------------------------------------------
# 4. The endpoint
# ---------------------------------------------------------------------------

test_that("the CDF is flat past T_max and F(Inf) is the finishing mass", {
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    at <- p_bawf(Tm, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                 rho = cs$rho)
    beyond <- p_bawf(c(Tm * 1.5, Tm * 50, Inf), cs$A, cs$b, cs$p1, cs$p2,
                     cs$k, launch = cs$launch, rho = cs$rho)
    expect_equal(beyond, rep(at, 3), tolerance = 1e-12)
    expect_lt(at, 1)                    # always defective
    expect_equal(d_bawf(c(Tm, Tm * 1.5), cs$A, cs$b, cs$p1, cs$p2, cs$k,
                        launch = cs$launch, rho = cs$rho), c(0, 0))
  }
})

test_that("the density shuts down quadratically with A > 0, linearly at A = 0", {
  k <- 1.3; b <- 1.2; Tm <- 1 / k
  eps <- c(1e-2, 5e-3, 2.5e-3, 1.25e-3)
  ord <- function(A) {
    fv <- d_bawf(Tm - eps, A, b, 1.0, 0.5, k)
    log(fv[-1] / fv[-length(fv)]) / log(eps[-1] / eps[-length(eps)])
  }
  expect_equal(ord(0), rep(1, 3), tolerance = 0.02)
  expect_equal(ord(0.4), rep(2, 3), tolerance = 0.02)
})

# ---------------------------------------------------------------------------
# 5. Monte Carlo
# ---------------------------------------------------------------------------

test_that("Monte Carlo agrees for both launch distributions", {
  skip_on_cran()
  set.seed(42)
  N <- 1000
  for (cs in bawf_cases[c(1, 3, 4)]) {
    V <- if (cs$launch == 1L) rlnorm(N, cs$p1, cs$p2)
         else msm::rtnorm(N, cs$p1, cs$p2, lower = 0)
    z <- runif(N, 0, cs$A)
    hit <- vapply(seq_len(N), function(i)
      ref_fp_bawf(V[i], z[i], cs$b, cs$k, cs$rho), numeric(1))
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    us <- Tm * c(0.25, 0.5, 0.75, 0.95)
    emp <- vapply(us, function(u) mean(hit <= u), numeric(1))
    got <- p_bawf(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                  rho = cs$rho)
    expect_lt(max(abs(got - emp)), 4 * sqrt(0.25 / N) + 1e-3)
    expect_lt(abs(mean(is.finite(hit)) -
                    p_bawf(Inf, cs$A, cs$b, cs$p1, cs$p2, cs$k,
                           launch = cs$launch, rho = cs$rho)),
              4 * sqrt(0.25 / N) + 1e-3)
  }
})

# ---------------------------------------------------------------------------
# 6. Numerical hygiene
# ---------------------------------------------------------------------------

test_that("the log output equals the log of the natural output", {
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    us <- Tm * c(0.01, 0.3, 0.6, 0.95, 0.999)
    for (fn in list(d_bawf, p_bawf)) {
      nat <- fn(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                rho = cs$rho)
      lg <- fn(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
               rho = cs$rho, log_out = TRUE)
      keep <- nat > 0
      expect_equal(lg[keep], log(nat[keep]), tolerance = 1e-9)
    }
  }
})

test_that("the CDF is monotone and the density non-negative across the seams", {
  for (cs in bawf_cases) {
    Tm <- t_max(cs$A, cs$b, cs$k, cs$rho)
    # Probe the support seams and representative interior points.  A dense
    # development sweep does not add package-level regression coverage here.
    us <- Tm * c(1e-4, 0.01, 0.1, 0.5, 0.9, 0.99, 0.9999)
    Fv <- p_bawf(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                 rho = cs$rho)
    fv <- d_bawf(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, launch = cs$launch,
                 rho = cs$rho)
    expect_false(anyNA(c(Fv, fv)))
    expect_true(all(diff(Fv) >= -1e-12))
    expect_true(all(fv >= 0))
  }
})

test_that("boundary parameters produce no NaN", {
  grid <- expand.grid(A = c(0, 1e-12, 0.5), k = c(0, 1e-12, 1e3),
                      p2 = c(1e-4, 3), u = c(1e-10, 0.3, 1e6, Inf))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    for (launch in c(0L, 1L)) {
      expect_false(is.na(d_bawf(g$u, g$A, 1, 1, g$p2, g$k, launch = launch)))
      expect_false(is.na(p_bawf(g$u, g$A, 1, 1, g$p2, g$k, launch = launch)))
    }
  }
})

# ---------------------------------------------------------------------------
# 7. The constructor
# ---------------------------------------------------------------------------

test_that("the constructor wires drift_distribution and rho consistently", {
  expect_equal(BAwF()$c_name, "BAwF_LOGN")
  expect_equal(names(BAwF()$p_types)[1:6],
               c("mean", "cv", "B", "A", "t0", "k"))
  expect_equal(BAwF("normal")$c_name, "BAwF")
  expect_equal(names(BAwF("normal")$p_types)[1:6],
               c("v", "sv", "B", "A", "t0", "k"))
  expect_equal(BAwF("normal", posdrift = FALSE)$c_name, "BAwFIO")
  expect_equal(BAwF(rho = 2)$c_name, "BAwF_LOGN_RHO2")
  expect_equal(BAwF("normal", rho = 4)$c_name, "BAwF_RHO4")
  expect_error(BAwF(rho = 1), "no finite endpoint")
  expect_error(BAwF(rho = 3), "must be one of")
  expect_error(BAwF("lognormal", posdrift = FALSE), "posdrift only applies")
  # k = 0 must stay exactly reachable: it is the LBA limit.
  expect_equal(BAwF()$bound$exception[["k"]], 0)
})

test_that("Ttransform reports the endpoint and the omission boundary", {
  m <- BAwF()
  pars <- cbind(mu = c(1, 1), sigma = c(0.5, 0.5), B = c(0.9, 1.8),
                A = c(0.3, 0.3), t0 = c(0.2, 0.2), k = c(1.25, 1.25))
  out <- m$Ttransform(pars, NULL)
  expect_equal(out[, "b"], pars[, "B"] + pars[, "A"])
  # Doubling B leaves T_max alone and moves V_c(0) -- the identification claim.
  expect_equal(out[, "Tmax"], rep(1 / 1.25, 2), tolerance = 1e-12)
  expect_equal(out[, "rt_max"], 0.2 + out[, "Tmax"])
  expect_equal(out[, "Vcrit"], exp(1) * 1.25 * out[, "b"], tolerance = 1e-12)
})

test_that("dfun/pfun use the same launch distribution as the c_name", {
  pars <- cbind(mean = 1, cv = 0.5, v = 1, sv = 0.5, b = 1.2, A = 0.3,
                t0 = 0.15, k = 1.1)
  rt <- 0.4
  expect_equal(BAwF()$dfun(rt, pars),
               d_bawf(rt - 0.15, 0.3, 1.2, 1, 0.5, 1.1, launch = 1L))
  expect_equal(BAwF("normal")$dfun(rt, pars),
               d_bawf(rt - 0.15, 0.3, 1.2, 1, 0.5, 1.1, launch = 0L))
  expect_equal(BAwF(rho = 2)$pfun(rt, pars),
               p_bawf(rt - 0.15, 0.3, 1.2, 1, 0.5, 1.1, launch = 1L, rho = 2))
})

# ---------------------------------------------------------------------------
# 8. The compiled race likelihood
# ---------------------------------------------------------------------------

ref_race_ll_bawf <- function(dadm, pars, launch, rho = Inf,
                             min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    ref_F_bawf(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
               pars[i, "k"], rho, launch)
  }
  fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    ref_f_bawf(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
               pars[i, "k"], rho, launch)
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

bawf_mk <- function(model, form, consts, dat) {
  des <- design(data = dat, model = model, matchfun = function(d) d$S == d$lR,
                formula = form, constants = consts)
  list(des = des, emc = suppressMessages(
    make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
             rt_resolution = NULL)))
}

bawf_ll <- function(fx, p, p_types_override = NULL) {
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

bawf_dat <- function(n = 60) {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat[seq_len(n), ]
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  dat <- bawf_dat()
  p <- c(mu = 0.9, sigma = log(0.6), v = 3, sv = log(1), B = log(0.8),
         A = log(0.3), t0 = log(0.15), k = log(0.8))
  fx <- list(
    ln = suppressMessages(bawf_mk(BAwF,
      list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), NULL, dat)),
    no = suppressMessages(bawf_mk(function() BAwF("normal"),
      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), c(sv = log(1)), dat)))
  for (which in c("ln", "no")) {
    launch <- if (which == "ln") 1L else 0L
    pp <- p[names(sampled_pars(fx[[which]]$des))]
    got <- bawf_ll(fx[[which]], pp)
    expect_true(is.finite(got))
    dadm <- fx[[which]]$emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(pp, dadm, fx[[which]]$emc[[1]]$model())
    expect_equal(got, ref_race_ll_bawf(dadm, pars, launch), tolerance = 1e-6)
  }
})

test_that("the compiled BAwF likelihood reproduces the LBA at k = 0", {
  # A second, independent route to the k = 0 oracle: two models, two adapter
  # branches, two column layouts, one answer.  This is what catches a wiring
  # error the kernel tests cannot see (wrong column order, wrong t0_index, a
  # stray posdrift default).
  skip_on_cran()
  dat <- bawf_dat()
  fl <- suppressMessages(bawf_mk(LBA, list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                                 c(sv = log(1)), dat))
  ff <- suppressMessages(bawf_mk(function() BAwF("normal"),
                                 list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                                 c(sv = log(1), k = log(0)), dat))
  pv <- c(v = 3, B = log(0.9), A = log(0.3), t0 = log(0.15))
  l_lba <- bawf_ll(fl, pv[names(sampled_pars(fl$des))])
  l_bawf <- bawf_ll(ff, pv[names(sampled_pars(ff$des))])
  expect_true(is.finite(l_lba))
  expect_equal(l_bawf, l_lba, tolerance = 1e-8)
})

test_that("a p_types reordering is caught by the column contract", {
  skip_on_cran()
  dat <- bawf_dat()
  fx <- suppressMessages(bawf_mk(BAwF,
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), NULL, dat))
  p <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(0.8))
  swapped <- names(BAwF()$p_types)
  swapped[1:2] <- swapped[2:1]
  expect_error(bawf_ll(fx, p[names(sampled_pars(fx$des))],
                       p_types_override = swapped),
               "BAwF_LOGN kernels expect parameter column")
})

test_that("omissions past t0 + T_max stay well posed", {
  skip_on_cran()
  dat <- bawf_dat(80)
  dat$rt[1:8] <- Inf                 # intrinsic omissions
  dat$R[1:8] <- NA                   # with an unknown response
  fx <- suppressMessages(bawf_mk(BAwF,
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), NULL, dat))
  p <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(0.8))
  got <- bawf_ll(fx, p[names(sampled_pars(fx$des))])
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p[names(sampled_pars(fx$des))], dadm,
                                    fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll_bawf(dadm, pars, 1L), tolerance = 1e-6)
  # A finite response beyond t0 + T_max has exactly zero density and must
  # floor rather than return -Inf or NaN.
  dat2 <- dat
  dat2$rt[1:8] <- 0.15 + 1 / 0.8 + 5
  fx2 <- suppressMessages(bawf_mk(BAwF,
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1), NULL, dat2))
  expect_true(is.finite(bawf_ll(fx2, p[names(sampled_pars(fx2$des))])))
})

# ---------------------------------------------------------------------------
# 9. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree with each other and with the CDF", {
  skip_on_cran()
  pars <- cbind(mu = 1.0, sigma = 0.5, b = 1.2, A = 0.4, t0 = 0.15, k = 1.2)
  n <- 1000
  pm <- pars[rep(1, 2 * n), , drop = FALSE]
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  set.seed(7)
  cpp <- EMC2:::.rfun_BAwF(lR, pm, rep(TRUE, 2 * n), launch = 1L)
  set.seed(7)
  withr::with_options(list(emc2.cpp_rfun = FALSE), {
    rr <- EMC2:::.rfun_BAwF(lR, pm, rep(TRUE, 2 * n), launch = 1L)
  })
  # The two paths consume the RNG stream differently (rlnorm vs
  # exp(mu + sigma * norm_rand)), so only the distributions may be compared.
  fin <- is.finite(cpp$rt)
  expect_lt(abs(mean(fin) - mean(is.finite(rr$rt))), 0.06)
  # The winner's RT is the min of two iid first passages, so its CDF is
  # 1 - (1 - F)^2 with F the single-accumulator CDF.
  qs <- stats::quantile(cpp$rt[fin], c(0.25, 0.5, 0.75))
  emp <- vapply(qs, function(q) mean(cpp$rt <= q, na.rm = TRUE), numeric(1))
  Fq <- p_bawf(as.numeric(qs) - 0.15, 0.4, 1.2, 1.0, 0.5, 1.2, launch = 1L)
  expect_lt(max(abs(as.numeric(emp) - (1 - (1 - Fq)^2))), 0.06)
  expect_true(all(cpp$rt[fin] <= 0.15 + 1 / 1.2 + 1e-8))   # the hard endpoint
})

test_that("make_data produces omissions the design can be fit back through", {
  skip_on_cran()
  dat <- bawf_dat(40)
  des <- suppressMessages(design(
    data = dat, model = BAwF, matchfun = function(d) d$S == d$lR,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1)))
  p <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(0.8))
  set.seed(3)
  sim <- make_data(p, design = des, n_trials = 60)
  expect_true(any(is.infinite(sim$rt)))          # intrinsic omissions
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] <= 0.15 + 1 / 0.8 + 1e-8))
  emc <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(bawf_ll(list(emc = emc, des = des), p)))
})
