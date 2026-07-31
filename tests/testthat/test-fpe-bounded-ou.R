# Bounded OU (Smith & Ratcliff 2004): the DDM with leak, solved by the
# Fokker-Planck solver with an absorbing boundary at BOTH ends of the domain.
#
# The load-bearing test is the first one: at beta = 0 the model IS the Wiener
# diffusion, so the package's own Navarro-Fuss DDM is an exact oracle -- and the
# bounded model reaches it down a completely different code path (PDE march vs
# series expansion), so agreement is real evidence rather than a tautology.

fpe_bou <- EMC2:::fpe_bou_fht_pdf_cdf_vec

ddm_oracle <- function(rt, resp, v, a, Z, s = 1) {
  R <- factor(rep(resp, length(rt)), levels = c("lower", "upper"))
  pars <- cbind(a = rep(a, length(rt)), v = v, t0 = 0, s = s,
                Z = Z, SZ = 0, sv = 0, st0 = 0)
  list(d = EMC2:::dDDM(rt, R, pars, precision = 1e-12),
       p = EMC2:::pDDM(rt, R, pars, precision = 1e-12))
}

test_that("beta = 0 reproduces the Wiener DDM on both boundaries", {
  tt <- seq(0.05, 2.5, by = 0.02)
  grid <- expand.grid(v = c(-1.5, 0, 1, 2.5), Z = c(0.3, 0.5, 0.7), a = c(0.8, 1.5))

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    r <- fpe_bou(tt, v = g$v, beta = 0, a = g$a, Z = g$Z, sigma = 1,
                 nx = 512, nt = 6000, grade = 1, tgrade = 32)
    up <- ddm_oracle(tt, "upper", g$v, g$a, g$Z)
    lo <- ddm_oracle(tt, "lower", g$v, g$a, g$Z)

    lab <- sprintf("v=%g Z=%g a=%g", g$v, g$Z, g$a)
    # Defective cdfs are the tighter check: they are integrals, so they do not
    # inherit the density's sensitivity on the rising flank.
    expect_lt(max(abs(r$cdf_upper - up$p)), 2e-3, label = paste("cdf upper", lab))
    expect_lt(max(abs(r$cdf_lower - lo$p)), 2e-3, label = paste("cdf lower", lab))

    # Densities, judged where a likelihood actually lives (central 99.5% of each
    # response's mass).  Outside that window the density is < 1e-7 of its peak
    # and its relative error says nothing about the fit.
    for (side in c("upper", "lower")) {
      ref <- if (side == "upper") up$d else lo$d
      got <- if (side == "upper") r$pdf_upper else r$pdf_lower
      cm <- cumsum(ref); if (max(cm) <= 0) next
      cm <- cm / max(cm)
      sel <- cm > 0.0025 & cm < 0.9975 & ref > 1e-6
      if (sum(sel) < 5) next
      expect_lt(max(abs(log(got[sel]) - log(ref[sel]))), 2e-2,
                label = paste("log pdf", side, lab))
    }
  }
})

test_that("the two boundaries are resolved alike (symmetry)", {
  # v = 0 and Z = 0.5 makes the problem exactly symmetric, so any difference
  # between the two response densities is pure discretisation asymmetry.  This
  # is what catches a mesh graded toward one end only.
  tt <- seq(0.05, 2.5, by = 0.02)
  for (b in c(0, 2, 4, 8)) {
    for (gr in c(1, 8)) {
      r <- fpe_bou(tt, v = 0, beta = b, a = 1, Z = 0.5, sigma = 1,
                   nx = 384, nt = 3000, grade = gr, tgrade = 32)
      expect_lt(max(abs(r$pdf_upper - r$pdf_lower)), 1e-10,
                label = sprintf("pdf symmetry beta=%g grade=%g", b, gr))
      expect_lt(max(abs(r$cdf_upper - r$cdf_lower)), 1e-10,
                label = sprintf("cdf symmetry beta=%g grade=%g", b, gr))
    }
  }
})

test_that("the two defective cdfs and the survivor account for all the mass", {
  # cdf_upper is formed as (1 - mass) - cdf_lower precisely so this identity is
  # exact rather than approximate: only the SPLIT between the two responses
  # carries quadrature error, never the total.
  tt <- seq(0.05, 3.0, by = 0.02)
  for (b in c(0, 4, 8)) {
    r <- fpe_bou(tt, v = 1.2, beta = b, a = 1.2, Z = 0.4, sigma = 1,
                 nx = 384, nt = 3000, grade = 1, tgrade = 32)
    expect_lt(max(abs(r$cdf_upper + r$cdf_lower + r$surv - 1)), 1e-12)
    expect_lt(r$mismatch, 1e-6)
    expect_false(any(diff(r$cdf_upper) < 0))
    expect_false(any(diff(r$cdf_lower) < 0))
    expect_true(all(r$pdf_upper >= 0) && all(r$pdf_lower >= 0))
  }
})

test_that("leak changes the model", {
  # Guards against a silent no-op: if beta were dropped somewhere in the affine
  # drift the tests above would all still pass.
  tt <- seq(0.05, 3.0, by = 0.02)
  r0 <- fpe_bou(tt, v = 1.2, beta = 0, a = 1.2, Z = 0.4, sigma = 1,
                nx = 384, nt = 3000, grade = 1, tgrade = 32)
  r8 <- fpe_bou(tt, v = 1.2, beta = 8, a = 1.2, Z = 0.4, sigma = 1,
                nx = 384, nt = 3000, grade = 1, tgrade = 32)
  # Decay toward a start point below the midpoint pulls mass away from the
  # upper barrier, so the upper response gets less likely.
  expect_lt(tail(r8$cdf_upper, 1), tail(r0$cdf_upper, 1) - 0.05)
  expect_gt(max(abs(r0$pdf_upper - r8$pdf_upper)), 0.1)
})

# ---------------------------------------------------------------------------
# Across-trial variability, and the simulator.
#
# These exercise fpe_bou.h (quadrature + solve cache) and bou_diffusion.cpp
# rather than the solver core above.
# ---------------------------------------------------------------------------

test_that("across-trial variability matches the DDM oracle at beta = 0", {
  tt <- seq(0.25, 2.0, by = 0.05)
  Z <- 0.45; a <- 1.2; v <- 1.5; t0 <- 0.15

  # NB the SZ contract differs between the two sides and it is easy to get
  # wrong: bou_pdf_cdf_vec takes SZ RAW and widens it internally (as
  # d_DDM_Wien_raw does, model_DDM.h:122), whereas the R-level dDDM()/pDDM()
  # expect the ALREADY-widened value because Ttransform normally ran first.
  widen <- function(SZ) 2 * SZ * min(Z, 1 - Z)

  cases <- list(c(sv = 0, SZ = 0, st0 = 0), c(sv = 1.0, SZ = 0, st0 = 0),
                c(sv = 0, SZ = 0.3, st0 = 0), c(sv = 0, SZ = 0, st0 = 0.1),
                c(sv = 1.0, SZ = 0.3, st0 = 0.1))
  for (cs in cases) {
    for (side in c("lower", "upper")) {
      Rf <- factor(rep(side, length(tt)), levels = c("lower", "upper"))
      Ri <- rep(if (side == "upper") 2L else 1L, length(tt))
      pm <- cbind(a = rep(a, length(tt)), v = v, t0 = t0, s = 1, Z = Z,
                  SZ = widen(cs[["SZ"]]), sv = cs[["sv"]], st0 = cs[["st0"]])
      b <- EMC2:::bou_pdf_cdf_vec(tt, Ri, v, a, Z, cs[["sv"]], cs[["SZ"]],
                                  t0, cs[["st0"]], 1, 0,
                                  nx = 512, dt_target = 2.5e-4, n_st0 = 15)
      lab <- sprintf("sv=%g SZ=%g st0=%g %s", cs[["sv"]], cs[["SZ"]],
                     cs[["st0"]], side)
      expect_lt(max(abs(b$pdf - EMC2:::dDDM(tt, Rf, pm, precision = 1e-10))),
                1e-3, label = paste("pdf", lab))
      expect_lt(max(abs(b$cdf - EMC2:::pDDM(tt, Rf, pm, precision = 1e-10))),
                1e-3, label = paste("cdf", lab))
    }
  }
})

test_that("st0 uses the DDM's lower-edge convention", {
  # t0 ~ U(t0, t0+st0), NOT U(t0-st0/2, t0+st0/2).  Reconstruct an st0 > 0
  # density from point-t0 densities under each convention: the wrong one is off
  # by ~0.34 where the right one agrees to ~2e-3, so this cannot pass by luck.
  tt <- seq(0.3, 1.5, by = 0.05)
  Rf <- factor(rep("upper", length(tt)), levels = c("lower", "upper"))
  d <- function(t0, st0) EMC2:::dDDM(tt, Rf,
    cbind(a = rep(1.2, length(tt)), v = 1.5, t0 = t0, s = 1, Z = .45,
          SZ = 0, sv = 0, st0 = st0), precision = 1e-10)
  ref <- d(0.15, 0.1)
  u <- seq(0, 1, length.out = 201)
  lower <- rowMeans(sapply(u, function(x) d(0.15 + 0.1 * x, 0)))
  centred <- rowMeans(sapply(u - 0.5, function(x) d(0.15 + 0.1 * x, 0)))
  expect_lt(max(abs(ref - lower)), 1e-2)
  expect_gt(max(abs(ref - centred)), 1e-1)
})

test_that("the quadrature costs exactly n_sv * n_sz solves", {
  # Guards the solve cache against two regressions that are invisible in the
  # density: a horizon that grows per row (which would re-solve every node as a
  # sorted rt vector is walked), and an st0 loop that triggers solves instead of
  # reusing them.
  tt <- seq(0.2, 2.0, by = 0.02)
  Ri <- rep(2L, length(tt))
  none <- EMC2:::bou_pdf_cdf_vec(tt, Ri, 1.5, 1.2, .45, 0, 0, 0.1, 0, 1, 0)
  expect_equal(none$n_solves, 1)
  # st0 alone shifts the query time only: no extra solves.
  st0only <- EMC2:::bou_pdf_cdf_vec(tt, Ri, 1.5, 1.2, .45, 0, 0, 0.1, 0.1, 1, 0)
  expect_equal(st0only$n_solves, 1)
  both <- EMC2:::bou_pdf_cdf_vec(tt, Ri, 1.5, 1.2, .45, 1, 0.3, 0.1, 0.1, 1, 0,
                                 n_sv = 5, n_sz = 3)
  expect_equal(both$n_solves, 15)
})

test_that("the simulator agrees with the solver", {
  set.seed(4)
  for (beta in c(0, 5)) {
    N <- 20000
    sim <- EMC2:::rbou_cpp(N, 1.2, 1.0, 0.5, 0, 0, 0.1, 0, 1, beta,
                           dt = 1e-4, t_max = 30)
    sim <- sim[is.finite(sim$rt) & !is.na(sim$R), ]
    p_sim <- mean(sim$R == 2)
    se <- sqrt(p_sim * (1 - p_sim) / nrow(sim))
    tg <- seq(0.1005, 6, by = 0.005)
    Fg <- EMC2:::bou_pdf_cdf_vec(tg, rep(2L, length(tg)), 1.2, 1.0, 0.5,
                                 0, 0, 0.1, 0, 1, beta,
                                 nx = 384, dt_target = 1e-3)$cdf
    expect_lt(abs(p_sim - max(Fg)), 4 * se)
    qsim <- unname(quantile(sim$rt[sim$R == 2], c(.1, .5, .9)))
    qsol <- suppressWarnings(approx(Fg / max(Fg), tg, xout = c(.1, .5, .9))$y)
    expect_lt(max(abs(qsim - qsol)), 0.02)
  }
})

test_that("BOU reaches the C++ DDM likelihood path and beta = 0 gives the DDM", {
  # End-to-end through design() -> make_emc() -> calc_ll_manager, i.e. the same
  # route the sampler takes.  This is what checks the wiring rather than the
  # kernels: column order against emc2col::bou::spec(), the DDMAdapter swap in
  # c_log_likelihood_DDM_pt, and the per-particle cache reset.
  dat <- droplevels(forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ])
  f <- list(v ~ 0 + S, a ~ E, t0 ~ 1, s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1)

  d_ddm <- design(data = dat, model = DDM, formula = f, constants = c(s = log(1)))
  d_bou <- design(data = dat, model = BOU, formula = c(f, list(beta ~ 1)),
                  constants = c(s = log(1), beta = log(0)))
  e_ddm <- make_emc(dat, d_ddm, compress = FALSE, rt_resolution = NULL, type = "single")
  e_bou <- make_emc(dat, d_bou, compress = FALSE, rt_resolution = NULL, type = "single")

  pn <- names(sampled_pars(d_ddm))
  set.seed(99)
  P <- matrix(rnorm(4 * length(pn), 0, 0.25), ncol = length(pn),
              dimnames = list(NULL, pn))
  P[, "t0"] <- log(0.15); P[, "a"] <- log(1.2); P[, "sv"] <- log(0.6)
  P[, "SZ"] <- qnorm(0.2); P[, "st0"] <- log(0.05)

  l_ddm <- EMC2:::calc_ll_manager(P, e_ddm[[1]]$data[[1]], e_ddm[[1]]$model)
  l_bou <- EMC2:::calc_ll_manager(P, e_bou[[1]]$data[[1]], e_bou[[1]]$model)
  expect_true(all(is.finite(l_bou)))
  # ~1e-5 per trial over ~800 trials; the solver's error, not a wiring fault.
  expect_lt(max(abs(l_ddm - l_bou)), 0.05)

  # And the leak must be live end to end, not silently dropped by the column
  # mapping -- which is exactly the failure a positional column order invites.
  d_b <- design(data = dat, model = BOU, formula = c(f, list(beta ~ 1)),
                constants = c(s = log(1)))
  e_b <- make_emc(dat, d_b, compress = FALSE, rt_resolution = NULL, type = "single")
  pn2 <- names(sampled_pars(d_b))
  P2 <- matrix(0, nrow = 2, ncol = length(pn2), dimnames = list(NULL, pn2))
  for (j in pn) P2[, j] <- P[1:2, j]
  P2[, "beta"] <- log(1e-12)
  l0 <- EMC2:::calc_ll_manager(P2, e_b[[1]]$data[[1]], e_b[[1]]$model)
  P2[, "beta"] <- log(4)
  l4 <- EMC2:::calc_ll_manager(P2, e_b[[1]]$data[[1]], e_b[[1]]$model)
  expect_gt(min(abs(l0 - l4)), 1)
})

# ---------------------------------------------------------------------------
# Collapsing bounds.
#
# Both barriers close on the midpoint a/2, so the separation runs from a down to
# aInf.  There is no closed form for a diffusion between two converging
# boundaries, so the evidence here is of three independent kinds: a degenerate
# case that must reduce EXACTLY to the fixed-bound model (which is itself
# validated against Navarro-Fuss), agreement between two independent
# implementations of the time-varying march, and agreement with a simulator that
# shares none of the solver's numerics.
# ---------------------------------------------------------------------------

test_that("a collapse to the starting separation is the fixed-bound model", {
  # Every collapse form degenerates when aInf == a, and FPE_Boundary::set_kind
  # reports fixed = true for it, which restores the one-time factorisation.  If
  # the moving-frame term were being added unconditionally, or the midpoint
  # geometry disagreed with the old (0, a) one by a rounding step, this would
  # not be bit-exact -- and it is the test that anchors the collapse code to the
  # already-validated fixed-bound model.
  tt <- seq(0.05, 3.0, by = 0.02)
  r_fix <- fpe_bou(tt, v = 1.0, beta = 2, a = 1.2, Z = 0.45, sigma = 1,
                   nx = 384, nt = 3000, grade = 1, tgrade = 32)
  for (bk in 1:4) {
    r_deg <- fpe_bou(tt, v = 1.0, beta = 2, a = 1.2, Z = 0.45, sigma = 1,
                     bkind = bk, aInf = 1.2, tau = 0.5, pw = 1.5,
                     nx = 384, nt = 3000, grade = 1, tgrade = 32)
    expect_identical(r_deg$pdf_upper, r_fix$pdf_upper)
    expect_identical(r_deg$pdf_lower, r_fix$pdf_lower)
  }
})

test_that("collapsing bounds keep the mass identity and monotone cdfs", {
  tt <- seq(0.05, 3.0, by = 0.02)
  for (bk in 1:4) {
    r <- fpe_bou(tt, v = 1.0, beta = 2, a = 1.2, Z = 0.45, sigma = 1,
                 bkind = bk, aInf = 0.3, tau = 0.6, pw = 1.5,
                 nx = 384, nt = 3000, grade = 1, tgrade = 32)
    expect_lt(max(abs(r$cdf_upper + r$cdf_lower + r$surv - 1)), 1e-12)
    expect_lt(r$mismatch, 1e-6)
    expect_false(any(diff(r$cdf_upper) < 0))
    expect_false(any(diff(r$cdf_lower) < 0))
    expect_true(all(r$pdf_upper >= 0) && all(r$pdf_lower >= 0))
  }
})

test_that("a symmetric collapse stays symmetric", {
  # v = 0, Z = 0.5 and barriers that close on the midpoint: the problem is
  # symmetric at every t, so the two densities must agree.  This is what catches
  # a collapse applied to one barrier only, which is the easiest way to get the
  # midpoint geometry wrong.
  tt <- seq(0.05, 3.0, by = 0.02)
  for (bk in 1:4) {
    r <- fpe_bou(tt, v = 0, beta = 3, a = 1, Z = 0.5, sigma = 1,
                 bkind = bk, aInf = 0.2, tau = 0.5, pw = 2,
                 nx = 384, nt = 3000, grade = 1, tgrade = 32)
    expect_lt(max(abs(r$pdf_upper - r$pdf_lower)), 1e-10,
              label = sprintf("collapsing symmetry bkind=%d", bk))
  }
})

test_that("the lane batch and the scalar march agree under collapse", {
  # The two paths implement the moving operator INDEPENDENTLY: fpe_solve rebuilds
  # an FPE_Op and refactorises a FPE_Tri, while the batch runs the interleaved
  # bou_build_op_lanes and its own Thomas elimination.  With fixed bounds both
  # build the operator once, so this comparison only becomes load-bearing once
  # it is rebuilt every step -- which is exactly the code the collapse adds.
  # They use different horizons, hence different time schedules, so this is
  # agreement to solver tolerance rather than bit equality.
  #
  # linear_additive gets its own, much looser tolerance, and it is not a fudge:
  # its b'(t) jumps at t = tau, Crank-Nicolson loses an order at the kink, and
  # the two schedules straddle it differently. The measured gap is the point --
  # 5e-5 for the three smooth forms against 6e-3 for the kinked one. Holding all
  # four to the loose bound would stop the smooth cases testing anything.
  tol <- c(5e-4, 5e-4, 2e-2, 5e-4)   # weibull, exponential, lin_add, lin_mult
  tt <- seq(0.1, 2.5, by = 0.02)
  for (bk in 1:4) {
    rs <- fpe_bou(tt, v = 1.0, beta = 2, a = 1.2, Z = 0.45, sigma = 1,
                  bkind = bk, aInf = 0.3, tau = 0.6, pw = 1.5,
                  nx = 384, nt = 6000, grade = 1, tgrade = 32)
    up <- EMC2:::bou_pdf_cdf_vec(tt, rep(2L, length(tt)), 1.0, 1.2, 0.45,
                                 0, 0, 0, 0, 1, 2,
                                 bkind = bk, aInf = 0.3, tau = 0.6, pw = 1.5,
                                 nx = 384, dt_target = 1 / 2000, tgrade = 32)
    lo <- EMC2:::bou_pdf_cdf_vec(tt, rep(1L, length(tt)), 1.0, 1.2, 0.45,
                                 0, 0, 0, 0, 1, 2,
                                 bkind = bk, aInf = 0.3, tau = 0.6, pw = 1.5,
                                 nx = 384, dt_target = 1 / 2000, tgrade = 32)
    ok <- rs$pdf_upper > 1e-4
    expect_lt(max(abs(up$pdf[ok] / rs$pdf_upper[ok] - 1)), tol[bk],
              label = sprintf("lane vs scalar, upper, bkind=%d", bk))
    expect_lt(max(abs(lo$pdf[ok] / rs$pdf_lower[ok] - 1)), tol[bk],
              label = sprintf("lane vs scalar, lower, bkind=%d", bk))
  }
})

test_that("collapsing bounds finish the race", {
  tt <- seq(0.05, 3.0, by = 0.02)
  r0 <- fpe_bou(tt, v = 0.8, beta = 1, a = 1.2, Z = 0.5, sigma = 1,
                nx = 384, nt = 3000, grade = 1, tgrade = 32)
  rc <- fpe_bou(tt, v = 0.8, beta = 1, a = 1.2, Z = 0.5, sigma = 1,
                bkind = 2L, aInf = 0.1, tau = 0.4,
                nx = 384, nt = 3000, grade = 1, tgrade = 32)
  # Bringing the barriers in can only absorb mass sooner, at every t.
  expect_true(all(rc$surv <= r0$surv + 1e-12))
  # And with the separation driven down to the floor, essentially nothing is
  # still running -- the collapse is a deadline, which is the point of it.
  expect_lt(tail(rc$surv, 1), 1e-6)
  expect_gt(tail(r0$surv, 1), tail(rc$surv, 1))
})

test_that("the simulator agrees with the solver under collapse", {
  # The simulator steps the exact OU transition against linearly interpolated
  # barriers and shares none of the PDE's numerics, so it is the only
  # independent check available for a moving boundary.  It is used as a check,
  # never as a calibration target.
  set.seed(11)
  N <- 40000
  for (bk in c(2L, 3L)) {
    sim <- EMC2:::rbou_cpp(N, 1.0, 1.2, 0.45, 0, 0, 0.1, 0, 1, 2,
                           bkind = bk, aInf = 0.3, tau = 0.6,
                           dt = 5e-5, t_max = 30)
    sim <- sim[is.finite(sim$rt) & !is.na(sim$R), ]
    p_sim <- mean(sim$R == 2)
    se <- sqrt(p_sim * (1 - p_sim) / nrow(sim))
    tg <- seq(0.1005, 5, by = 0.005)
    Fg <- EMC2:::bou_pdf_cdf_vec(tg, rep(2L, length(tg)), 1.0, 1.2, 0.45,
                                 0, 0, 0.1, 0, 1, 2,
                                 bkind = bk, aInf = 0.3, tau = 0.6,
                                 nx = 384, dt_target = 5e-4, tgrade = 32)$cdf
    expect_lt(abs(p_sim - max(Fg)), 4 * se)
    qsim <- unname(quantile(sim$rt[sim$R == 2], c(.1, .5, .9)))
    qsol <- suppressWarnings(approx(Fg / max(Fg), tg, xout = c(.1, .5, .9))$y)
    expect_lt(max(abs(qsim - qsol)), 0.02)
  }
})

test_that("BOU(boundary_collapse=) reaches the C++ likelihood with the collapse live", {
  # The wiring the collapse adds on top of the fixed-bound path: the c_name
  # suffix that carries the FORM, the optional aInf/tau/pw columns declared
  # after the nine required ones, and bou_row's gate on ContextForDDMModels'
  # bnd_kind.  Getting any of those wrong shows up as an error, a floored
  # likelihood, or -- worst and quietest -- a collapse that does nothing.
  dat <- droplevels(forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ])
  f <- list(v ~ 0 + S, a ~ E, t0 ~ 1, s ~ 1, Z ~ 1, beta ~ 1,
            aInf ~ 1, tau ~ 1)
  des <- design(data = dat,
                model = function() BOU(boundary_collapse = "exponential"),
                formula = f, constants = c(s = log(1)))
  emc <- make_emc(dat, des, compress = FALSE, rt_resolution = NULL,
                  type = "single")
  pn <- names(sampled_pars(des))
  P <- matrix(0, nrow = 2, ncol = length(pn), dimnames = list(NULL, pn))
  P[, "t0"] <- log(0.15); P[, "a"] <- log(1.5); P[, "beta"] <- log(1)
  P[, "tau"] <- log(0.5)
  P[, grep("^v", pn)] <- 1

  P[, "aInf"] <- log(1.5)                     # degenerate: no collapse
  l_none <- EMC2:::calc_ll_manager(P, emc[[1]]$data[[1]], emc[[1]]$model)
  P[, "aInf"] <- log(0.4)                     # a real collapse
  l_coll <- EMC2:::calc_ll_manager(P, emc[[1]]$data[[1]], emc[[1]]$model)
  expect_true(all(is.finite(l_none)) && all(is.finite(l_coll)))
  expect_gt(min(abs(l_none - l_coll)), 1)

  # tau must be live too -- it is the one collapse parameter a wrong column
  # offset would silently swap with aInf.
  P[, "tau"] <- log(0.05)
  l_fast <- EMC2:::calc_ll_manager(P, emc[[1]]$data[[1]], emc[[1]]$model)
  expect_gt(min(abs(l_fast - l_coll)), 1)

  # The fixed-bound generator must be untouched by any of this.
  des0 <- design(data = dat, model = BOU,
                 formula = list(v ~ 0 + S, a ~ E, t0 ~ 1, s ~ 1, Z ~ 1,
                                beta ~ 1),
                 constants = c(s = log(1)))
  expect_false("aInf" %in% names(sampled_pars(des0)))
  expect_identical(BOU()$c_name, "BOU")
  expect_identical(BOU(boundary_collapse = "weibull")$c_name, "BOU_BWEIB")
})
