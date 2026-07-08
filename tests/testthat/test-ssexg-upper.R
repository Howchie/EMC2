# Regression test for wiring the `upper` argument into the plain-EXG
# stop-success integral (ss_exg_stop_success_lpdf). Previously `upper` was
# accepted but ignored: the function always integrated the fixed heuristic
# window [muS - halfw, muS + halfw]. It now caps the upper edge at `upper`.

# pars layout (see test-stop_success_gl.R): go muG=0 sigG=1 tauG=2;
# stop muS=3 sigS=4 tauS=5; (6,7 unused) lbG=8 lbS=9.
make_texg_pars <- function(muG, sigG, tauG, lbG, muS, sigS, tauS, lbS) {
  matrix(c(muG, sigG, tauG, muS, sigS, tauS, 0, 0, lbG, lbS),
         nrow = 1, ncol = 10)
}

SSD <- 0.15
muG <- 0.50; sigG <- 0.05; tauG <- 0.08
muS <- 0.20; sigS <- 0.03; tauS <- 0.05
k_sigma <- 8; k_tau <- 16
halfw <- k_sigma * sigS + k_tau * tauS   # 1.04
a     <- muS - halfw                     # -0.84  (window lower edge)
b_win <- muS + halfw                     # 1.24   (window upper edge)

test_that("ss_exg_stop_success honors a finite upper deadline", {
  p <- make_texg_pars(muG, sigG, tauG, lbG = 0, muS, sigS, tauS, lbS = 0)

  # Auto (upper <= 0 -> Inf): integrate the full heuristic window [a, b_win].
  v_full <- EMC2:::ss_exg_stop_success_value(SSD, p, upper = -1,
                                             k_sigma = k_sigma, k_tau = k_tau)
  # Deadline inside the window: strictly less stop-success mass.
  u_mid  <- 0.5
  v_cap  <- EMC2:::ss_exg_stop_success_value(SSD, p, upper = u_mid,
                                             k_sigma = k_sigma, k_tau = k_tau)
  # Deadline beyond the window: capped back to b_win -> no-op.
  v_wide <- EMC2:::ss_exg_stop_success_value(SSD, p, upper = b_win + 5,
                                             k_sigma = k_sigma, k_tau = k_tau)

  expect_true(is.finite(v_full) && v_full > 0)
  # The canary: pre-fix `upper` was ignored so v_cap == v_full.
  expect_lt(v_cap, v_full)
  expect_equal(v_wide, v_full, tolerance = 1e-8)
})

test_that("plain-EXG matches TEXG when the stop truncation is negligible", {
  # Truncation always starts at 0 (an RT < 0 is impossible), so lbS=lbG=0 is the
  # physical "no extra truncation" case: with muS=0.20, sigS=0.03 the EXG mass
  # below 0 is ~6.7 sigma into the tail (negligible), so 1-F(0) ~ 1 and dtexg
  # reduces to dexg (ptexg to pexg). Both then integrate essentially the same
  # integrand up to `upper` (the plain-EXG window's sub-zero portion contributes
  # negligibly), so the two stop-success values must agree.
  p <- make_texg_pars(muG, sigG, tauG, lbG = 0, muS, sigS, tauS, lbS = 0)
  upper <- 0.5

  v_plain <- EMC2:::ss_exg_stop_success_value(SSD, p, upper = upper,
                                              k_sigma = k_sigma, k_tau = k_tau)
  v_texg  <- EMC2:::ss_texg_stop_success_value(SSD, p, method = "integrate",
                                               upper = upper,
                                               k_sigma = k_sigma, k_tau = k_tau)
  expect_equal(v_plain, v_texg, tolerance = 1e-4)
})
