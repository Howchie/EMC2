# Extracted from test-stop_success_gl.R:143

# prequel ----------------------------------------------------------------------
make_texg_pars <- function(muG, sigG, tauG, lbG, muS, sigS, tauS, lbS) {
  # column layout expected by ss_texg_stop_success_*: muG=0 sigG=1 tauG=2
  # muS=3 sigS=4 tauS=5 (6,7 unused) lbG=8 lbS=9   (0-based -> R 1-based below)
  matrix(c(muG, sigG, tauG, muS, sigS, tauS, 0, 0, lbG, lbS),
         nrow = 1, ncol = 10)
}
make_rdex_pars <- function(v, B, A, t0, s, muS, sigS, tauS, lbS) {
  # go: v=0 B=1 A=2 t0=3 s=4 ; stop: muS=5 sigS=6 tauS=7 (8,9 unused) lbS=10
  matrix(c(v, B, A, t0, s, muS, sigS, tauS, 0, 0, lbS),
         nrow = 1, ncol = 11)
}

# test -------------------------------------------------------------------------
p   <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0.05,
                        muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0.05)
for (SSD in c(0, .1, .2, .3, .5)) {
    expect_identical(ss_texg_stop_success_auto_branch(SSD, p), "analytic_trunc")
    ana <- ss_texg_stop_success_value(SSD, p, method = "analytic")
    # tight adaptive reference (small tolerances, wide window)
    ref <- ss_texg_stop_success_value(SSD, p, method = "integrate",
                                      max_subdiv = 200, abs_tol = 1e-12,
                                      rel_tol = 1e-10, k_sigma = 10,
                                      k_tau = 20)
    expect_equal(ana, ref, tolerance = 1e-7, info = sprintf("SSD=%.2f", SSD))
  }
