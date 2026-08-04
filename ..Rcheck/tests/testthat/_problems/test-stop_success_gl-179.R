# Extracted from test-stop_success_gl.R:179

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
for (SSD in c(0, .15, .4)) {
    p <- make_texg_pars(0.50, 0.05, 0.08, 0.05, 0.20, 0.03, 0.05, 0.05)
    cpp <- ss_texg_stop_success_value(SSD, p, method = "analytic")
    r   <- stop_success_texg_analytic1_R(mu = c(.2, .5), sigma = c(.03, .05),
                                         tau = c(.05, .08), lb = c(.05, .05),
                                         SSD = SSD)
    expect_equal(r, cpp, tolerance = 1e-12)
  }
