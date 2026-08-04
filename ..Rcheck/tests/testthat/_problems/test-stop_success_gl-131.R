# Extracted from test-stop_success_gl.R:131

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
grid <- expand.grid(sigma = c(.05, .03, .02, .01, .005, .002),
                      tau   = c(.04, .08, .12),
                      SSD   = c(0, .15, .3))
for (i in seq_len(nrow(grid))) {
    p <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0.05,
                        muS = 0.20, sigS = grid$sigma[i], tauS = grid$tau[i],
                        lbS = 0.05)
    ref  <- ss_texg_stop_success_value(grid$SSD[i], p, method = "integrate")
    auto <- ss_texg_stop_success_value(grid$SSD[i], p, method = "auto")
    expect_equal(auto, ref, tolerance = 1e-4,
                 info = sprintf("sigma=%.3f tau=%.3f SSD=%.2f",
                                grid$sigma[i], grid$tau[i], grid$SSD[i]))
  }
