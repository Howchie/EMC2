# Extracted from test-stop_success_gl.R:119

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
on.exit(emc2_set_stop_method("auto", 64L), add = TRUE)
SSD <- 0.15
p <- make_texg_pars(0.50, 0.05, 0.08, 0, 0.20, 0.02, 0.04, 0)
emc2_set_stop_method("integrate")
expect_identical(emc2_get_stop_method()$method, "integrate")
expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "integrate"))
emc2_set_stop_method("gl", 128L)
expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "gl",
                                              n_nodes = 128L))
emc2_set_stop_method("gl", 4L)
gl4 <- ss_texg_stop_success_value(SSD, p, method = "live")
emc2_set_stop_method("gl", 128L)
gl128 <- ss_texg_stop_success_value(SSD, p, method = "live")
ref <- ss_texg_stop_success_value(SSD, p, method = "integrate")
expect_gt(abs(gl4 - ref), abs(gl128 - ref))
expect_equal(gl128, ref, tolerance = 1e-5)
emc2_set_stop_method("auto", 64L)
expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "auto"))
pr <- make_rdex_pars(v = 2.5, B = 0.6, A = 0.2, t0 = 0.15, s = 1.0,
                       muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0)
emc2_set_stop_method("integrate")
expect_identical(ss_rdex_stop_success_value(SSD, pr, method = "live"),
                   ss_rdex_stop_success_value(SSD, pr, method = "integrate"))
emc2_set_stop_method("gl", 96L)
expect_identical(ss_rdex_stop_success_value(SSD, pr, method = "live"),
                   ss_rdex_stop_success_value(SSD, pr, method = "gl",
                                              n_nodes = 96L))
