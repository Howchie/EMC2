# Extracted from test-stop_success_gl.R:290

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
skip_if_not_installed("EMC2")
on.exit(emc2_set_stop_method("auto", 64L), add = TRUE)
make_dadm <- function(stop_method, stop_n_nodes = 64L) {
    des <- design(
      model    = SSEXG(stop_method = stop_method, stop_n_nodes = stop_n_nodes),
      factors  = list(subjects = 1, S = c("left", "right")),
      Rlevels  = c("left", "right"),
      matchfun = function(d) as.character(d$S) == as.character(d$lR),
      functions = list(
        lI  = function(d) factor(rep(2, nrow(d)), levels = 1:2),
        SSD = function(d) rep(Inf, nrow(d))),
      formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                      muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
      verbose  = FALSE
    )
    # two stop trials (one successful stop, one failed) + one go trial
    dat <- data.frame(
      subjects = factor(rep("1", 3)),
      S        = factor(c("left", "left", "right"),
                        levels = c("left", "right")),
      R        = factor(c(NA, "left", "right"),
                        levels = c("left", "right")),
      rt       = c(NA, 0.55, 0.61),
      SSD      = c(0.2, 0.25, Inf),
      trials   = 1:3
    )
    list(des = des, dadm = design_model(dat, des, verbose = FALSE))
  }
