# Extracted from test-stop_success_gl.R:335

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
pv_for <- function(des) {
    pv <- sampled_pars(des, doMap = FALSE)
    pv["mu_lMFALSE"] <- log(0.65); pv["mu_lMTRUE"] <- log(0.50)
    pv["sigma"] <- log(0.05); pv["tau"] <- log(0.15)
    pv["muS"] <- log(0.20); pv["sigmaS"] <- log(0.04); pv["tauS"] <- log(0.08)
    pv["gf"] <- qnorm(0.05); pv["tf"] <- qnorm(0.05)
    pv
  }
ll_for <- function(stop_method, stop_n_nodes = 64L) {
    x <- make_dadm(stop_method, stop_n_nodes)
    pv <- pv_for(x$des)
    proposals <- matrix(pv, nrow = 2, ncol = length(pv), byrow = TRUE)
    colnames(proposals) <- names(pv)
    model_fun <- attr(x$dadm, "model")
    lls <- calc_ll_manager(proposals, x$dadm, model_fun)
    # the model list carries the user's choice
    expect_identical(model_fun()$stop_method, stop_method)
    lls[1]
  }
ll_int  <- ll_for("integrate")
expect_identical(emc2_get_stop_method()$method, "integrate")
ll_gl   <- ll_for("gl")
expect_identical(emc2_get_stop_method()$method, "gl")
ll_auto <- ll_for("auto")
expect_identical(emc2_get_stop_method()$method, "auto")
ll_gl4  <- ll_for("gl", stop_n_nodes = 4L)
expect_identical(emc2_get_stop_method()$n_nodes, 4L)
expect_equal(ll_gl, ll_int, tolerance = 1e-5)
expect_equal(ll_auto, ll_int, tolerance = 1e-5)
expect_gt(abs(ll_gl4 - ll_int), abs(ll_gl - ll_int))
x <- make_dadm("gl", 128L)
pv <- pv_for(x$des)
r_gl  <- calc_ll_R(pv, attr(x$dadm, "model")(), x$dadm)
x2 <- make_dadm("integrate")
r_int <- calc_ll_R(pv, attr(x2$dadm, "model")(), x2$dadm)
expect_equal(r_gl, r_int, tolerance = 1e-5)
