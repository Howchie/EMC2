rm(list=ls())
library(EMC2)

check_dat <- function(dat) {
  datok <- dat[!is.na(dat$R),]
  if (length(levels(dat$R))>1) print(c(Accuracy=mean(datok$S==datok$R)))
  print(c(MRT=mean(datok$rt)))
  pStopSignal=mean(is.finite(dat$SSD))
  print(c(pStopSignal=pStopSignal,
          pStop=mean(is.finite(dat$SSD)&is.na(dat$R))/pStopSignal))
  print(c(SSD=sort(unique(dat$SSD))))
}


# NB: all datasets below are (re)simulated WITH finite UC censoring, so the old
# uncensored EXGintegrate.RData is intentionally NOT loaded here.
# load("WorkingTests/EXGintegrate.RData")

# "st" level is stop-triggered (lI=1), the rest are go (lI=2)
lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
p_ssd     = c(0.25)


# ===========================================================================
# DIRECT CHECK (fast, no fits) -- run after devtools::load_all(".")
# Confirms, for UC-censored stop-signal data, that:
#   (1) a finite upper limit (UC - SSD) routes through the GL fallback
#       (branch "gl_finite_upper"), never the analytic closed form, and the
#       result agrees with the reference route; with upper = Inf the n_go = 1
#       closed form is used (analytic_fullline / analytic_trunc);
#   (2) the LIVE likelihood (calc_ll) agrees across stop_method =
#       integrate / gl / auto / analytic on a UC-censored staircase dataset;
#   (3) the staircase treats an rt beyond UC as a successful stop, so a
#       censored stop trial steps the ladder UP (both simulator paths).
# This block is self-contained; it does not depend on the fits below.
#
# Reference notes:
#   * The C++ qags ("integrate") route cannot start from lbS = -Inf, so for the
#     UNTRUNCATED case (exg_lb = -Inf) the R integrate route (stats::integrate)
#     is the oracle. The TRUNCATED case (finite exg_lb) uses C++ integrate.
#   * Under a finite UC every no-response stop trial is a finite-upper
#     evaluation, so "auto"/"analytic" both fall back to GL (analytic only
#     fires for the upper = Inf, n_go = 1 case).
# ===========================================================================
cat("\n=========== DIRECT CHECK: finite UC -> GL fallback (EXG) ===========\n")
dc_pass <- TRUE
dc_report <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) dc_pass <<- FALSE
}

## summed C++ log-likelihood at one parameter vector, honouring stop_method
## exactly as calc_ll_manager() does (process-global config).
dc_summed_ll <- function(dat, design, p_vector, model_fn, method) {
  dadm    <- EMC2:::design_model(dat, design, verbose = FALSE)
  model_d <- attr(dadm, "model")()
  EMC2:::set_stop_method_from_model(model_fn(stop_method = method))
  on.exit(EMC2:::set_stop_method_from_model(model_fn()), add = TRUE)
  p_matrix <- matrix(p_vector, nrow = 1); colnames(p_matrix) <- names(p_vector)
  designs <- list()
  for (p in names(model_d$p_types))
    designs[[p]] <- attr(dadm, "designs")[[p]][
      attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  sum(EMC2:::calc_ll(p_matrix, dadm, constants, designs,
                     model_d$c_name, model_d$bound, model_d$transform,
                     model_d$pre_transform, names(model_d$p_types),
                     log(1e-10), model_d$trend))
}

## ---- (1) integral level: branch dispatch + agreement --------------------
local({
  # pars matrix layout expected by ss_texg_stop_success_* (see make_pars in
  # WorkingTests/stop_success_methods.R): col 1:3 = muG/sigG/tauG per go racer,
  # row1 cols 4:6 = muS/sigS/tauS, col 9 = lbG, row1 col 10 = lbS.
  mk_pars <- function(muS, sigS, tauS, lbS, muG, sigG, tauG, lbG) {
    n <- length(muG); m <- matrix(0, n, 10)
    m[, 1] <- muG; m[, 2] <- sigG; m[, 3] <- tauG
    m[1, 4] <- muS; m[1, 5] <- sigS; m[1, 6] <- tauS
    m[, 9] <- lbG; m[1, 10] <- lbS; m
  }
  SSD <- .25; uc_eff <- .55 - SSD                  # finite upper = UC - SSD

  ## untruncated (lbS = lbG = -Inf): closed form is the full-line 4-pnorm;
  ## C++ qags cannot start from -Inf, so the R integrate route is the oracle.
  Rint <- function(up) EMC2:::stop_success_texg_R(
    c(.2, .5), c(.05, .1), c(.1, .25), c(-Inf, -Inf), SSD, upper = up,
    method = "integrate")
  p_un <- mk_pars(.2, .05, .1, -Inf, .5, .1, .25, -Inf)
  b_un_inf <- ss_texg_stop_success_auto_branch(SSD, p_un)
  b_un_fin <- ss_texg_stop_success_auto_branch(SSD, p_un, uc_eff)
  a_un_inf <- ss_texg_stop_success_value(SSD, p_un, "auto")
  a_un_fin <- ss_texg_stop_success_value(SSD, p_un, "auto", upper = uc_eff)
  o_un_inf <- Rint(Inf); o_un_fin <- Rint(uc_eff)
  r_un_inf <- abs(a_un_inf - o_un_inf) / o_un_inf
  r_un_fin <- abs(a_un_fin - o_un_fin) / o_un_fin
  dc_report(sprintf("untrunc upper=Inf branch = %s", b_un_inf),
            b_un_inf == "analytic_fullline")
  dc_report(sprintf("untrunc finite upper branch = %s", b_un_fin),
            b_un_fin == "gl_finite_upper")
  dc_report(sprintf("untrunc auto vs oracle, upper=Inf    (rel %.1e)", r_un_inf),
            r_un_inf < 1e-4)
  dc_report(sprintf("untrunc auto vs oracle, finite upper (rel %.1e)", r_un_fin),
            r_un_fin < 1e-4)
  dc_report("untrunc finite upper reduces P(stop) vs Inf", a_un_fin < a_un_inf)

  ## truncated (finite lb): C++ qags is valid -> compare auto vs integrate.
  p_tr <- mk_pars(.2, .05, .1, .2, .5, .1, .25, .2)
  b_tr_inf <- ss_texg_stop_success_auto_branch(SSD, p_tr)
  b_tr_fin <- ss_texg_stop_success_auto_branch(SSD, p_tr, uc_eff)
  ta_inf <- ss_texg_stop_success_value(SSD, p_tr, "auto")
  ti_inf <- ss_texg_stop_success_value(SSD, p_tr, "integrate")
  ta_fin <- ss_texg_stop_success_value(SSD, p_tr, "auto", upper = uc_eff)
  ti_fin <- ss_texg_stop_success_value(SSD, p_tr, "integrate", upper = uc_eff)
  r_tr_inf <- abs(ta_inf - ti_inf) / ti_inf
  r_tr_fin <- abs(ta_fin - ti_fin) / ti_fin
  dc_report(sprintf("trunc upper=Inf branch = %s", b_tr_inf),
            b_tr_inf == "analytic_trunc")
  dc_report(sprintf("trunc finite upper branch = %s", b_tr_fin),
            b_tr_fin == "gl_finite_upper")
  dc_report(sprintf("trunc auto vs integrate, upper=Inf    (rel %.1e)", r_tr_inf),
            r_tr_inf < 1e-4)
  dc_report(sprintf("trunc auto vs integrate, finite upper (rel %.1e)", r_tr_fin),
            r_tr_fin < 1e-4)
})

## ---- (2) end-to-end calc_ll on a UC-censored staircase dataset -----------
## Truncated (finite exg_lb) so the C++ integrate route is valid; all four
## stop_methods must agree.
dc_design <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1), Rlevels = 1,
  TC = list(UC = 0.7),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st, SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                  gf ~ 1, tf ~ 1, exg_lb ~ 1, exgS_lb ~ 1),
  constants = c(exg_lb = .2, exgS_lb = .2))
dc_pv <- sampled_pars(dc_design, doMap = FALSE)
dc_pv[] <- c(log(.5), log(.1), log(.25), log(.2), log(.05), log(.1), qnorm(.08), qnorm(.08))
set.seed(2024)
dc_dat <- make_data(dc_pv, dc_design, n_trials = 800,
                    staircase = list(SSD0 = .25, stairstep = .05, stairmin = 0, stairmax = Inf))
is_stop <- is.finite(dc_dat$SSD)
cat(sprintf("  UC=0.7 staircase data: %d trials, %.0f%% stop signal, %.0f%% stop-success\n",
            nrow(dc_dat), 100 * mean(is_stop),
            100 * mean(is.na(dc_dat$R[is_stop]))))
local({
  lls <- sapply(c("integrate", "gl", "auto", "analytic"),
                function(m) dc_summed_ll(dc_dat, dc_design, dc_pv, SSEXG, m))
  cat("  summed log-likelihood by stop_method:\n")
  print(round(lls, 6))
  d <- max(abs(lls - lls["integrate"]))
  dc_report(sprintf("calc_ll: gl/auto/analytic match integrate (max |diff| %.1e)", d),
            d < 1e-3)
})

## ---- (3) staircase x UC: rt > UC counts as a stop success (ladder UP) -----
## (a) unconditional simulator path (stair_advance: explicit rt > UC -> NA)
local({
  sc <- list(SSD0 = .3, stairstep = .05, stairmin = 0, stairmax = Inf, UC = .5)
  d0 <- data.frame(subjects = factor("1"),
                   lR = factor("go", levels = c("go", "st")),
                   lI = factor(2, levels = 1:2))
  # the simulators always read the current SSD (stair_ssd_value, which also
  # initialises the ladder) before stepping it (stair_advance); mirror that.
  mk <- function() { s <- EMC2:::stair_state_init(d0, sc)
                     EMC2:::stair_ssd_value(s, 1, "all"); s }
  # go response AFTER UC -> recoded stop success -> ladder UP
  s_up <- mk(); ssd0 <- EMC2:::stair_ssd_value(s_up, 1, "all")
  EMC2:::stair_advance(s_up, 1, "all", Rlab = "go", rt = 0.60, data_row = d0)
  up <- EMC2:::stair_ssd_value(s_up, 1, "all")
  # go response BEFORE UC -> genuine go response -> ladder DOWN
  s_dn <- mk(); EMC2:::stair_advance(s_dn, 1, "all", Rlab = "go", rt = 0.40, data_row = d0)
  dn <- EMC2:::stair_ssd_value(s_dn, 1, "all")
  dc_report(sprintf("stair_advance: rt>UC steps UP (%.2f->%.2f), rt<UC steps DOWN (->%.2f)",
                    ssd0, up, dn), up > ssd0 && dn < ssd0)
})
## (b) conditional/block simulator path: the rfun censors go finishing times
## beyond UC to Inf; staircase_function must then score that trial as a stop
## (Ri==1) and step the ladder UP.
local({
  sc <- list(SSD0 = .3, stairstep = .05, stairmin = 0, stairmax = Inf, labels = NULL)
  # rows: 1 = stop accumulator, 2 = go accumulator. Trial 1's go is censored
  # to Inf (as the rfun does when min(go finish) > UC); trial 2 responds.
  dts <- matrix(c(0.40, Inf,     # trial 1: go censored -> stop success -> UP
                  0.40, 0.20),   # trial 2: go responds before stop -> DOWN
                nrow = 2)
  res <- EMC2:::staircase_function(dts, sc)
  dc_report(sprintf("staircase_function: censored go -> stop (sR NA), SSD %.2f->%.2f UP",
                    res$SSD[1], res$SSD[2]),
            is.na(res$sR[1]) && res$SSD[2] > res$SSD[1])
})

cat(sprintf("=========== DIRECT CHECK %s ===========\n\n",
            if (dc_pass) "PASSED" else "FAILED"))


##### n = 1 ----
# For n = 1 used the following, n>1 case uses the defaults
rt_resolution <- stairstep <- .2
UCfun <- function(d) quantile(d$rt, probs = .9, na.rm = TRUE)

#### n = 1, untruncated ----
# NB: the C++ "integrate"/qags route cannot start from lbS = -Inf, so for the
# untruncated model the valid routes are "auto" (exact analytic at upper = Inf,
# GL fallback under UC) and "gl". The integrate fit is included for reference
# but is expected to be biased on stop-success trials.

designSSEXG1f <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
p_vector <- sampled_pars(designSSEXG1f, doMap = FALSE)
p_vector[] <- c(log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=400
dat1 <- make_data(p_vector, designSSEXG1f, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=stairstep,stairmin=0,stairmax=Inf))
check_dat(dat1)

emcR1f <- make_emc(dat1,designSSEXG1f,type="single",rt_resolution = rt_resolution)
tAf <- system.time({emcR1f <- fit(emcR1f)})[3]
recovery(emcR1f,p_vector,main="n=1, auto untruncated")

designSSEXG1gl <- design(model = SSEXG("gl"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
emcR1gl <- make_emc(dat1,designSSEXG1gl,type="single",rt_resolution = rt_resolution)
tAfGL <- system.time({emcR1gl <- fit(emcR1gl)})[3]
recovery(emcR1gl,p_vector,main="n=1, GL untruncated")

designSSEXG1i <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
emcR1i <- make_emc(dat1,designSSEXG1i,type="single",rt_resolution = rt_resolution)
tAfi <- system.time({emcR1i <- fit(emcR1i)})[3]
recovery(emcR1i,p_vector,main="n=1, integrate untruncated (NB: qags cannot handle lb=-Inf)")


#### n = 1, truncated ----

designSSEXG1t <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
p_vector <- sampled_pars(designSSEXG1t, doMap = FALSE)
p_vector[] <- c(log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=400
dat1t <- make_data(p_vector, designSSEXG1t, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=stairstep,stairmin=0,stairmax=Inf))
check_dat(dat1t)

emcR1t <- make_emc(dat1t,designSSEXG1t,type="single",rt_resolution = rt_resolution)
tAt <- system.time({emcR1t <- fit(emcR1t)})[3]
recovery(emcR1t,p_vector,main="n=1, auto truncated")

designSSEXG1tgl <- design(model = SSEXG("gl"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
emcR1tGL <- make_emc(dat1t,designSSEXG1tgl,type="single",rt_resolution = rt_resolution)
tAtGL <- system.time({emcR1tGL <- fit(emcR1tGL)})[3]
recovery(emcR1tGL,p_vector,main="n=1, GL truncated")

designSSEXG1ti <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
emcR1ti <- make_emc(dat1t,designSSEXG1ti,type="single",rt_resolution = rt_resolution)
tAti <- system.time({emcR1ti <- fit(emcR1ti)})[3]
recovery(emcR1ti,p_vector,main="n=1, integrate truncated")


#### 2 choice EXG ----
designSSEXG2 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSEXG2, doMap = FALSE)
p_vector[] <- c(log(.69),log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=200
dat2 <- make_data(p_vector, designSSEXG2, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat2)

emcR2 <- make_emc(dat2,designSSEXG2,type="single")
t2 <- system.time({emcR2 <- fit(emcR2)})[3]
recovery(emcR2,p_vector,main="n=2, auto")

designSSEXG2i <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR2i <- make_emc(dat2,designSSEXG2i,type="single")
t2i <- system.time({emcR2i <- fit(emcR2i)})[3]
recovery(emcR2i,p_vector,main="n=2, integrate")

#### 4 choice EXG ----
designSSEXG4 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSEXG4, doMap = FALSE)
p_vector[] <- c(log(.8),log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=100
dat4 <- make_data(p_vector, designSSEXG4, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat4)

emcR4 <- make_emc(dat4,designSSEXG4,type="single")
t4=system.time({emcR4 <- fit(emcR4)})[3]
recovery(emcR4,p_vector,main="n=4, auto")

designSSEXG4i <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=UCfun,verbose=TRUE),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR4i <- make_emc(dat4,designSSEXG4i,type="single")
t4i=system.time({emcR4i <- fit(emcR4i)})[3]
recovery(emcR4i,p_vector,main="n=4, integrate")

# save(dat1,dat1t,dat2,dat4,file="WorkingTests/EXGintegrate_UC.RData")

#### print timing results
print(c(tAf=tAf,tAfGL=tAfGL,tAfi=tAfi))
print(c(tAt=tAt,tAtGL=tAtGL,tAti=tAti))
print(c(t2=t2,t2i=t2i))
print(c(t4=t4,t4i=t4i))

times <- c(tAf=tAf,tAfGL=tAfGL,tAfi=tAfi,tAt=tAt,tAtGL=tAtGL,tAti=tAti,
       t2=t2,t2i=t2i,t4=t4,t4i=t4i)
save(times,file="tEXG.RData")
