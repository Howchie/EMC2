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


# NB: all datasets below are (re)simulated WITH the finite UC censoring, so the
# old uncensored RDEXintegrate.RData is intentionally NOT loaded here.
# load("WorkingTests/RDEXintegrate.RData")

# "st" level is stop-triggered (lI=1), the rest are go (lI=2)
lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
p_ssd     = c(0.25)


# ===========================================================================
# DIRECT CHECK (fast, no fits) -- run after devtools::load_all(".")
# Confirms, for UC-censored stop-signal data, that:
#   (1) a finite upper limit (UC - SSD) routes through the GL fallback and
#       agrees with the untouched integrate/qags route (RDEX has no analytic
#       form, so "auto" == GL everywhere; the test is finite vs Inf upper);
#   (2) the LIVE likelihood (calc_ll) agrees across stop_method =
#       integrate / gl / auto on a UC-censored staircase dataset;
#   (3) the staircase treats an rt beyond UC as a successful stop, so a
#       censored stop trial steps the ladder UP (both simulator paths).
# This block is self-contained; it does not depend on the fits below.
# ===========================================================================
cat("\n=========== DIRECT CHECK: finite UC -> GL fallback (RDEX) ===========\n")
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

## ---- (1) integral level: finite upper (= UC - SSD) agrees with integrate --
local({
  muS <- .2; sigS <- .05; tauS <- .1; lbS <- 0     # stop ex-Gaussian
  v <- 1.5; B <- 1; A <- 0; t0 <- .2; s <- 1        # single go (RDM) racer
  SSD <- .25; uc_eff <- .55 - SSD                   # finite upper = UC - SSD
  pin <- function(up, m) EMC2:::stop_success_rdex_R(
    n_acc = 1, mu = muS, sigma = sigS, tau = tauS, lb = lbS,
    v = v, B = B, A = A, t0 = t0, s = s, SSD = SSD, upper = up, method = m)
  inf_i <- pin(Inf, "integrate");  inf_a <- pin(Inf, "auto")
  fin_i <- pin(uc_eff, "integrate"); fin_a <- pin(uc_eff, "auto")
  r_inf <- abs(inf_a - inf_i) / inf_i
  r_fin <- abs(fin_a - fin_i) / fin_i
  dc_report(sprintf("auto==integrate, upper=Inf    (rel %.1e)", r_inf), r_inf < 1e-4)
  dc_report(sprintf("auto==integrate, finite upper (rel %.1e)", r_fin), r_fin < 1e-4)
  dc_report("finite upper reduces P(stop) vs Inf", fin_a < inf_a)
})

## ---- (2) end-to-end calc_ll on a UC-censored staircase dataset -----------
dc_design <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1), Rlevels = 1,
  TC = list(UC = 0.6),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st, SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 1, B ~ 1, t0 ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
dc_pv <- sampled_pars(dc_design, doMap = FALSE)
dc_pv[] <- c(log(1.5), log(1), log(0.2), log(0.2), log(0.05), log(0.1), qnorm(0.08), qnorm(0.08))
set.seed(2024)
dc_dat <- make_data(dc_pv, dc_design, n_trials = 800,
                    staircase = list(SSD0 = .25, stairstep = .05, stairmin = 0, stairmax = Inf))
is_stop <- is.finite(dc_dat$SSD)
cat(sprintf("  UC=0.6 staircase data: %d trials, %.0f%% stop signal, %.0f%% stop-success\n",
            nrow(dc_dat), 100 * mean(is_stop),
            100 * mean(is.na(dc_dat$R[is_stop]))))
local({
  lls <- sapply(c("integrate", "gl", "auto"),
                function(m) dc_summed_ll(dc_dat, dc_design, dc_pv, SSRDEX, m))
  cat("  summed log-likelihood by stop_method:\n")
  print(round(lls, 6))
  d <- max(abs(lls - lls["integrate"]))
  dc_report(sprintf("calc_ll: gl/auto match integrate (max |diff| %.1e)", d), d < 1e-3)
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


#### 1 choice RDEX ----
designSSRDEX1 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 1, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX1, doMap = FALSE)
p_vector[] <- c(log(1.5),log(1),log(0.2),log(0.2),
                log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=400

dat1 <- make_data(p_vector, designSSRDEX1, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat1)


emcR1 <- make_emc(dat1,designSSRDEX1,type="single")
t1=system.time({emcR1 <- fit(emcR1)})[3]
recovery(emcR1,p_vector,main="n=1, auto")

designSSRDEX1i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 1, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR1i <- make_emc(dat1,designSSRDEX1i,type="single")
t1i=system.time({emcR1i <- fit(emcR1i)})[3]
recovery(emcR1i,p_vector,main="n=1, integrate")

#### 2 choice RDEX ----
designSSRDEX2 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX2, doMap = FALSE)
p_vector[] <- c(log(0.25),log(1.5),log(1),log(0.2),
                log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=200
dat2 <- make_data(p_vector, designSSRDEX2, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat2)

emcR2 <- make_emc(dat2,designSSRDEX2,type="single")
t2 <- system.time({emcR2 <- fit(emcR2)})[3]
recovery(emcR2,p_vector,main="n=2, auto")

designSSRDEX2i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR2i <- make_emc(dat2,designSSRDEX2i,type="single")
t2i <- system.time({emcR2i <- fit(emcR2i)})[3]
recovery(emcR2i,p_vector,main="n=2, integrate")

#### 4 choice RDEX ----
designSSRDEX4 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX4, doMap = FALSE)
p_vector[] <- c(log(0.1),log(3),log(2),log(0.2),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=100
dat4 <- make_data(p_vector, designSSRDEX4, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat4)

emcR4 <- make_emc(dat4,designSSRDEX4,type="single")
t4=system.time({emcR4 <- fit(emcR4)})[3]
# Time difference of 3.419095 mins
recovery(emcR4,p_vector,main="n=4, auto")

designSSRDEX4i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR4i <- make_emc(dat4,designSSRDEX4i,type="single")
t4i=system.time({emcR4i <- fit(emcR4i)})[3]
recovery(emcR4i,p_vector,main="n=4, integrate")

# save(dat1,dat2,dat4,file="WorkingTests/RDEXintegrate.RData")

#### print timing results
print(c(t1=t1,t1i=t1i))
print(c(t2=t2,t2i=t2i))
print(c(t4=t4,t4i=t4i))

times=c(t1=t1,t1i=t1i,t2=t2,t2i=t2i,t4=t4,t4i=t4i)

save(times,file="tRDEX.RData")

