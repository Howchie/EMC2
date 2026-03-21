library(EMC2)
rm(list=ls())
source("WorkingTests/test_likelihood_plotfuns_ah.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

RUN_FITS         <- TRUE
RUN_PROFILE_PLOTS <- TRUE
#hist(EMC2:::rexG(1000,.5,.1,.35))
# ── helpers ───────────────────────────────────────────────────────────────────

# All accumulators are go
lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)

# "st" level is stop-triggered (lI=1), the rest are go (lI=2)
lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)

# Three-level mu factor so the ST accumulator gets its own parameter rather
# than sharing the mismatch mu.
AccTypeFun <- function(d) {
  is_st    <- as.character(d$lR) == "st"
  is_match <- !is_st & (as.character(d$S) == as.character(d$lR))
  factor(ifelse(is_st, "st", ifelse(is_match, "match", "mismatch")),
         levels = c("mismatch", "match", "st"))
}

summarise_ss_data <- function(dat, label) {
  is_stop <- is.finite(dat$SSD)
  resp    <- as.character(dat$R); resp[is.na(resp)] <- "NR"
  cat("\n───", label, "───\n")
  cat("Trial type counts:\n");  print(table(trial_type = ifelse(is_stop, "stop", "go")))
  cat("Response counts:\n");    print(table(response = resp))
  cat("Response × trial type:\n")
  print(table(trial_type = ifelse(is_stop, "stop", "go"), response = resp))
  cat("NR rate by trial type:\n")
  print(tapply(is.na(dat$R), ifelse(is_stop, "stop", "go"), mean))
  cat("Mean RT (finite only):\n")
  print(tapply(dat$rt, dat$S, function(x) mean(x[is.finite(x)], na.rm = TRUE)))
}

# core demo function
run_ss_demo <- function(
    include_stop_triggered = FALSE,
    n_trials  = 1000,
    UC        = Inf,
    ssd_vals  = c(0.20, 0.40),
    p_ssd     = c(0.25, 0.25),
    label     = NULL,
    cores_per_chain = 3,
    cores_for_chains = 3
) {

  if (!include_stop_triggered) {

    designSS <- design(
      model    = SSEXG,
      factors  = list(subjects = 1, S = c("left", "right")),
      Rlevels  = c("left", "right"),
      matchfun = function(d) as.character(d$S) == as.character(d$lR),
      functions = list(lI = lIfun_go, SSD = function(d) SSD_function(d, SSD = ssd_vals, pSSD = p_ssd)),
      formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                      muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    )

    p_vector <- sampled_pars(designSS, doMap = FALSE)
    p_vector["mu_lMFALSE"] <- log(0.65)   # incorrect accumulator
    p_vector["mu_lMTRUE"]  <- log(0.50)   # correct accumulator
    p_vector["sigma"]      <- log(0.15)
    p_vector["tau"]        <- log(0.5)
    p_vector["muS"]        <- log(0.2)
    p_vector["sigmaS"]     <- log(0.08)
    p_vector["tauS"]       <- log(0.1)
    p_vector["gf"]         <- qnorm(0.08)
    p_vector["tf"]         <- qnorm(0.08)

  } else {

    # The ST accumulator races from SSD time (like the stop process) and
    # produces an observable overt response.  Its mu is the finish time from
    # SSD onset, so effective RT ≈ SSD + mu_st + tau ≈ 0.28 + 0.38 + 0.15 ≈ 0.81 s.
    # The stop process still races and can inhibit go accumulators; when stop
    # wins, only the ST accumulator remains to produce a response.
    designSS <- design(
      model    = SSEXG,
      factors  = list(subjects = 1, S = c("left", "right")),
      Rlevels  = c("left", "right", "st"),
      matchfun = function(d) as.character(d$S) == as.character(d$lR),
      functions = list(lI = lIfun_st, SSD = function(d) SSD_function(d, SSD = ssd_vals, pSSD = p_ssd), AccType = AccTypeFun),
      formula  = list(mu ~ 0 + AccType, sigma ~ 1, tau ~ 1,
                      muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    )

    p_vector <- sampled_pars(designSS, doMap = FALSE)
    p_vector["mu_AccTypemismatch"] <- log(0.65)
    p_vector["mu_AccTypematch"]    <- log(0.50)
    p_vector["mu_AccTypest"]       <- log(0.38)  # finish time from SSD
    p_vector["sigma"]      <- log(0.15)
    p_vector["tau"]        <- log(0.5)
    p_vector["muS"]        <- log(0.20)
    p_vector["sigmaS"]     <- log(0.08)
    p_vector["tauS"]       <- log(0.12)
    p_vector["gf"]         <- qnorm(0.05)
    p_vector["tf"]         <- qnorm(0.05)
  }

  TC  <- list(UC = UC, LC = 0, LT = 0, UT = Inf, verbose = FALSE)
  dat <- make_data(p_vector, designSS, n_trials = n_trials, TC = TC)

  lbl <- if (!is.null(label)) label else
    paste0(if (include_stop_triggered) "ST" else "no-ST",
           ", UC=", UC)

  summarise_ss_data(dat, lbl)

  if (isTRUE(RUN_PROFILE_PLOTS)) {
    profile_plot_test(dat, designSS, p_vector,
                      n_cores = 1,
                      use_c   = TRUE,
                      figure_title = paste("C++:", lbl))
  }

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designSS, type = "single")
    emc <- fit(emc,fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(lbl)), ".RData"),stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain=cores_per_chain, cores_for_chains = cores_for_chains), max_tries=30)
    pp  <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = pp, n_post = 20, true_pars = p_vector)
    return(invisible(list(data = dat, design = designSS, true_pars = p_vector, emc=emc, pp=pp)))
  }

  invisible(list(data = dat, design = designSS, true_pars = p_vector))
}

# ── Test scenarios ────────────────────────────────────────────────────────────

# 1. Baseline: no censoring, standard stop-signal
res_no_st <- run_ss_demo(
  include_stop_triggered = FALSE,
  UC    = Inf,
  label = "SSexG no-ST UC=Inf"
)
if (RUN_FITS) print(recovery(res_no_st$emc, true_pars = res_no_st$true_pars))

# 2. Deadline censoring, no ST.
# Observed stop wins become part of the censored-NR mass instead of true NR.
# The C++ integrates over P(stop wins before UC) + P(go > UC) for deadline-NR rows.
res_no_st_cens <- run_ss_demo(
  include_stop_triggered = FALSE,
  UC    = 1.2,
  label = "SSexG no-ST UC=1.2"
)
if (RUN_FITS) print(recovery(res_no_st_cens$emc, true_pars = res_no_st_cens$true_pars))

# 3. Stop-triggered accumulator, no censoring.
# When stop wins, go accumulators are inhibited but the ST accumulator continues
# and produces an overt "st" response.  NR can only happen if BOTH go failure
# AND trigger failure occur simultaneously (≈ 0.25% of stop trials here).
res_st <- run_ss_demo(
  include_stop_triggered = TRUE,
  UC    = Inf,
  label = "SSexG ST UC=Inf"
)
if (RUN_FITS) print(recovery(res_st$emc, true_pars = res_st$true_pars))

# 4. Stop-triggered accumulator + deadline censoring.
# The corrected C++ censoring formula for the ST case is:
#   P(no resp by UC | triggered) = S_ST(UC−SSD) × [gf + (1−gf)×core_trig]
# where core_trig = P(stop wins before UC) + S_go(UC)×S_stop(UC−SSD).
res_st_cens <- run_ss_demo(
  include_stop_triggered = TRUE,
  UC    = 1.2,
  label = "SSexG ST UC=1.2"
)
if (RUN_FITS) print(recovery(res_st_cens$emc, true_pars = res_st_cens$true_pars))
save.image("demo_SSEXG.RData")