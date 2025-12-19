## Stop-signal likelihood demo for EMC2: censoring/truncation and R vs C comparisons
# - Simulates SSexG / SShybrid data under various censoring/truncation settings
# - Compares R likelihood (old vs censored) against C likelihood

## Couple of notes: 
# 1. The default R definitions have NOT been updated. This function gives helpers to use the new likelihoods.
# Once happy, simply (a) update the default R likelihoods in model_SS.R to use the new functions,
# and (b) add the c_name exactly as written in this file.
# 2. The likelihoods do not currently support "false-alarm" stops (stops on go trials that aren't purely censored)
# This was not in the original EMC2 SS models so I did not build it in.

library(EMC2)
source("Censoring Test/test_likelihood_plotfuns.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- FALSE
RUN_PROFILE_PLOTS <- TRUE
# Helper functions
Rtfun <- function(d){ d$rt[d$rt==Inf]=NA; d$rt}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
make_ss_model <- function(model_type = c("SSexG", "SShybrid"),
                          use_censored_ll = TRUE,
                          normalise_trunc = FALSE) {
  model_type <- match.arg(model_type)
  function() {
    m <- if (model_type == "SSexG") SSexG() else SShybrid()
    m$c_name <- model_type

    if (use_censored_ll) {
      m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
        EMC2:::log_likelihood_race_ss_cens_trunc(
          pars, dadm, model, min_ll = min_ll, normalise_trunc = normalise_trunc
        )
      }
    } else {
      m$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
        EMC2:::log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll)
      }
    }
    m
  }
}

set_ss_pars <- function(p_vector, model_type, include_stop_triggered) {
  set_if_present <- function(name, value) {
    if (name %in% names(p_vector)) p_vector[[name]] <<- value
  }

  if (model_type == "SSexG") {
    for (nm in names(p_vector)) {
      if (grepl("^mu_", nm)) p_vector[[nm]] <- log(.55)
    }
    if (include_stop_triggered && "mu_lRst" %in% names(p_vector)) {
      p_vector[["mu_lRst"]] <- log(.75)
    }
    set_if_present("muS", log(.60))
    set_if_present("sigma", log(.18))
    set_if_present("tau", log(.20))
    set_if_present("sigmaS", log(.10))
    set_if_present("tauS", log(.12))
    set_if_present("gf", qnorm(.10))
    set_if_present("tf", qnorm(.10))
  } else {
    set_if_present("v", log(1.0))
    set_if_present("v_lMTRUE", log(1.2))
    set_if_present("v_lMFALSE", log(.4))
    set_if_present("B", log(1.2))
    set_if_present("A", log(.3))
    set_if_present("t0", log(.20))
    set_if_present("s", log(1.0))
    set_if_present("muS", log(.50))
    set_if_present("sigmaS", log(.12))
    set_if_present("tauS", log(.10))
    set_if_present("gf", qnorm(.10))
    set_if_present("tf", qnorm(.10))
  }

  p_vector
}

summarise_ss_data <- function(d, include_stop_triggered, Rlevels) {
  is_stop <- is.finite(d$SSD)
  is_nr <- is.na(d$R) | is.infinite(d$rt)
  resp <- as.character(d$R)
  resp[is.na(resp)] <- "NR"

  cat("\ninclude_stop_triggered =", include_stop_triggered, "\n")
  cat("Rlevels =", paste(Rlevels, collapse = ", "), "\n\n")

  cat("Counts by trial type (finite SSD => stop):\n")
  print(table(trial_type = ifelse(is_stop, "stop", "go")))

  cat("\nCounts by response (NR = no response):\n")
  print(table(response = resp))

  cat("\nCounts by trial type x response:\n")
  print(table(trial_type = ifelse(is_stop, "stop", "go"), response = resp))

  cat("\nNon-response rate by trial type:\n")
  print(tapply(is_nr, ifelse(is_stop, "stop", "go"), mean))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

}

run_ss_demo <- function(model_type = c("SSexG", "SShybrid"),
                        use_censored_ll = TRUE,
                        include_stop_triggered = FALSE,
                        n_trials = 500,
                        UC = Inf,
                        LT = 0,
                        UT = Inf,
                        normalise_trunc = FALSE,
                        make_deadline_censor = FALSE,
                        ssd_values = c(.15, .25),
                        p_ssd = c(.25, .25),
                        label = NULL) {
  model_type <- match.arg(model_type)
  if (normalise_trunc && LT == 0 && is.infinite(UT)) {
    stop("normalise_trunc=TRUE requires LT > 0 or finite UT.")
  }
  if (make_deadline_censor && !is.finite(UC)) {
    stop("make_deadline_censor=TRUE requires finite UC.")
  }

  st_label <- "st"
  Rlevels <- if (include_stop_triggered) c("left", "right", st_label) else c("left", "right")

  matchfun <- function(d) as.character(d$S) == as.character(d$lR)
  lIfun <- function(d) {
    if (!include_stop_triggered) {
      return(factor(rep(2, nrow(d)), levels = 1:2))
    }
    factor(ifelse(as.character(d$lR) == st_label, 1, 2), levels = 1:2)
  }

  mySSD_function <- function(d) SSD_function(d, SSD = ssd_values, pSSD = p_ssd)

  designSS <- design(
    model = make_ss_model(model_type, use_censored_ll = use_censored_ll,
                          normalise_trunc = normalise_trunc),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = Rlevels,
    matchfun = matchfun,
    functions = list(lI = lIfun, SSD = mySSD_function),
    formula = if (model_type == "SSexG") {
      list(mu ~ 0 + lR, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    } else {
      list(v ~ 0 + lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    }
  )

  p_vector <- sampled_pars(designSS, doMap = FALSE)
  p_vector <- set_ss_pars(p_vector, model_type, include_stop_triggered)

  dat <- make_data(
    p_vector, designSS, n_trials = n_trials,
    UC = UC, UCresponse = FALSE,
    LT = LT, UT = UT, LC = LT,
    return_Ffunctions = TRUE
  )

  if (make_deadline_censor) {
    dat$rt[is.na(dat$rt)] <- Inf
  }

  cat("\n--- SS demo ---\n")
  if (!is.null(label)) cat(label, "\n")
  cat("model_type =", model_type, "\n")
  cat("use_censored_ll =", use_censored_ll, " normalise_trunc =", normalise_trunc, "\n")
  cat("UC =", UC, " LT =", LT, " UT =", UT, "\n\n")
  summarise_ss_data(dat, include_stop_triggered, Rlevels)

  if (isTRUE(RUN_PROFILE_PLOTS)) {
    profile_plot_test(dat, designSS, p_vector, n_cores = 1, layout = c(3, 3),
                      use_c = FALSE, figure_title = paste("R Likelihood:", label))
    profile_plot_test(dat, designSS, p_vector, n_cores = 1, layout = c(3, 3),
                      use_c = TRUE, figure_title = paste("C Likelihood:", label))
  }

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designSS, type = "single")
    emc <- fit(emc, cores_for_chains = 3, fileName = "samples_ss_cens.RData")
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, n_post = 20, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designSS, true_pars = p_vector))
}

## First two tests demonstrate equivalence of new censored C likelihood vs existing R likelihood when no censoring present
# Then equivalence of the new R function with the C function
# NB I believe the tiny visible movements between R and C are just differences between the integration routines
# --- Test 1: SSexG baseline (no censoring/truncation) ------------------------
res_exg_baseline_old <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = FALSE,
  include_stop_triggered = FALSE,
  UC = Inf,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = FALSE,
  label = "SSexG baseline (old R likelihood)"
)
res_exg_baseline_cens <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = TRUE,
  include_stop_triggered = FALSE,
  UC = Inf,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = FALSE,
  label = "SSexG baseline (censored R likelihood)"
)

## Then show new C++ code does much better when censoring is introduced
# and that new R code matches C++ code
# NB I set the censor very low because I wasn't sure what parameters made sense to generate longer RTs
# --- Test 2: SSexG with deadline censoring -----------------------------------
res_exg_cens_old <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = FALSE,
  include_stop_triggered = FALSE,
  UC = 0.8,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = TRUE,
  label = "SSexG censoring only (old R likelihood)"
)
res_exg_cens <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = TRUE,
  include_stop_triggered = FALSE,
  UC = 1.5,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = TRUE,
  label = "SSexG censoring only (censored R likelihood)"
)

## Again test equivalence of new censored C likelihood vs existing R likelihood when no censoring present
# using the stop_triggered accumulator this time
# Then equivalence of the new R function with the C function
# --- Test 3: SSexG with stop-triggered response accumulator ------------------
res_exg_st_old <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = FALSE,
  include_stop_triggered = TRUE,
  UC = Inf,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = FALSE,
  label = "SSexG + ST accumulator (old R likelihood)"
)
res_exg_st_cens <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = TRUE,
  include_stop_triggered = TRUE,
  UC = .8,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = TRUE,
  label = "SSexG + ST accumulator (censored R likelihood)"
)

# --- Test 4: SSexG with truncation (LT/UT) -----------------------------------
res_exg_trunc_old <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = FALSE,
  include_stop_triggered = FALSE,
  UC = .8,
  LT = 0.15,
  UT = 2.0,
  normalise_trunc = TRUE,
  make_deadline_censor = FALSE,
  label = "SSexG truncation (old R likelihood)"
)
res_exg_trunc_cens <- run_ss_demo(
  model_type = "SSexG",
  use_censored_ll = TRUE,
  include_stop_triggered = FALSE,
  UC = .8,
  LT = 0.15,
  UT = 2.0,
  normalise_trunc = TRUE,
  make_deadline_censor = FALSE,
  label = "SSexG truncation (censored R likelihood)"
)

## Couple of tests with SShybrid model. Depending on UC value drift seems poorly identified but probably just my choice of parameters
# --- Test 5: SShybrid with deadline censoring --------------------------------
res_hyb_cens_old <- run_ss_demo(
  model_type = "SShybrid",
  use_censored_ll = FALSE,
  include_stop_triggered = FALSE,
  UC = 1.0,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = TRUE,
  label = "SShybrid censoring only (old R likelihood)"
)
res_hyb_cens <- run_ss_demo(
  model_type = "SShybrid",
  use_censored_ll = TRUE,
  include_stop_triggered = FALSE,
  UC = 1.0,
  LT = 0,
  UT = Inf,
  normalise_trunc = FALSE,
  make_deadline_censor = TRUE,
  label = "SShybrid censoring only (censored R likelihood)"
)

# --- Test 6: SShybrid with truncation (LT/UT) --------------------------------
res_hyb_trunc_cens <- run_ss_demo(
  model_type = "SShybrid",
  use_censored_ll = TRUE,
  include_stop_triggered = FALSE,
  UC = .9,
  LT = 0.15,
  UT = 2.0,
  normalise_trunc = TRUE,
  make_deadline_censor = FALSE,
  label = "SShybrid truncation (censored R likelihood)"
)
