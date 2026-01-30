library(EMC2)
library(dplyr)
source("Censoring Test/test_likelihood_plotfuns_ah.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- FALSE
RUN_PROFILE_PLOTS <- TRUE
# Helper functions
Rtfun <- function(d){ d$rt[d$rt==Inf]=NA; d$rt}
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)

set_ss_pars <- function(p_vector, model_type, include_stop_triggered) {
  set_if_present <- function(name, value) {
    for (nm in names(p_vector)) {
      if (grepl(name,nm)) p_vector[[nm]] <<- value
    }
  }

  if (model_type()$c_name == "SSEXG") {
    set_if_present("^mu_.*TRUE*", log(.45))
    set_if_present("^mu_.*FALSE", log(.6))
    set_if_present("^mu_RateST", log(.3))
    set_if_present("^muS", log(.4))
    set_if_present("^sigma", log(.4))
    set_if_present("^tau", log(.5))
    set_if_present("^sigmaS", log(.55))
    set_if_present("^tauS", log(.25))
    set_if_present("^gf", qnorm(.1))
    set_if_present("^tf", qnorm(.05))
  } else {
    set_if_present("^v", log(1.0))
    set_if_present("^v_.*TRUE", log(2))
    set_if_present("^v_.*FALSE", log(1))
    set_if_present("^v_RateST", log(2.5))
    set_if_present("^B", log(1))
    set_if_present("^A", log(.4))
    set_if_present("^t0", log(.15))
    set_if_present("^s", log(1.0))
    set_if_present("^muS", log(.40))
    set_if_present("^sigmaS", log(.55))
    set_if_present("^tauS", log(.25))
    set_if_present("^gf", qnorm(.1))
    set_if_present("^tf", qnorm(.05))
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
  print(tapply(Rtfun(d), d$S, function(x) mean(x, na.rm = TRUE)))

}

run_ss_demo <- function(model_type = SSEXG,
                        include_stop_triggered = FALSE,
                        n_trials = 1000,
                        UC = Inf,
                        LT = 0,
                        UT = Inf,
                        ssd_values = c(.15, .25),
                        p_ssd = c(.25, .25),
                        label = NULL) {

  st_label <- "st"
  Rlevels <- if (include_stop_triggered) c("left", "right", st_label) else c("left", "right")
  Slevels <-  c("left", "right")
  
  # Use contrast to set drift/mu parameters
  drift_df <- expand.grid(
    R = Rlevels,
    S = Slevels,
    stringsAsFactors = FALSE
  )
  drift_df$cell_label <- apply(drift_df, 1, function(x) paste(x, collapse = "."))
  drift_cells <- drift_df$cell_label
  rate_names  =  c("lMTRUE_Left", "lMTRUE_Right", "lMFALSE", "ST")
  for (p in rate_names) { drift_df[[p]] <- 0}
  drift_df = drift_df %>%
    mutate(lM=tolower(R)==tolower(S),
           lMTRUE_Left=case_when(R=="left" & lM==TRUE ~ 1, TRUE ~ 0),
           lMTRUE_Right=case_when(R=="right" & lM==TRUE ~ 1, TRUE ~ 0),
           lMFALSE=case_when(!(R==st_label) & lM==FALSE ~ 1, TRUE ~ 0),
           ST=case_when(R==st_label ~ 1, TRUE ~ 0)
    )

  rate_design = drift_df %>%
    select(all_of(rate_names)) %>%
    as.matrix()
  rate_design <- rate_design[, colSums(rate_design != 0) > 0, drop = FALSE]
  rate_names <- colnames(rate_design)
  rownames(rate_design) <- drift_df$cell_label
  RateFun <- function(d)factor(paste(d$lR,d$S,sep="."))
  matchfun <- function(d) as.character(d$S) == as.character(d$lR)
  lIfun <- function(d) {
    if (!include_stop_triggered) {
      return(factor(rep(2, nrow(d)), levels = 1:2))
    }
    factor(ifelse(as.character(d$lR) == st_label, 1, 2), levels = 1:2)
  }
  
  mySSD_function <- function(d) SSD_function(d, SSD = ssd_values, pSSD = p_ssd)

  designSS <- design(
    model = model_type,
    factors = list(subjects = 1, S = Slevels),
    Rlevels = Rlevels,
    matchfun = matchfun,
    functions = list(lI = lIfun, SSD = mySSD_function,Rate=RateFun),
    contrasts = if (model_type()$c_name == "SSEXG") {
      contrasts=list(mu=list(Rate=rate_design))
    } else { 
      contrasts=list(v=list(Rate=rate_design))
    },
    formula = if (model_type()$c_name == "SSEXG") {
      list(mu ~ Rate, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    } else {
      list(v ~ Rate, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
    },
    constants = if (model_type()$c_name == "SSEXG") {
      c(mu=0)
      } else {
        c(v = log(1))
      }
  )
  p_vector <- sampled_pars(designSS, doMap = FALSE)
  p_vector <- set_ss_pars(p_vector, model_type, include_stop_triggered)
  TC = list(UC = UC,LT = LT, UT = UT, LC = LT, verbose=TRUE)
  dat <- make_data(
    p_vector, designSS, n_trials = n_trials,
    TC=TC
  )
  dadm = EMC2:::design_model(dat,designSS)
  cat("\n--- SS demo ---\n")
  if (!is.null(label)) cat(label, "\n")
  cat("model_type =", model_type()$c_name, "\n")
  cat("UC =", UC, " LT =", LT, " UT =", UT, "\n\n")
  summarise_ss_data(dat, include_stop_triggered, Rlevels)

  if (isTRUE(RUN_PROFILE_PLOTS)) {
    # SS R sucks so it's commented out
    #profile_plot_test(dat, designSS, p_vector, n_cores = 1, layout = c(3, 3),
    #                  use_c = FALSE, figure_title = paste("R Likelihood:", label))
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

## Test profile plot of C++ code with no censoring
# --- Test 1: SSexG baseline (no censoring) ------------------------
res_exg <- run_ss_demo(
  model_type = SSEXG,
  include_stop_triggered = FALSE,
  UC = Inf,
  label = "SSexG UC=Inf"
)

mean(res_exg$data$R==res_exg$data$S,na.rm=TRUE)
# --- Test 2: SSexG with deadline censoring -----------------------------------
res_exg_cens <- run_ss_demo(
  model_type = SSEXG,
  include_stop_triggered = FALSE,
  UC = 1.5,
  label = "SSexG UC=1.5"
)

# --- Test 3: SSexG with ST -----------------------------------
res_exg_st <- run_ss_demo(
  model_type = SSEXG,
  include_stop_triggered = TRUE,
  UC = Inf,
  label = "SSexG ST"
)

# --- Test 4: SSexG with deadline censoring and ST-----------------------------------
res_exg_cens_st <- run_ss_demo(
  model_type = SSEXG,
  include_stop_triggered = TRUE,
  UC = 1.5,
  label = "SSexG ST UC=1.5"
)

## Couple of tests with SShybrid model. 
# --- Test 5: SShybrid with deadline censoring --------------------------------
res_hyb <- run_ss_demo(
  model_type = SSRDEX,
  include_stop_triggered = FALSE,
  label = "SShybrid"
)
res_hyb_cens <- run_ss_demo(
  model_type = SSRDEX,
  include_stop_triggered = FALSE,
  UC = 1.5,
  label = "SShybrid UC=1.5"
)

# Tail of Stop Accumulator seems difficult to recover if ST accumulator is present (probably expected)
res_hyb_st <- run_ss_demo(
  model_type = SSRDEX,
  include_stop_triggered = TRUE,
  label = "SShybrid ST"
)
res_hyb_cens_st <- run_ss_demo(
  model_type = SSRDEX,
  include_stop_triggered = TRUE,
  UC = 1.5,
  LT = 0,
  UT = Inf,
  label = "SShybrid ST UC=1.5"
)
