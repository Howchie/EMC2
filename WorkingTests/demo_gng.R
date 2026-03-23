
# Demo for Go/No-Go Models
rm(list=ls())
library(EMC2)
source("WorkingTests/test_likelihood_plotfuns_ah.R")


# Whether to run the (slow) MCMC recovery fits
RUN_FITS <- TRUE

# Helper functions
Rtfun <- function(d) {
  out <- d$rt
  out[is.infinite(out)] <- NA
  out
}
Cfun = function(d) dplyr::case_when(d$S=="go" & d$R =="go"~ TRUE,d$S=="nogo" & is.na(d$R) ~ TRUE,
                                          d$S=="go" & is.na(d$R) ~ FALSE,d$S=="nogo" & d$R =="go" ~ FALSE)


#### LBA ----

run_lba_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 500, UC = 3, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,cores_per_chain=3,
                         print_stats=TRUE,
                         label = NULL, n_subj = 1, load = c("L")) {

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designLBA <- design(
    factors = list(subjects = seq(1:n_subj), S = c("go", "nogo", "nogo", "nogo"), L=load),
    Rlevels = c("go", "nogo"),
    matchfun = matchfun,
    functions = list(
      match = function(d) ifelse(d$lM==TRUE, 1, 0),
      mismatch = function(d) ifelse(d$lM==TRUE, 0, 1)
    ),
    model=LBA,UC=3,
    formula = c(
      list(B ~ 0+lR), if(length(load)>1) list(v~0+mismatch + S:match:L) else list(v~0+mismatch + S:match), list(A ~ 1, t0 ~ 1, sv ~ 1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(sv = log(1),A=log(0.4))
  )

  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLBA, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B_lRgo" %in% names(p_vector)) p_vector[["B_lRgo"]] <- log(.75)
  if ("B_lRnogo" %in% names(p_vector)) p_vector[["B_lRnogo"]] <- log(.6) # slight bias for nogo given design
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(.4)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("v_mismatch" %in% names(p_vector)) p_vector[["v_mismatch"]] <- .5
  if ("v_Sgo:match" %in% names(p_vector)) p_vector[["v_Sgo:match"]] <- 2.1
  if ("v_Snogo:match" %in% names(p_vector)) p_vector[["v_Snogo:match"]] <- 2.5
  if ("v_Sgo:match:LL" %in% names(p_vector)) p_vector[["v_Sgo:match:LL"]] <- 2.1
  if ("v_Snogo:match:LL" %in% names(p_vector)) p_vector[["v_Snogo:match:LL"]] <- 2.3
  if ("v_Sgo:match:LM" %in% names(p_vector)) p_vector[["v_Sgo:match:LM"]] <- 1.65
  if ("v_Snogo:match:LM" %in% names(p_vector)) p_vector[["v_Snogo:match:LM"]] <- 1.95
  if ("v_Sgo:match:LH" %in% names(p_vector)) p_vector[["v_Sgo:match:LH"]] <- 1.1
  if ("v_Snogo:match:LH" %in% names(p_vector)) p_vector[["v_Snogo:match:LH"]] <- 1.65

  # Sv for match accumulator
  if ("sv" %in% names(p_vector)) p_vector[["sv"]] <- log(1)
  if ("sv_mismatch" %in% names(p_vector)) p_vector[["sv_mismatch"]] <- log(1)
  if ("sv_match" %in% names(p_vector)) p_vector[["sv_match"]] <- log(.8)
  if ("sv_lMTRUE" %in% names(p_vector)) p_vector[["sv_lMTRUE"]] <- log(.8)
  print(p_vector)
  # Contaminant omissions: pContaminant is on probit scale (transformed via pnorm)
  if (estimate_contaminant) {
    if (!("pContaminant" %in% names(p_vector))) {
      stop("estimate_contaminant=TRUE but pContaminant not found in sampled parameters.")
    }
    if (p_contaminant <= 0 || p_contaminant >= 1) stop("p_contaminant must be in (0,1) when estimated.")
    p_vector[["pContaminant"]] <- qnorm(p_contaminant)
  } else {
    # Fixed to 0 contaminant probability (i.e., pnorm(-Inf) = 0)
    if ("pContaminant" %in% names(p_vector)) p_vector[["pContaminant"]] <- qnorm(0)
  }
  
  if (n_subj>1) {
    # Hierarchical simulation for recovery should use an explicit covariance matrix.
    # The default make_random_effects() scales variances by abs(mean), which can
    # create very large variances and exact zero variances (if a mean is 0).
    # I modified it to baseline at natural mean on the transformed scale zero (e.g. for the default log(1)=0, we'd get a variance of variance_prop * exp(0))
    
    tmp <- make_random_effects(designLBA,p_vector,n_subj=n_subj, variance_proportion = .1)
    dat <- make_data(tmp,designLBA,n_trials = n_trials,
                     TC=list(UC = UC, UT = UT, LC = LC, LT = LT))
  } else {
      dat <- make_data(
      p_vector, designLBA,
      n_trials = n_trials,
      TC=list(UC = UC, UT = UT, LC = LC, LT = LT))
  }

  if (print_stats) {

  cat("\n--- LBA demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)
  
  rtct <- system.time({rtc <- profile_plot_test(dat, designLBA, p_vector, n_cores = cores_for_chains, range=range,
                                                layout = NULL, use_c = TRUE,  figure_title = paste("C++:", label), natural=natural)})
  flat_selection = ifelse(n_subj>1, c("alpha", "subj_ll"), c("alpha", "subj_ll","theta_mu"))
  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBA, type = ifelse(n_subj>1,"standard","single"))
    emc <- fit(emc,stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = flat_selection,
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain=cores_per_chain, cores_for_chains = cores_for_chains), max_tries = 30)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designLBA, true_pars = p_vector,emc=emc,pp=post_predict)))
  }

  invisible(list(data = dat, design = designLBA, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)


# Test 1: Mild upper Censor (GNG must have a finite UC)
# With current design +  parameters burn is a struggle and sampling takes a while
# I think sv:lM, v_match:S, B_lR it's hard to identify, but does converge. Parameters aren't "right" but the posterior predicted data is.
res_lba_gng <- run_lba_demo(
  n_trials = 10000,
  UC = 3,
  load = c("L","M"),
  label = "GNG-LBA UC=3"
)


if (RUN_FITS) print(recovery(res_lba_gng$emc, true_pars = res_lba_gng$true_pars))
if (RUN_FITS) plot_cdf(res_lba_gng$dat, post_predict=res_lba_gng$pp, functions=list(Correct=Cfun), defective_factor = "Correct", factors="S")
if (RUN_FITS) plot_stat(res_lba_gng$dat, post_predict=res_lba_gng$pp, factors="S", stat_name = "MeanCorrect",
                        stat_fun = function(d){mean(d$Correct)}, functions=list(Correct=Cfun))
res_lba_gng_lc <- run_lba_demo(
  n_trials = 10000,
  UC = 1.8,
  LC = .7,
  label = "GNG-LBA UC=1.8, LC=.7"
)

res_lba_gng_hierarchical <- run_lba_demo(
  n_trials = 200,
  UC = 3,
  n_subj=30,
  Load = c("L","M","H"),
  label = "GNG-LBA Hierarchical + Load UC=3"
)

#### RDM ----

run_rdm_demo <- function(p_contaminant = 0,estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,cores_per_chain=3,
                         print_stats=TRUE,
                         label = NULL, n_subj=1, load=c("L")) {

  matchfun <- function(d) tolower(d$S) == tolower(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  designRDM <- design(
    factors = list(subjects = 1, S = c("nogo", "nogo","nogo","go")), # showing the order doesn't matter
    Rlevels = c("go", "nogo"),
    matchfun = matchfun,
    functions = list(
      match = function(d) ifelse(d$lM==TRUE, 1, 0)
    ),
    model = RDM,
    formula = c(
      list(B ~ lR, v ~ match + S, A ~ 1, t0 ~ 1, s ~ 1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(s = log(1), A = log(.4))
  )
  
  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designRDM, doMap = FALSE)

  # Threshold / start-point / non-decision time
  if ("B" %in% names(p_vector)) p_vector[["B"]] <- log(1.6)
  if ("B_lRnogo" %in% names(p_vector)) p_vector[["B_lRnogo"]] <- log(.85)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(.4)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("v" %in% names(p_vector)) p_vector[["v"]] <- log(1)
  if ("v_match" %in% names(p_vector)) p_vector[["v_match"]] <- log(2.5)
  if ("v_Sgo" %in% names(p_vector)) p_vector[["v_Sgo"]] <- log(.9)

  # s for match accumulator
  if ("s" %in% names(p_vector)) p_vector[["s"]] <- log(1)
  if ("s_lMTRUE" %in% names(p_vector)) p_vector[["s_lMTRUE"]] <- log(.8)

  # Contaminant omissions: pContaminant is on probit scale (transformed via pnorm)
  if (estimate_contaminant) {
    if (!("pContaminant" %in% names(p_vector))) {
      stop("estimate_contaminant=TRUE but pContaminant not found in sampled parameters.")
    }
    if (p_contaminant <= 0 || p_contaminant >= 1) stop("p_contaminant must be in (0,1) when estimated.")
    p_vector[["pContaminant"]] <- qnorm(p_contaminant)
  } else {
    # Fixed to 0 contaminant probability (i.e., pnorm(-Inf) = 0)
    if ("pContaminant" %in% names(p_vector)) p_vector[["pContaminant"]] <- qnorm(0)
  }
  print(p_vector)
  dat <- make_data(
    p_vector, designRDM,
    n_trials = n_trials,
    TC=list(UC = UC,UT = UT,LC = LC,LT = LT)
  )

  if (print_stats) {

  cat("\n--- RDM demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  
  #  C _race_cens_trunc
  rtct <- system.time({rtc <- profile_plot_test(dat, designRDM, p_vector, n_cores = cores_for_chains, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood RDM",natural=natural)})
  
  flat_selection = ifelse(n_subj>1, c("alpha", "subj_ll"), c("alpha", "subj_ll","theta_mu"))
  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designRDM, type = ifelse(n_subj>1,"standard","single"))
    emc <- fit(emc,stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = flat_selection,
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      ),cores_per_chain=cores_per_chain, cores_for_chains = cores_for_chains), max_tries = 30)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designRDM, true_pars = p_vector,emc=emc,pp=post_predict)))
  }
  

  invisible(list(data = dat, design = designRDM, true_pars = p_vector))
}

## Parameter estimates aren't accurate but predicted RT x accuracy is

res_rdm_gng <- run_rdm_demo(
  n_trials = 10000,
  UC = 3,
  load = c("L","M"),
  label = "GNG-RDM UC=3"
)


if (RUN_FITS) print(recovery(res_rdm_gng$emc, true_pars = res_rdm_gng$true_pars))
if (RUN_FITS) plot_cdf(res_rdm_gng$dat, post_predict=res_rdm_gng$pp, functions=list(Correct=Cfun), defective_factor = "Correct", factors="S")
if (RUN_FITS) plot_stat(res_rdm_gng$dat, post_predict=res_rdm_gng$pp, factors="S", stat_name = "MeanCorrect",
                        stat_fun = function(d){mean(d$Correct)}, functions=list(Correct=Cfun))

## Heavy censoring seems to make v_match:Snogo unidentifiable because there's no RT shape to distinguish it.
# LC makes gng even slower because there's numerical integration required
res_rdm_gng_lc <- run_rdm_demo(
  n_trials = 10000,
  UC = 4,
  LC = .2,
  label = "GNG-RDM UC=1.8, LC=.7"
)

if (RUN_FITS) print(recovery(res_rdm_gng_lc$emc, true_pars = res_rdm_gng_lc$true_pars))
if (RUN_FITS) plot_cdf(res_rdm_gng_lc$dat, post_predict=res_rdm_gng_lc$pp, functions=list(Correct=Cfun), defective_factor = "Correct", factors="S")
if (RUN_FITS) plot_stat(res_rdm_gng_lc$dat, post_predict=res_rdm_gng_lc$pp, factors="S", stat_name = "MeanCorrect",
                        stat_fun = function(d){mean(d$Correct)}, functions=list(Correct=Cfun))

res_rdm_gng_hierarchical <- run_rdm_demo(
  n_trials = 200,
  UC = 3,
  n_subj=30,
  Load = c("L","M","H"),
  label = "GNG-RDM Hierarchical + Load UC=3"
)

#### LNR ----
run_LNR_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 500, UC = NULL, UT=NULL, LC=NULL, LT=NULL,range=1,
                         cores_for_chains = 3,layout=c(2,3),natural=TRUE,n_cores=1,
                         print_stats=FALSE,
                         sample_file = "samples_LNR.RData") {

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  # Base 2-choice design; lM is automatically constructed from matchfun.
  capture.output(suppressMessages(designLNR <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LNR,
    formula = c(
      list(m~lM,s~lM,t0~1),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
  )),file=NULL)

  # Set simulation parameters (on the transformed scale expected by p_types)
  p_vector <- sampled_pars(designLNR, doMap = FALSE)

  # non-decision time
   if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)

  # Drift means for mismatch vs match accumulators
  if ("m" %in% names(p_vector)) p_vector[["m"]] <- log(1.5)
  if ("m_lMTRUE" %in% names(p_vector)) p_vector[["m_lMTRUE"]] <- log(.45)

  # s for match accumulator
  if ("s" %in% names(p_vector)) p_vector[["s"]] <- log(1)
  if ("s_lMTRUE" %in% names(p_vector)) p_vector[["s_lMTRUE"]] <- log(.8)

  # Contaminant omissions: pContaminant is on probit scale (transformed via pnorm)
  if (estimate_contaminant) {
    if (!("pContaminant" %in% names(p_vector))) {
      stop("estimate_contaminant=TRUE but pContaminant not found in sampled parameters.")
    }
    if (p_contaminant <= 0 || p_contaminant >= 1) stop("p_contaminant must be in (0,1) when estimated.")
    p_vector[["pContaminant"]] <- qnorm(p_contaminant)
  } else {
    # Fixed to 0 contaminant probability (i.e., pnorm(-Inf) = 0)
    if ("pContaminant" %in% names(p_vector)) p_vector[["pContaminant"]] <- qnorm(0)
  }

  dat <- make_data(
    p_vector, designLNR,
    n_trials = n_trials,
    TC=list(UC = UC,UT = UT,LC = LC,LT = LT)
  )

  if (print_stats) {

  cat("\n--- LNR demo ---\n")

  cat("Counts by response (NA = omission):\n")
  print(table(dat$R, useNA = "ifany"))

  cat("\nMean RT by stimulus (finite RTs only):\n")
  print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))

  cat("\nMean Acc by stimulus (defined responses only):\n")
  print(tapply(Cfun(dat), dat$S, function(x) mean(x,na.rm=TRUE)))

  cat("\nOmission rate by stimulus:\n")
  print(tapply(is.na(dat$R), dat$S, mean))

  cat("\nTruncation rate by stimulus:\n")
  print(1-tapply(dat$R, dat$S, length)/n_trials)

  cat("\n")
  }

  # Profile plots. Use n_cores=1 on Windows.
  library(parallel)
  par(mfrow=layout)

  #  C _race_cens_trunc
  rtct <- system.time({rtc <- profile_plot_test(dat, designLNR, p_vector, n_cores = n_cores, range=range,
  layout = NULL, use_c = TRUE, figure_title = "C Likelihood",natural=natural)})

  print(rtc)
  cat(paste0("C likelihood: ",round(rtct[3],2),"\n"))

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLNR, type = "single")
    emc <- fit(emc, cores_for_chains = cores_for_chains, fileName = sample_file)
    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
  }

  invisible(list(data = dat, design = designLNR, true_pars = p_vector))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Test 1: Equivalence of R and C++ function no censor/trunc/contam ----

res <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  sample_file = "samples_LNR_old.RData"
)

# Test 2: Censored RTs ----
res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.4,
  sample_file = "samples_LNR_cens.RData"
)

res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LC = .45,
  n_cores=9,
  sample_file = "samples_LNR_cens.RData"
)


res_cens <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  n_cores=9,
  UC = 1.6,
  LC = .4,
  sample_file = "samples_LNR_cens.RData"
)


# Test 3: Truncated RTs ------------

res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.4,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)


res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  LT = .45,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)


res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UT = 1.5,
  LT = .4,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)

# Test 4: Truncation and censoring ------------

res_contam <- run_LNR_demo(
  p_contaminant = 0,
  estimate_contaminant = FALSE,
  n_trials = 10000,
  UC = 1.5,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_LNR_cens_trunc.RData"
)

# Test 5: Contaminant omissions  ------------

res_contam <- run_LNR_demo(
  p_contaminant = 0.15,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  sample_file = "samples_LNR_contam.RData"
)


res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.5,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

res_contam <- run_LNR_demo(
  p_contaminant = .05,
  estimate_contaminant = TRUE,
  n_trials = 10000,
  UC = 1.5,
  LC = .425,
  UT = 1.7,
  LT = .375,
  n_cores=9,
  sample_file = "samples_rdm_cens_trunc.RData"
)

