## Run seed=123 with rt_resolution=1/60 AND informed prior on v_lMTRUE ~ N(2,1)
## Does NOT source the full demo (which has 15 inline fit calls).
## Inlines only the run_lba_demo definition from demo_cens_trunc.R.

rm(list = ls())
library(EMC2)
library(parallel)
source("WorkingTests/test_likelihood_plotfuns_ah.R")

RUN_FITS <- TRUE

Rtfun <- function(d) { out <- d$rt; out[is.infinite(out)] <- NA; out }
Cfun  <- function(d) as.numeric(d$S) == as.numeric(d$R)

run_lba_demo <- function(p_contaminant = 0, estimate_contaminant = FALSE,
                         n_trials = 500, UC = Inf, UT = Inf, LC = 0, LT = 0, range = 1,
                         cores_for_chains = 3, layout = c(2, 3), natural = TRUE,
                         cores_per_chain = 3,
                         rt_resolution = 1/60,
                         prior = NULL,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         force_direction = TRUE, force_response = TRUE,
                         print_stats = TRUE,
                         label = NULL, posdrift = TRUE, seed = NULL) {

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  designLBA <- design(
    factors  = list(subjects = 1, S = c("left", "right")),
    Rlevels  = c("left", "right"),
    matchfun = matchfun,
    model    = LBA(posdrift = posdrift),
    formula  = c(
      list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ 0 + lM),
      if (estimate_contaminant) list(pContaminant ~ 1) else list()
    ),
    constants = c(A = log(0.4))
  )

  p_vector <- sampled_pars(designLBA, doMap = FALSE)
  if ("B"         %in% names(p_vector)) p_vector[["B"]]         <- log(1.2)
  if ("A"         %in% names(p_vector)) p_vector[["A"]]         <- log(0.4)
  if ("t0"        %in% names(p_vector)) p_vector[["t0"]]        <- log(0.15)
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.4
  if ("v_lMTRUE"  %in% names(p_vector)) p_vector[["v_lMTRUE"]]  <- 1.2
  if ("sv_lMFALSE" %in% names(p_vector)) p_vector[["sv_lMFALSE"]] <- log(1.2)
  if ("sv_lMTRUE"  %in% names(p_vector)) p_vector[["sv_lMTRUE"]]  <- log(1)

  if (estimate_contaminant) {
    if (!("pContaminant" %in% names(p_vector)))
      stop("estimate_contaminant=TRUE but pContaminant not found in sampled parameters.")
    if (p_contaminant <= 0 || p_contaminant >= 1)
      stop("p_contaminant must be in (0,1) when estimated.")
    p_vector[["pContaminant"]] <- qnorm(p_contaminant)
  } else {
    if ("pContaminant" %in% names(p_vector)) p_vector[["pContaminant"]] <- qnorm(0)
  }

  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  dat <- make_data(
    p_vector, designLBA,
    n_trials = n_trials,
    TC = list(UC = UC, UT = UT, LC = LC, LT = LT,
              LCresponse = LCresponse, UCresponse = UCresponse,
              LCdirection = LCdirection, UCdirection = UCdirection)
  )

  if (print_stats) {
    cat("\n--- LBA demo ---\n")
    cat("Counts by response (NA = omission):\n")
    print(table(dat$R, useNA = "ifany"))
    cat("\nMean RT by stimulus (finite RTs only):\n")
    print(tapply(Rtfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))
    cat("\nMean Acc by stimulus (defined responses only):\n")
    print(tapply(Cfun(dat), dat$S, function(x) mean(x, na.rm = TRUE)))
    cat("\nOmission rate by stimulus:\n")
    print(tapply(is.na(dat$R), dat$S, mean))
    cat("\nTruncation rate by stimulus:\n")
    print(1 - (tapply(dat$R, dat$S, length) / n_trials))
    cat("\n")
  }

  par(mfrow = layout)
  rtct <- system.time({
    rtc <- profile_plot_test(dat, designLBA, p_vector,
                             n_cores = cores_for_chains, range = range,
                             layout = NULL, use_c = TRUE,
                             figure_title = paste("C++:", label), natural = natural)
  })

  if (isTRUE(RUN_FITS)) {
    emc <- make_emc(dat, designLBA, type = "single",
                    rt_resolution = rt_resolution, prior = prior)

    emc <- fit(emc,
      fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),
      stop_criteria = list(
        adapt = list(
          min_unique = 1500,   # 2.5x default (600) for longer adaptation
          max_gd     = 1.10
        ),
        sample = list(
          iter            = 1000,
          max_gd          = 1.10,
          max_flat_loc    = 0.5,
          flat_selection  = c("alpha", "subj_ll"),
          flat_p1         = 1/3,
          flat_p2         = 1/3,
          max_sample_iter = 5000
        )
      ),
      cores_per_chain  = cores_per_chain,
      cores_for_chains = cores_for_chains,
      max_tries = 30
    )

    post_predict <- predict(emc, n_post = 50)
    plot_pars(emc, post_predict = post_predict, true_pars = p_vector)
    return(invisible(list(data = dat, design = designLBA, true_pars = p_vector,
                          emc = emc, pp = post_predict)))
  }

  invisible(list(data = dat, design = designLBA, true_pars = p_vector))
}

# ── Build informed prior: v_lMTRUE centred on 2 instead of 0 ─────────────────
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
designLBA <- design(
  factors  = list(subjects = 1, S = c("left", "right")),
  Rlevels  = c("left", "right"),
  matchfun = matchfun,
  model    = LBA(posdrift = TRUE),
  formula  = list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1, sv ~ 0 + lM),
  constants = c(A = log(0.4))
)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

res <- run_lba_demo(
  n_trials         = 10000,
  label            = "no_cens_trunc_8chains",
  seed             = 123,
  rt_resolution    = 1/60,
  prior            = NULL,   # default N(0,1)
  print_stats      = FALSE,
  cores_for_chains = 8,
  cores_per_chain  = 1
)

if (isTRUE(RUN_FITS)) {
  free_pars <- res$true_pars[names(res$true_pars) != "A"]
  cat("\n=== RECOVERY: seed=123, rt_resolution=1/60, 8 chains, long adapt ===\n")
  print(recovery(res$emc, true_pars = free_pars))
}
