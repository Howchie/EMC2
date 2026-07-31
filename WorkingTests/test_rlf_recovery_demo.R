## RLF Race Model Simulation Recovery Test
## Testing low alpha (1.1) and high alpha (1.8) single-subject single-condition recovery

rm(list = ls())
library(EMC2)

run_rlf_recovery <- function(alpha_mapped = 1.1, n_trials = 1000, label = "rlf_test", cores_for_chains = 3, cores_per_chain = 1) {
  cat(sprintf("\n========================================================\n"))
  cat(sprintf("   Starting RLF Simulation Recovery Test: %s (alpha = %.2f)\n", label, alpha_mapped))
  cat(sprintf("========================================================\n\n"))
  
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  designRLF <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = RLF(),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, alpha ~ 1),
    constants = c(s = log(1))
  )

  p_vector <- sampled_pars(designRLF, doMap = FALSE)
  p_vector["B"] <- log(1.5)
  p_vector["A"] <- log(0.5)
  p_vector["t0"] <- log(0.2)
  p_vector["v"] <- log(1)
  p_vector["v_lMTRUE"] <- log(2) - log(1)
  
  # For alpha: mapped_alpha = 1 + (2-1)*pnorm(alpha_sampled)
  # So alpha_sampled = qnorm(alpha_mapped - 1)
  p_vector["alpha"] <- qnorm(alpha_mapped - 1)
  
  set.seed(42 + round(alpha_mapped * 100))
  
  cat("Generating simulated dataset...\n")
  t_sim <- system.time({
    dat <- make_data(p_vector, designRLF, n_trials = n_trials)
  })
  cat(sprintf("Data generation time: %.3f sec\n", t_sim["elapsed"]))
  
  cat("Data summary:\n")
  cat("Counts by response:\n")
  print(table(dat$R, useNA = "ifany"))
  cat("Mean RT by stimulus:\n")
  print(tapply(dat$rt, dat$S, function(x) mean(x[is.finite(x)], na.rm = TRUE)))
  
  emc <- make_emc(dat, designRLF, type = "single")
  
  cat("\nFitting RLF model (MCMC)...\n")
  t_fit <- system.time({
    fit_res <- fit(
      emc,
      fileName = paste0("samples_ss_", gsub("[^a-z0-9]", "_", tolower(label)), ".RData"),
      cores_for_chains = cores_for_chains,
      cores_per_chain = cores_per_chain,
      stop_criteria = list(
        sample = list(
          iter = 1000,
          max_gd = 1.10,
          max_flat_loc = 0.5,
          flat_selection = c("alpha", "subj_ll"),
          flat_p1 = 1/3,
          flat_p2 = 1/3,
          max_sample_iter = 4000
        )
      ),
      max_tries = 30
    )
  })
  
  cat(sprintf("\nFit Completed! System time for fit (%s):\n", label))
  print(t_fit)
  
  cat("\n--- Parameter Recovery (Sampled Scale) ---\n")
  rec_stats <- recovery(fit_res, true_pars = p_vector)
  print(rec_stats)
  
  # Calculate mapped scale parameter recovery summary
  cat("\n--- Mapped (Natural Scale) Parameter Comparison ---\n")
  post_arr <- get_pars(fit_res, selection = "alpha", merge_chains = TRUE, return_mcmc = FALSE)
  post_samples <- t(post_arr[, 1, ])
  
  mapped_true <- c(
    v_mismatch = 1.0,
    v_match = 2.0,
    B = 1.5,
    A = 0.5,
    t0 = 0.2,
    alpha = alpha_mapped
  )
  
  mapped_samples <- matrix(NA, nrow = nrow(post_samples), ncol = 6,
                           dimnames = list(NULL, names(mapped_true)))
  mapped_samples[, "v_mismatch"] <- exp(post_samples[, "v"])
  mapped_samples[, "v_match"]    <- exp(post_samples[, "v"] + post_samples[, "v_lMTRUE"])
  mapped_samples[, "B"]          <- exp(post_samples[, "B"])
  mapped_samples[, "A"]          <- exp(post_samples[, "A"])
  mapped_samples[, "t0"]         <- exp(post_samples[, "t0"])
  mapped_samples[, "alpha"]      <- 1 + pnorm(post_samples[, "alpha"])
  
  mapped_summary <- data.frame(
    True = mapped_true,
    Mean = colMeans(mapped_samples),
    Median = apply(mapped_samples, 2, median),
    SD = apply(mapped_samples, 2, sd),
    q2.5 = apply(mapped_samples, 2, quantile, probs = 0.025),
    q97.5 = apply(mapped_samples, 2, quantile, probs = 0.975),
    Bias = colMeans(mapped_samples) - mapped_true
  )
  print(round(mapped_summary, 4))
  
  invisible(list(
    data = dat,
    design = designRLF,
    true_pars = p_vector,
    fit = fit_res,
    system_time = t_fit,
    mapped_summary = mapped_summary
  ))
}

cat("Script loaded successfully.\n")
