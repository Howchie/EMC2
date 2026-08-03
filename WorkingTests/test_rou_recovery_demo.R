## RLF Race Model Simulation Recovery Test
## Testing low alpha (1.1) and high alpha (1.8) single-subject single-condition recovery

rm(list = ls())
library(EMC2)

leak = 1.5
n_trials = 10000
label = "rou_test"
cores_for_chains = 3
cores_per_chain = 8
cat(sprintf("\n========================================================\n"))
cat(sprintf("   Starting RLF Simulation Recovery Test: %s (alpha = %.2f)\n", label, leak))
cat(sprintf("========================================================\n\n"))

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designROU <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = ROU(parameterization="equilibrium"),
  formula = list(theta ~ lM, B ~ 1, A ~ 1, t0 ~ 1, chi ~ 1, tk ~ 1),
  constants = c(B = log(1),A=log(0))
)

designRDM <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RDMSWTN(),
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
  constants = c(s = log(1),A=log(0))
)

p_vector <- sampled_pars(designROU, doMap = FALSE)
p_vector["t0"] <- log(0.15)
p_vector["chi"] <- log(1)
p_vector["theta"] <- log(.5)
p_vector["theta_lMTRUE"] <- log(2.5)
p_vector["tk"] <- log(.5)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

cat("Generating simulated dataset...\n")
t_sim <- system.time({
  dat <- make_data(p_vector, designROU, n_trials = n_trials)
})
cat(sprintf("Data generation time: %.3f sec\n", t_sim["elapsed"]))

cat("Data summary:\n")
cat("Counts by response:\n")
print(table(dat$R, dat$S,useNA = "ifany"))
cat("Mean RT by stimulus:\n")
print(tapply(dat$rt, dat$S, function(x) mean(x[is.finite(x)], na.rm = TRUE)))
plot(density(dat$rt))
emc <- make_emc(dat, designROU, type = "single")
emc2 <- make_emc(dat, designRDM, type = "single", rt_resolution = 1/60)
cat("\nFitting ROU model (MCMC)...\n")
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
cat("\nFitting RDM model (MCMC)...\n")
t_fit <- system.time({
  fit_res2 <- fit(
    emc2,
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

predData = predict(fit_res)
predData2 = predict(fit_res2)
pdf("fit-ROU.pdf")
plot_cdf(dat, post_predict=predData,
         functions=list(Correct=matchfun), defective_factor = "R", factors = "S")
dev.off()

invisible(list(
  data = dat,
  pred = predData,
  design = designRLF,
  true_pars = p_vector,
  fit = fit_res,
  system_time = t_fit,
  mapped_summary = mapped_summary
))

