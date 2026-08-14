## RLF Race Model Simulation Recovery Test
rm(list = ls())
library(EMC2)

alpha_mapped = 1.5
n_trials = 1000
label = "rlf_test"
cores_for_chains = 3
cores_per_chain = 10
cat(sprintf("\n========================================================\n"))
cat(sprintf("   Starting RLF Simulation Recovery Test: %s (alpha = %.2f)\n", label, alpha_mapped))
cat(sprintf("========================================================\n\n"))

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designRLF <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RLF(),
  formula = list(v ~ S*lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, alpha ~ 1),
  constants = c(s = log(1),A=log(0),`v_Sright`=0),
  TC = list(UC=2)
)

p_vector <- sampled_pars(designRLF, doMap = FALSE)
p_vector["B"] <- log(2)
p_vector["t0"] <- log(0.2)
p_vector["v"] <- log(1)
p_vector["v_lMTRUE"] <- log(1.5)
p_vector["v_Sright:lMTRUE"] <- log(0.8)
# For alpha: mapped_alpha = 1 + (2-1)*pnorm(alpha_sampled)
# So alpha_sampled = qnorm(alpha_mapped - 1)
p_vector["alpha"] <- qnorm(alpha_mapped - 1)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

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


predData = predict(fit_res)
pdf("fit-RLF.pdf")
plot_cdf(dat, post_predict=predData,
         functions=list(Correct=matchfun), defective_factor = "R", factors = "S")
dev.off()


