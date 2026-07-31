## Run RLF simulation recovery test for alpha = 1.8

source("WorkingTests/test_rlf_recovery_demo.R")

cat("\n=================================================================\n")
cat(" RUNNING RLF SIMULATION RECOVERY TEST 2: HIGH ALPHA = 1.8\n")
cat("=================================================================\n\n")

time_high_alpha <- system.time({
  res_alpha_1_8 <- run_rlf_recovery(
    alpha_mapped = 1.8,
    n_trials = 1000,
    label = "rlf_alpha_1_8"
  )
})

cat("\nTotal elapsed time for alpha = 1.8 test (data gen + MCMC fit + recovery):\n")
print(time_high_alpha)

save(res_alpha_1_8, time_high_alpha, file = "WorkingTests/rlf_alpha_1_8_results.RData")
cat("\nResults saved to WorkingTests/rlf_alpha_1_8_results.RData\n")
