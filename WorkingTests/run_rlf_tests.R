## Main script to run RLF simulation recovery tests for low alpha (1.1) and high alpha (1.8)

source("WorkingTests/test_rlf_recovery_demo.R")

cat("\n=================================================================\n")
cat(" RUNNING RLF SIMULATION RECOVERY TEST 1: LOW ALPHA = 1.1\n")
cat("=================================================================\n\n")

time_low_alpha <- system.time({
  res_alpha_1_1 <- run_rlf_recovery(
    alpha_mapped = 1.1,
    n_trials = 2000,
    label = "rlf_alpha_1_1"
  )
})

cat("\nTotal elapsed time for alpha = 1.1 test (data gen + MCMC fit + recovery):\n")
print(time_low_alpha)


cat("\n=================================================================\n")
cat(" RUNNING RLF SIMULATION RECOVERY TEST 2: HIGH ALPHA = 1.8\n")
cat("=================================================================\n\n")

time_high_alpha <- system.time({
  res_alpha_1_8 <- run_rlf_recovery(
    alpha_mapped = 1.8,
    n_trials = 2000,
    label = "rlf_alpha_1_8"
  )
})

cat("\nTotal elapsed time for alpha = 1.8 test (data gen + MCMC fit + recovery):\n")
print(time_high_alpha)

cat("\n=================================================================\n")
cat(" SUMMARY COMPARISON OF RECOVERY RESULTS\n")
cat("=================================================================\n")

cat("\n--- Test 1: Alpha = 1.1 Timing & Recovery ---\n")
cat("System Time:\n")
print(time_low_alpha)
cat("\nMapped Parameter Recovery (alpha = 1.1):\n")
print(round(res_alpha_1_1$mapped_summary, 4))

cat("\n--- Test 2: Alpha = 1.8 Timing & Recovery ---\n")
cat("System Time:\n")
print(time_high_alpha)
cat("\nMapped Parameter Recovery (alpha = 1.8):\n")
print(round(res_alpha_1_8$mapped_summary, 4))

save(res_alpha_1_1, res_alpha_1_8, time_low_alpha, time_high_alpha,
     file = "WorkingTests/rlf_recovery_results.RData")
cat("\nResults saved to WorkingTests/rlf_recovery_results.RData\n")
