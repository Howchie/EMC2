# Extracted from test-sampling-rejection.R:10

# test -------------------------------------------------------------------------
samples <- list(
    theta_mu = matrix(1:6, nrow = 2),
    theta_var = array(1:12, dim = c(2, 2, 3)),
    last_theta_var_inv = diag(2),
    stage = c("init", "burn", "burn"),
    idx = 2
  )
out <- reject_sample_iteration(samples, 3)
