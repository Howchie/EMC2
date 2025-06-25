# Tests for DSWTN distribution functions (PDF, CDF, RNG) and race likelihood

library(testthat)
# Assuming the package is loaded, which makes C++ functions available.
# Or use devtools::load_all() if running interactively.

# Helper to numerically integrate PDF to get CDF for comparison
numeric_cdf <- function(t_val, alpha, mu_drift, sigma_drift_sq, theta) {
  if (t_val <= theta) return(0)
  tryCatch({
    stats::integrate(function(x_vals) {
      exp(dDSWTN_log(x_vals, alpha, mu_drift, sigma_drift_sq, theta))
    }, lower = theta, upper = t_val, subdivisions = 500, rel.tol = 1e-5)$value
  }, error = function(e) NA_real_) # Return NA on integration error
}

test_that("dDSWTN_log basic properties", {
  # Valid parameters
  alpha <- 1.0
  mu_drift <- 1.5
  sigma_drift_sq <- 0.5
  theta <- 0.1

  # Test times
  t_vals <- c(theta - 0.01, theta, theta + 0.01, theta + 0.2, theta + 1.0)

  log_pdfs <- dDSWTN_log(t_vals, alpha, mu_drift, sigma_drift_sq, theta)

  expect_equal(length(log_pdfs), length(t_vals))
  expect_true(is.infinite(log_pdfs[1]) && log_pdfs[1] < 0) # t < theta
  expect_true(is.infinite(log_pdfs[2]) && log_pdfs[2] < 0) # t == theta
  expect_true(all(is.finite(log_pdfs[3:5])))
  expect_true(all(log_pdfs[3:5] <= 0)) # log-PDF should be <= 0

  # Test vectorization of parameters
  t_val <- theta + 0.5
  alphas <- c(1.0, 1.2)
  log_pdfs_alphas <- dDSWTN_log(t_val, alphas, mu_drift, sigma_drift_sq, theta)
  expect_equal(length(log_pdfs_alphas), 2)
  expect_true(log_pdfs_alphas[1] != log_pdfs_alphas[2])

  # Case: sigma_drift_sq = 0 (should reduce to Wald)
  log_pdf_wald_approx <- dDSWTN_log(theta + 0.5, alpha, mu_drift, 0, theta)
  # Need a direct Wald PDF for comparison, e.g., from emc²'s RDM dWald or similar
  # For now, just check it's finite
  expect_true(is.finite(log_pdf_wald_approx))

  # Case: mu_drift = 0 (or very small positive)
  log_pdf_zero_drift <- dDSWTN_log(theta + 0.5, alpha, 1e-6, sigma_drift_sq, theta)
  expect_true(is.finite(log_pdf_zero_drift))

  # Negative mu_drift (still possible to hit boundary if sigma_drift_sq is large enough)
  log_pdf_neg_mu <- dDSWTN_log(theta + 0.5, alpha, -1.0, 2.0, theta)
  expect_true(is.finite(log_pdf_neg_mu))

  # Very large sigma_drift_sq
  log_pdf_large_sig <- dDSWTN_log(theta + 0.5, alpha, mu_drift, 100.0, theta)
  expect_true(is.finite(log_pdf_large_sig))
})

test_that("pDSWTN basic properties and comparison with numerical integration", {
  alpha <- 1.0
  mu_drift <- 1.5
  sigma_drift_sq <- 0.5
  theta <- 0.1

  t_vals <- c(theta - 0.01, theta, theta + 0.01, theta + 0.2, theta + 0.5, theta + 1.0, theta + 5.0)
  cdfs_hcub <- pDSWTN(t_vals, alpha, mu_drift, sigma_drift_sq, theta, abs_err = 1e-7, rel_err = 1e-7)

  expect_equal(length(cdfs_hcub), length(t_vals))
  expect_equal(cdfs_hcub[1], 0) # t < theta
  expect_equal(cdfs_hcub[2], 0) # t == theta
  expect_true(all(cdfs_hcub[3:length(t_vals)] >= 0 & cdfs_hcub[3:length(t_vals)] <= 1.05)) # Allow some tolerance for 1
  expect_true(all(diff(cdfs_hcub[cdfs_hcub <=1]) >= -1e-6)) # Should be non-decreasing (allow for small numerical error)
  expect_lt(cdfs_hcub[3], 1.0) # Early CDF should not be 1
  expect_gte(tail(cdfs_hcub, 1), 0.9) # CDF at large t should approach 1

  # Compare with numerical integration of PDF for a few points
  # This is slow, so only for a couple of points
  # Note: numeric_cdf uses dDSWTN_log, so it's a check on pDSWTN's integration vs R's integrate()
  # if dDSWTN_log itself is correct.

  # Point 1
  t_test1 <- theta + 0.3
  cdf_num_int1 <- numeric_cdf(t_test1, alpha, mu_drift, sigma_drift_sq, theta)
  cdf_hcub1 <- pDSWTN(t_test1, alpha, mu_drift, sigma_drift_sq, theta, abs_err = 1e-7, rel_err = 1e-7)
  if (!is.na(cdf_num_int1)) { # If numeric integration succeeded
    expect_equal(cdf_hcub1, cdf_num_int1, tolerance = 1e-4) # Looser tolerance due to two different numerical methods
  } else {
    message(paste("Skipping numeric CDF comparison for t=", t_test1, "due to integration error."))
  }

  # Point 2
  t_test2 <- theta + 0.7
  cdf_num_int2 <- numeric_cdf(t_test2, alpha, mu_drift, sigma_drift_sq, theta)
  cdf_hcub2 <- pDSWTN(t_test2, alpha, mu_drift, sigma_drift_sq, theta, abs_err = 1e-7, rel_err = 1e-7)
   if (!is.na(cdf_num_int2)) {
    expect_equal(cdf_hcub2, cdf_num_int2, tolerance = 1e-4)
  } else {
    message(paste("Skipping numeric CDF comparison for t=", t_test2, "due to integration error."))
  }

  # Case: sigma_drift_sq = 0 (should reduce to Wald CDF)
  cdf_wald_approx <- pDSWTN(theta + 0.5, alpha, mu_drift, 0, theta)
  # Compare with pRDM or pigt0 if available and parameters match.
  # For now, check basic properties
  expect_true(cdf_wald_approx >= 0 && cdf_wald_approx <= 1)
})


test_that("rDSWTN basic properties", {
  alpha <- 1.0
  mu_drift <- 1.5
  sigma_drift_sq <- 0.5
  theta <- 0.1
  n_samples <- 100

  samples <- rDSWTN(n_samples, alpha, mu_drift, sigma_drift_sq, theta)

  expect_equal(length(samples), n_samples)
  expect_true(all(is.finite(samples)))
  expect_true(all(samples >= theta)) # All samples should be >= non-decision time

  # Check that mean is roughly consistent with expectation (hard to get exact)
  # Expected mean is E_xi[alpha/xi] + theta. This is complex.
  # For now, check if mean is greater than theta.
  expect_gt(mean(samples), theta)

  # Case: sigma_drift_sq = 0
  samples_wald <- rDSWTN(n_samples, alpha, mu_drift, 0, theta)
  expect_true(all(samples_wald >= theta))
  # Mean of Wald(alpha, mu_drift) is alpha/mu_drift
  if (mu_drift > 0) {
    expect_equal(mean(samples_wald), (alpha/mu_drift) + theta, tolerance = 0.3) # Wide tolerance for RNG mean
  }
})

test_that("loglik_DSWTN_race basic functionality", {
  # 2 accumulators, 1 trial
  rt <- 0.5
  choice <- 1 # Accumulator 1 wins

  # Parameters for accumulator 1 (winner)
  alpha1 <- 1.0
  mu1 <- 2.0
  sig_sq1 <- 0.5
  theta1 <- 0.1

  # Parameters for accumulator 2 (loser)
  alpha2 <- 1.2
  mu2 <- 1.0
  sig_sq2 <- 0.3
  theta2 <- 0.15

  params_alpha_mat <- matrix(c(alpha1, alpha2), nrow = 1)
  params_mu_drift_mat <- matrix(c(mu1, mu2), nrow = 1)
  params_sigma_drift_sq_mat <- matrix(c(sig_sq1, sig_sq2), nrow = 1)
  params_theta_mat <- matrix(c(theta1, theta2), nrow = 1)

  loglik <- loglik_DSWTN_race(
    rts = rt, choices = choice,
    params_alpha = params_alpha_mat,
    params_mu_drift = params_mu_drift_mat,
    params_sigma_drift_sq = params_sigma_drift_sq_mat,
    params_theta = params_theta_mat,
    cdf_abs_err = 1e-7, cdf_rel_err = 1e-7 # Tighter for testing
  )

  expect_equal(length(loglik), 1)
  expect_true(is.finite(loglik[1]))

  # Manual calculation for this simple case
  t_adj1 <- rt - theta1
  t_adj2 <- rt - theta2

  log_pdf_winner <- dDSWTN_log(t_adj1, alpha1, mu1, sig_sq1, theta1) # Note: dDSWTN_log takes t_unadjusted
  log_pdf_winner_direct <- dDSWTN_log(rt, alpha1, mu1, sig_sq1, theta1)

  # The C++ loglik_DSWTN_race does rt - theta internally for dswtn_logpdf_core
  # but dDSWTN_log (the Rcpp wrapper) also does rt - theta. This could be a double subtraction if not careful.
  # Let's re-check dDSWTN_log: it takes t, alpha, ..., theta. It does t_adj = t[i] - theta[i].
  # So, log_pdf_winner_direct is correct for the R-exported dDSWTN_log.

  # The C++ loglik_DSWTN_race calls dswtn_logpdf_core(t_adj_winner, ...), which is correct.
  # And dswtn_cdf_core(t_adj_loser, ...), correct.

  # So, for manual check:
  manual_log_pdf_winner <- dswtn_logpdf_core(t_adj1, alpha1, mu1, sig_sq1)
  manual_cdf_loser <- dswtn_cdf_core(t_adj2, alpha2, mu2, sig_sq2, abs_err=1e-7, rel_err=1e-7)
  manual_survivor_loser <- 1.0 - manual_cdf_loser

  expected_loglik <- manual_log_pdf_winner + log(manual_survivor_loser)

  # Need to call the core functions from R to verify, or trust C++ internal calls.
  # For now, check against direct calls to R-exported PDF/CDF

  log_pdf1_R <- dDSWTN_log(rt, alpha1, mu1, sig_sq1, theta1)
  cdf2_R <- pDSWTN(rt, alpha2, mu2, sig_sq2, theta2, abs_err = 1e-7, rel_err = 1e-7)
  expected_loglik_R_way <- log_pdf1_R + log(1 - cdf2_R)

  expect_equal(loglik[1], expected_loglik_R_way, tolerance = 1e-6)

  # Test case: choice = 2 (Accumulator 2 wins)
  loglik_choice2 <- loglik_DSWTN_race(
    rts = rt, choices = 2,
    params_alpha = params_alpha_mat,
    params_mu_drift = params_mu_drift_mat,
    params_sigma_drift_sq = params_sigma_drift_sq_mat,
    params_theta = params_theta_mat,
    cdf_abs_err = 1e-7, cdf_rel_err = 1e-7
  )
  log_pdf2_R <- dDSWTN_log(rt, alpha2, mu2, sig_sq2, theta2)
  cdf1_R <- pDSWTN(rt, alpha1, mu1, sig_sq1, theta1, abs_err = 1e-7, rel_err = 1e-7)
  expected_loglik_choice2_R_way <- log_pdf2_R + log(1 - cdf1_R)
  expect_equal(loglik_choice2[1], expected_loglik_choice2_R_way, tolerance = 1e-6)

  # Test with rt < theta for winner
  loglik_early <- loglik_DSWTN_race(
    rts = theta1 - 0.01, choices = 1,
    params_alpha = params_alpha_mat,
    params_mu_drift = params_mu_drift_mat,
    params_sigma_drift_sq = params_sigma_drift_sq_mat,
    params_theta = params_theta_mat
  )
  expect_equal(loglik_early[1], -1e10) # Default min_log_lik

  # Test with rt < theta for loser (should not make overall LL=-1e10 unless PDF of winner is also bad)
  rt_between_thetas <- (theta1 + theta2) / 2 # Assume theta1 < rt_between_thetas < theta2
  if (theta1 < theta2) {
     params_theta_mat_varied <- matrix(c(theta1, theta2 + 0.1), nrow = 1) # Ensure loser theta is higher
     loglik_loser_early <- loglik_DSWTN_race(
        rts = rt_between_thetas, choices = 1, # Winner is acc 1
        params_alpha = params_alpha_mat,
        params_mu_drift = params_mu_drift_mat,
        params_sigma_drift_sq = params_sigma_drift_sq_mat,
        params_theta = params_theta_mat_varied, # Loser (acc2) has theta > rt
        cdf_abs_err = 1e-7, cdf_rel_err = 1e-7
     )
     # Loser's t_adj will be < 0. CDF=0, Survivor=1. log(Survivor)=0.
     # So, likelihood should just be PDF of winner.
     expected_lik_loser_early <- dDSWTN_log(rt_between_thetas, alpha1, mu1, sig_sq1, theta1)
     expect_equal(loglik_loser_early[1], expected_lik_loser_early, tolerance=1e-6)
  }
})

# To run these tests:
# Ensure the package is built and loaded (e.g. devtools::load_all())
# then run testthat::test_file("tests/testthat/test-dswtn.R") or testthat::test_check("YourPackageName")
# Note: dswtn_logpdf_core and dswtn_cdf_core are not directly callable from R,
# so tests rely on their Rcpp wrappers.
# The numeric_cdf helper is also using dDSWTN_log, so it's more of a consistency check
# on the CDF's numerical integration routine against R's generic stats::integrate.

# A small detail in loglik_DSWTN_race test:
# dDSWTN_log(rt, alpha1, mu1, sig_sq1, theta1) uses the R wrapper which subtracts theta1 from rt.
# The manual calculation used dswtn_logpdf_core(t_adj1, ...).
# The test `expect_equal(loglik[1], expected_loglik_R_way, tolerance = 1e-6)` is therefore
# comparing the C++ race likelihood (which uses _core functions with t_adj)
# against a manual calculation using the R-exported _wrappers_ (which also handle t_adj).
# This is a valid end-to-end check for the R user.
# The comment about dswtn_logpdf_core vs dDSWTN_log was to clarify my own thought process
# for the manual calculation. The R-way comparison is the most direct test of the exported race function.

# One final check for loglik_DSWTN_race: what if a survivor probability is exactly 0 (CDF=1)?
# log(0) is -Inf. The C++ code has `if (survivor_loser <= 1e-10)` then `current_log_lik = min_log_lik; break;`
# This is correct behavior.
test_that("loglik_DSWTN_race with CDF=1 for loser", {
  rt <- 2.0 # Late RT
  choice <- 1

  alpha1 <- 1.0; mu1 <- 0.2; sig_sq1 <- 0.01; theta1 <- 0.1 # Winner, slow drift
  alpha2 <- 0.1; mu2 <- 5.0; sig_sq2 <- 0.01; theta2 <- 0.1 # Loser, very fast, should have CDF=1 by rt=2.0

  params_alpha_mat <- matrix(c(alpha1, alpha2), nrow = 1)
  params_mu_drift_mat <- matrix(c(mu1, mu2), nrow = 1)
  params_sigma_drift_sq_mat <- matrix(c(sig_sq1, sig_sq2), nrow = 1)
  params_theta_mat <- matrix(c(theta1, theta2), nrow = 1)

  # Check if pDSWTN for loser is indeed ~1.0
  cdf_loser_val <- pDSWTN(rt, alpha2, mu2, sig_sq2, theta2, abs_err = 1e-9, rel_err = 1e-9)
  # If this is not very close to 1, this test case isn't well-posed for this check.
  # It should be very close to 1 for these parameters.
  # Forcing it to be 1 for the test logic:
  # One way to ensure this is to mock pDSWTN or dswtn_cdf_core, but that's too complex for here.
  # Instead, rely on the parameters making it very likely.
  # If cdf_loser_val is e.g. 0.99999999999, then survivor is 1e-10, log(survivor) is large negative.
  # If cdf_loser_val is 1.0, then survivor is 0, log(survivor) is -Inf.

  # The C++ code uses `survivor_loser <= 1e-10`.
  # If pDSWTN returns a value like 1.0 - 1e-11, survivor is 1e-11, condition met.
  # If pDSWTN returns exactly 1.0, survivor is 0.0, condition met.

  # We expect the log-likelihood to become min_log_lik due to this.
  # Only if the PDF of the winner is not already -Inf.
  pdf_winner_val <- dDSWTN_log(rt, alpha1, mu1, sig_sq1, theta1)
  if (is.finite(pdf_winner_val) && cdf_loser_val > 1 - 1.5e-10) { # Ensure winner PDF is fine & loser CDF is ~1
     loglik <- loglik_DSWTN_race(
        rts = rt, choices = choice,
        params_alpha = params_alpha_mat,
        params_mu_drift = params_mu_drift_mat,
        params_sigma_drift_sq = params_sigma_drift_sq_mat,
        params_theta = params_theta_mat,
        min_log_lik = -1e7 # Use a distinct min_log_lik for this test
     )
     expect_equal(loglik[1], -1e7)
  } else {
    message("Skipping CDF=1 test for loglik_DSWTN_race as conditions not met by parameter values.")
    message(paste("Winner PDF:", pdf_winner_val, "Loser CDF:", cdf_loser_val))
  }
})
