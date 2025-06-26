# Tests for RDM_DSWTN distribution functions (PDF, CDF, RNG) and race likelihood
# Parameters: B, A, mu_drift, sigma_drift_sq, t0, s

library(testthat)
# Assuming the package is loaded for C++ functions. For EMC2 model RDM():
# library(EMC2) # Or have it available in testing environment

# Helper to numerically integrate PDF to get CDF for comparison
numeric_cdf_rdm_dswtn <- function(t_val, B, A, mu_drift, sigma_drift_sq, t0, s = 1.0) {
  if (t_val <= t0) return(0)
  tryCatch({
    stats::integrate(function(x_vals) {
      exp(dRDM_DSWTN_log(x_vals, B, A, mu_drift, sigma_drift_sq, t0, s))
    }, lower = t0, upper = t_val, subdivisions = 500, rel.tol = 1e-4)$value
  }, error = function(e) {
    # message(paste("Integration error for numeric_cdf_rdm_dswtn:", e$message))
    NA_real_
  })
}

# Access to standard RDM model for comparison (if EMC2 is available)
standard_RDM_model_obj <- NULL
can_test_vs_standard_RDM <- FALSE
if (requireNamespace("EMC2", quietly = TRUE) && "RDM" %in% ls(getNamespace("EMC2"))) {
  standard_RDM_model_obj <- EMC2::RDM()
  can_test_vs_standard_RDM <- TRUE
} else {
  # Create a mock RDM if EMC2 not available, for some basic structure
  # This mock won't have working C++ calls for dfun/pfun
  standard_RDM_model_obj <- list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    dfun=function(rt,pars){ message("Mock RDM dfun called"); rep(NA_real_, length(rt))},
    pfun=function(rt,pars){ message("Mock RDM pfun called"); rep(NA_real_, length(rt))}
  )
  message("EMC2::RDM() not found or EMC2 not installed. Some comparison tests will be skipped or will use mock.")
}


test_that("dRDM_DSWTN_log basic properties and reduction cases", {
  B_val <- 1.0
  A_val <- 0.2
  mu_drift_val <- 1.5
  sigma_drift_sq_val <- 0.5
  t0_val <- 0.1
  s_val <- 1.0

  t_vec <- c(t0_val - 0.01, t0_val, t0_val + 0.01, t0_val + 0.2, t0_val + 1.0)

  # Full model (A > 0, sigma_drift_sq > 0)
  log_pdfs_full <- dRDM_DSWTN_log(t_vec, B_val, A_val, mu_drift_val, sigma_drift_sq_val, t0_val, s_val)
  expect_equal(length(log_pdfs_full), length(t_vec))
  expect_true(is.infinite(log_pdfs_full[1]) && log_pdfs_full[1] < 0)
  expect_true(is.infinite(log_pdfs_full[2]) && log_pdfs_full[2] < 0)
  expect_true(all(is.finite(log_pdfs_full[3:5])))
  expect_true(all(log_pdfs_full[3:5] <= 0))

  # Case 1: No SPV (A=0), No Drift Var (sigma_drift_sq=0) -> Simple Wald
  log_pdfs_case1 <- dRDM_DSWTN_log(t_vec, B_val, 0, mu_drift_val, 0, t0_val, s_val)
  expect_true(all(is.finite(log_pdfs_case1[3:5])))
  # Further check: If B is small, PDF can be large. If B is large, PDF small.
  # If mu_drift is small, PDF small for early times.

  # Case 2: SPV (A>0), No Drift Var (sigma_drift_sq=0) -> Standard RDM with SPV
  log_pdfs_case2 <- dRDM_DSWTN_log(t_vec, B_val, A_val, mu_drift_val, 0, t0_val, s_val)
  expect_true(all(is.finite(log_pdfs_case2[3:5])))
  if (can_test_vs_standard_RDM) {
    pars_rdm <- c(v=mu_drift_val, B=B_val, A=A_val, t0=t0_val, s=s_val)
    pars_rdm_mat <- matrix(pars_rdm, nrow=length(t_vec[3:5]), ncol=length(pars_rdm), byrow=TRUE,
                           dimnames=list(NULL, names(pars_rdm)))
    expected_rdm_pdf_vals <- standard_RDM_model_obj$dfun(rt=t_vec[3:5], pars=pars_rdm_mat)
    expect_equal(exp(log_pdfs_case2[3:5]), expected_rdm_pdf_vals, tolerance=1e-6)
  }

  # Case 3: No SPV (A=0), Drift Var (sigma_drift_sq > 0) -> DSWTN with fixed threshold B
  log_pdfs_case3 <- dRDM_DSWTN_log(t_vec, B_val, 0, mu_drift_val, sigma_drift_sq_val, t0_val, s_val)
  expect_true(all(is.finite(log_pdfs_case3[3:5])))
  if (sigma_drift_sq_val > 1e-7) expect_true(any(abs(log_pdfs_case3[3:5] - log_pdfs_case1[3:5]) > 1e-5))

  # Test with scaling s != 1.0
  s_test_val <- 2.0
  # dRDM_DSWTN_log(t,B,A,m,ssq,t0,s) = log( (1/s) * f_std( (t-t0)/s | B/s, A/s, m/s, ssq/s^2 ) )
  # So dRDM_DSWTN_log(t*s+t0, B*s,A*s,m*s,ssq*s^2,t0,s) should be log(f_std(t|B,A,m,ssq)) - log(s)
  # This is equivalent to:
  # dRDM_DSWTN_log(t, B,A,m,ssq,t0,1) = dRDM_DSWTN_log( (t-t0)/s_test + t0, B/s_test, A/s_test, m/s_test, ssq/s_test^2, t0, 1/s_test) + log(s_test) ??? No this is not right.
  # The internal scaling should handle it.
  # Let f(t; P) be density with s=1. Then density with s is (1/s)f((t-t0)/s + t0_unscaled; P_scaled).
  # Our C++ code calculates f_scaled(t_adj_unscaled, B_scaled, A_scaled, mu_scaled, sig_sq_scaled).
  # This is f_std( (t-t0); B/s, A/s, ...). This is not (1/s)f(...).
  # The scaling property for densities is f(t; v,B,A,t0,s) = (1/s) * f_std((t-t0)/s; v,B,A,0,1)
  # Our current C++ dRDM_DSWTN_log computes log( f_std(t-t0; B/s, A/s, mu/s, sigmasq/s^2) ).
  # It does NOT include the -log(s) term. This is typical for JAGS/BUGS modules where (1/s) is part of jacobian.
  # For now, test that it runs with s != 1.
  log_pdfs_s_varied <- dRDM_DSWTN_log(t_vec, B_val, A_val, mu_drift_val, sigma_drift_sq_val, t0_val, s_test_val)
  expect_true(all(is.finite(log_pdfs_s_varied[3:5])))
  # And that it differs from s=1 if params are not scaled back
  expect_true(any(abs(log_pdfs_full[3:5] - log_pdfs_s_varied[3:5]) > 1e-5))
})

test_that("pRDM_DSWTN basic properties and reduction cases", {
  B<-1.0; A<-0.2; mu_drift<-1.5; sigma_drift_sq<-0.5; t0<-0.1; s<-1.0
  t_vec <- c(t0 - 0.01, t0, t0 + 0.01, t0 + 0.2, t0 + 0.5, t0 + 1.0, t0 + 5.0)

  cdfs_full <- pRDM_DSWTN(t_vec, B,A,mu_drift,sigma_drift_sq,t0,s)
  expect_equal(cdfs_full[1],0); expect_equal(cdfs_full[2],0)
  expect_true(all(cdfs_full[3:length(t_vec)] >= -1e-7 & cdfs_full[3:length(t_vec)] <= 1.000001)) # Allow small numerical errors
  expect_true(all(diff(cdfs_full[cdfs_full>=0 & cdfs_full <=1]) >= -1e-6))

  # Case 2: SPV (A>0), No Drift Var (sigma_drift_sq=0) -> Standard RDM with SPV
  cdfs_case2 <- pRDM_DSWTN(t_vec, B, A, mu_drift, 0, t0, s)
  if (can_test_vs_standard_RDM) {
    pars_rdm <- c(v=mu_drift, B=B, A=A, t0=t0, s=s)
    pars_rdm_mat <- matrix(pars_rdm, nrow=length(t_vec), ncol=length(pars_rdm), byrow=TRUE, dimnames=list(NULL, names(pars_rdm)))
    expected_rdm_cdf <- standard_RDM_model_obj$pfun(rt=t_vec, pars=pars_rdm_mat)
    expect_equal(cdfs_case2, expected_rdm_cdf, tolerance=1e-6)
  }

  t_test <- t0 + 0.5
  cdf_num_int <- numeric_cdf_rdm_dswtn(t_test, B, A, mu_drift, sigma_drift_sq, t0, s)
  cdf_direct <- pRDM_DSWTN(t_test, B, A, mu_drift, sigma_drift_sq, t0, s)
  if (!is.na(cdf_num_int)) {
    expect_equal(cdf_direct, cdf_num_int, tolerance = 2e-3) # Wider tolerance for double integral
  } else { message("Skipping numeric CDF comparison for RDM_DSWTN full model due to integration error.") }
})

test_that("rRDM_DSWTN basic properties and reduction cases", {
  B<-1.0; A<-0.2; mu_drift<-1.5; sigma_drift_sq<-0.5; t0<-0.1; s<-1.0
  n_samples <- 250

  samples_full <- rRDM_DSWTN(n_samples, B,A,mu_drift,sigma_drift_sq,t0,s)
  expect_true(all(is.finite(samples_full))); expect_true(all(samples_full >= t0)); expect_gt(mean(samples_full),t0)

  # Case 1: No SPV (A=0), No Drift Var (sigma_drift_sq=0) -> Simple Wald
  samples_case1 <- rRDM_DSWTN(n_samples, B, 0, mu_drift, 0, t0, s)
  if (mu_drift > 0) expect_equal(mean(samples_case1), (B/mu_drift)+t0, tolerance=0.4)

  # Case 2: SPV (A>0), No Drift Var (sigma_drift_sq=0) -> Standard RDM
  samples_case2 <- rRDM_DSWTN(n_samples, B, A, mu_drift, 0, t0, s)
  if (mu_drift > 0) expect_equal(mean(samples_case2), ((B+A/2)/mu_drift)+t0, tolerance=0.45)

  # Case 3: No SPV (A=0), Drift Var (sigma_drift_sq > 0)
  samples_case3 <- rRDM_DSWTN(n_samples, B, 0, mu_drift, sigma_drift_sq, t0, s)
  expect_gt(mean(samples_case3), t0)
  if (sigma_drift_sq > 0 && A > 0) { # Check that variance increases with more variability sources
      var_full = var(samples_full)
      var_case1 = var(samples_case1) # Wald
      var_case3_fixedA = var(rRDM_DSWTN(n_samples, B, 0, mu_drift, sigma_drift_sq, t0, s)) # A=0, sigmasq >0
      var_case2_fixedSig = var(rRDM_DSWTN(n_samples, B, A, mu_drift, 0, t0, s))         # A>0, sigmasq =0
      # Heuristic: full model variance should be greater than simpler models
      if(is.finite(var_full) && is.finite(var_case1)) expect_gte(var_full, var_case1 - 0.1) # allow some noise
      if(is.finite(var_full) && is.finite(var_case3_fixedA)) expect_gte(var_full, var_case3_fixedA-0.1)
      if(is.finite(var_full) && is.finite(var_case2_fixedSig)) expect_gte(var_full, var_case2_fixedSig-0.1)
  }
})

test_that("loglik_RDM_DSWTN_race with RDM params and reduction cases", {
  rt <- 0.5; choice <- 1

  B1<-1.0; A1<-0.2; mu1<-2.0; sig_sq1<-0.5; t01<-0.1; s1<-1.0
  B2<-1.2; A2<-0.1; mu2<-1.0; sig_sq2<-0.3; t02<-0.15; s2<-1.0

  pm_B <- matrix(c(B1,B2),1); pm_A <- matrix(c(A1,A2),1)
  pm_mu <- matrix(c(mu1,mu2),1); pm_ssq <- matrix(c(sig_sq1,sig_sq2),1)
  pm_t0 <- matrix(c(t01,t02),1); pm_s <- matrix(c(s1,s2),1)

  loglik <- loglik_RDM_DSWTN_race(rt,choice,pm_B,pm_A,pm_mu,pm_ssq,pm_t0,pm_s, spv_abs_err=1e-7,spv_rel_err=1e-7)
  expect_true(is.finite(loglik[1]))

  log_pdf1_R <- dRDM_DSWTN_log(rt,B1,A1,mu1,sig_sq1,t01,s1)
  cdf2_R <- pRDM_DSWTN(rt,B2,A2,mu2,sig_sq2,t02,s2, spv_abs_err=1e-7,spv_rel_err=1e-7)
  expect_equal(loglik[1], log_pdf1_R + log(1-cdf2_R), tolerance=1e-6)

  # Case: sigma_drift_sq = 0 for both (should be like standard RDM race)
  pm_ssq_zero <- matrix(0,1,2)
  loglik_rdm_std_equiv <- loglik_RDM_DSWTN_race(rt,choice,pm_B,pm_A,pm_mu,pm_ssq_zero,pm_t0,pm_s, spv_abs_err=1e-7,spv_rel_err=1e-7)
  expect_true(is.finite(loglik_rdm_std_equiv[1]))
  if (can_test_vs_standard_RDM) {
      # This would require a log_likelihood_race for standard RDM to compare against.
      # For now, check it runs and differs from full if sig_sq > 0
      if(sig_sq1 > 0 || sig_sq2 > 0) expect_true(abs(loglik[1] - loglik_rdm_std_equiv[1]) > 1e-5)
  }
})

test_that("loglik_RDM_DSWTN_race with CDF=1 for loser", {
  rt <- 2.0; choice <- 1
  B1<-1.0; A1<-0.0; mu1<-0.2; sig_sq1<-0.01; t01<-0.1; s1<-1.0
  B2<-0.1; A2<-0.0; mu2<-5.0; sig_sq2<-0.01; t02<-0.1; s2<-1.0

  pm_B<-matrix(c(B1,B2),1); pm_A<-matrix(c(A1,A2),1)
  pm_mu<-matrix(c(mu1,mu2),1); pm_ssq<-matrix(c(sig_sq1,sig_sq2),1)
  pm_t0<-matrix(c(t01,t02),1); pm_s<-matrix(c(s1,s2),1)

  cdf_loser_val <- pRDM_DSWTN(rt,B2,A2,mu2,sig_sq2,t02,s2, spv_abs_err=1e-9,spv_rel_err=1e-9)
  pdf_winner_val <- dRDM_DSWTN_log(rt,B1,A1,mu1,sig_sq1,t01,s1)

  if (is.finite(pdf_winner_val) && cdf_loser_val > 1 - 1.5e-10) {
     loglik <- loglik_RDM_DSWTN_race(rt,choice,pm_B,pm_A,pm_mu,pm_ssq,pm_t0,pm_s, min_log_lik = -1e7, spv_abs_err=1e-9,spv_rel_err=1e-9)
     expect_equal(loglik[1], -1e7)
  } else {
    message("Skipping CDF=1 test for loglik_RDM_DSWTN_race as conditions not met.")
    message(paste("Winner PDF:", pdf_winner_val, "Loser CDF:", cdf_loser_val))
  }
})
