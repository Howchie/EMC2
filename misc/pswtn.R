# You must have the 'mvtnorm' package installed.
# If not, run: install.packages("mvtnorm")
library(mvtnorm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Helper function: Classic Wald CDF (R version of your pwald_classic)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pwald_classic_R <- function(t, boundary_alpha, drift_rate_xi) {
  # Handle edge cases
  count=0;cdf_val=NULL
  for(xi in drift_rate_xi) {
    count=count+1
    if (t <= 0) return(0.0)
    if (boundary_alpha <= 0) return(ifelse(xi > 0, 1.0, 0.0))
    if (xi <= 1e-10) return(0.0)
    t_sqrt <- sqrt(t)
    term1_arg <- (xi * t - boundary_alpha) / t_sqrt
    term2_arg <- -(xi * t + boundary_alpha) / t_sqrt
  
    p1 <- pnorm(term1_arg)
    p2 <- 2.0 * boundary_alpha * xi + pnorm(term2_arg,log.p = TRUE)
  
    cdf_val[count] <- p1 + exp(p2)
  
    # Clamp values to [0, 1] for stability
    cdf_val[count] <- max(0, min(1, cdf_val[count]))
    
  }
  return(cdf_val)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Main Function: Analytic SWTN CDF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Analytic CDF for the Shifted Wald with Truncated Normal drift rate
#'
#' @param t_adj Time, already adjusted for non-decision time (t - t0)
#' @param alpha Boundary separation (threshold)
#' @param mu_drift Mean of the (untruncated) drift rate distribution
#' @param sigma_drift SD of the drift rate distribution
#' @return The cumulative probability, P(RT < t_adj)
pswtn_Genz <- function(t_adj, alpha, mu_drift, sigma_drift) {
  
  # --- 1. Handle edge cases ---
  if (t_adj <= 1e-10) return(0.0)
  if (alpha <= 1e-10) return(1.0)
  if (sigma_drift < 0) return(NaN)
  
  # --- 2. Handle special case: No drift variability -> standard Wald ---
  if (sigma_drift <= 1e-10) {
    return(pwald_classic_R(t = t_adj, boundary_alpha = alpha, drift_rate_xi = mu_drift))
  }
  
  # --- 3. Main Analytic Calculation (Corrected) ---
  
  # Convenience variables
  t <- t_adj
  mu <- mu_drift
  sigma <- sigma_drift
  sigma_sq <- sigma^2
  
  # Normalization constant for truncated normal drift rate: P(xi > 0)
  prob_xi_gt_0 <- pnorm(mu / sigma)
  if (prob_xi_gt_0 < 1e-100) return(0.0)
  
  # Pre-calculate terms for correlations and denominators
  denom <- sqrt(t * (1 + t * sigma_sq))
  
  # --- Term 1 Calculation ---
  # Corresponds to P(Z < (xit - alpha)/sqrt t, xi > 0)
  
  # Upper bounds for the standardized bivariate normal P(X<h, Y<k)
  h1 <- (mu * t - alpha) / denom
  k1 <- mu / sigma
  
  # Correlation for the first bivariate normal
  rho1 <- (sigma * t) / denom
  
  term1_prob <- pbvn_tvpack(h1, k1,rho1)
  
  # --- Term 2 Calculation ---
  # Corresponds to exp(...) * P(Z < -(xit+alpha)/sqrt t, xi > 0), where xi~N(mu',sigma^2)
  mu_prime <- mu + 2.0 * alpha * sigma_sq
  
  # Upper bounds for the second bivariate normal
  h2 <- (-mu_prime * t - alpha) / denom
  k2 <- mu_prime / sigma
  
  # Correlation for the second bivariate normal
  rho2 <- -rho1 # Correlation has opposite sign
  
  term2_prob = pbvn_tvpack(h2, k2,rho2)
  
  
  # --- Combine terms ---
  exp_term <- exp(2.0 * alpha * mu + 2.0 * alpha * alpha * sigma_sq)
  
  cdf_val <- (term1_prob + exp_term * term2_prob) / prob_xi_gt_0
  
  # Final checks for numerical stability
  cdf_val <- max(0, min(1, cdf_val))
  
  return(cdf_val)
}
pswtn_fast <- function(t_adj, alpha, mu_drift, sigma_drift) {
  
  # --- 1. Handle edge cases ---
  if (t_adj <= 1e-10) return(0.0)
  if (alpha <= 1e-10) return(1.0)
  if (sigma_drift < 0) return(NaN)
  
  # --- 2. Handle special case: No drift variability -> standard Wald ---
  if (sigma_drift <= 1e-10) {
    return(pwald_classic_R(t = t_adj, boundary_alpha = alpha, drift_rate_xi = mu_drift))
  }
  
  # --- 3. Main Analytic Calculation (Corrected) ---
  
  # Convenience variables
  t <- t_adj
  mu <- mu_drift
  sigma <- sigma_drift
  sigma_sq <- sigma^2
  
  # Normalization constant for truncated normal drift rate: P(xi > 0)
  prob_xi_gt_0 <- pnorm(mu / sigma)
  if (prob_xi_gt_0 < 1e-100) return(0.0)
  
  # Pre-calculate terms for correlations and denominators
  denom <- sqrt(t * (1 + t * sigma_sq))
  
  # --- Term 1 Calculation ---
  # Corresponds to P(Z < (xit - alpha)/sqrt t, xi > 0)
  
  # Upper bounds for the standardized bivariate normal P(X<h, Y<k)
  h1 <- (mu * t - alpha) / denom
  k1 <- mu / sigma
  
  # Correlation for the first bivariate normal
  rho1 <- (sigma * t) / denom
  term1_prob <- pbvn_drezner(h1, k1,rho1)
  
  # --- Term 2 Calculation ---
  # Corresponds to exp(...) * P(Z < -(xit+alpha)/sqrt t, xi > 0), where xi~N(mu',sigma^2)
  mu_prime <- mu + 2.0 * alpha * sigma_sq
  
  # Upper bounds for the second bivariate normal
  h2 <- (-mu_prime * t - alpha) / denom
  k2 <- mu_prime / sigma
  
  # Correlation for the second bivariate normal
  rho2 <- -rho1 # Correlation has opposite sign
  
  term2_prob = pbvn_drezner(h2, k2,rho2)
  
  
  # --- Combine terms ---
  exp_term <- exp(2.0 * alpha * mu + 2.0 * alpha * alpha * sigma_sq)
  
  cdf_val <- (term1_prob + exp_term * term2_prob) / prob_xi_gt_0
  
  # Final checks for numerical stability
  cdf_val <- max(0, min(1, cdf_val))
  
  return(cdf_val)
}
swt_cdf_GH <- function(t, a, mu, sigma, n = 40) {
  gh   <- gauss.quad.prob(n, "normal")          # N(0,1) nodes/weights
  x    <- mu + sigma * gh$nodes                 # drift grid
  keep <- x > 0                                 # **truncate here**
  if (!any(keep)) return(0)                     # all mass left of zero
  
  fw   <- wald_cdf(t, a, x[keep])               # Wald CDF at positive drifts
  w    <- gh$weights[keep]
  ppos <- pnorm(mu / sigma)                     # P(xi > 0)
  
  sum(w * fw) / ppos
}
pars <- expand.grid(mu = c( 1, 2, 4 ),
                    cv = c( 0.05,0.1,0.2,0.4, 0.6, 0.8, 0.99 ),
                    alpha = c( 0.5, 1, 1.5, 3 ),
                    t = c( .1, .3, .6, 1, 2 ))

f.ref  <- function(t, a, mu, sigma) pswtn_numeric_integral(t,a,mu,sigma)  # every time
f.fast <- function(t, a, mu, sigma) pswtn_Genz(t,a,mu,sigma) # with your switch

err <- mapply(function(t, alpha, mu, cv) {
  sigma <- cv * mu          #sigma = CV * mu
  f.fast(t, alpha, mu, sigma) -
    f.ref (t, alpha, mu, sigma)
},
pars$t, pars$alpha, pars$mu, pars$cv)
range(err)         # look for |err| > 1e-8
pars$err=err

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. For Verification: A Numerical Integration Version
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pswtn_numeric_R <- function(t_adj, alpha, mu_drift, sigma_drift) {
  
  # PDF of a N(mu, sigma) distribution truncated at 0
  dtnorm <- function(x, mean, sd) {
    norm_const <- pnorm(mean / sd) # P(X > 0)
    # Return 0 if x is non-positive or if no density above 0
    ifelse(x <= 0 | norm_const < 1e-100, 0, dnorm(x, mean, sd) / norm_const)
  }
  
  # The integrand is P(t|xi) * f(xi)
  integrand <- function(xi, t_val, alpha_val, mu_val, sigma_val) {
    pwald_classic_R(t = t_val, boundary_alpha = alpha_val, drift_rate_xi = xi) *
      dtnorm(x = xi, mean = mu_val, sd = sigma_val)
  }
  
  # Use R's built-in numerical integrator
  result <- integrate(
    integrand,
    lower = 0,
    upper = Inf,
    t_val = t_adj,
    alpha_val = alpha,
    mu_val = mu_drift,
    sigma_val = sigma_drift
  )
  
  return(result$value)
}

cdf_swtn <- function(t_adj, alpha, mu_xi, sigma_xi){
  s      <- sqrt(t_adj*(1 + t_adj*sigma_xi^2))
  rho    <-  (sigma_xi*t_adj)/s
  
  # first phi2
  h1 <-  (mu*t_adj - alpha)/s
  k1 <-  mu / sigma_xi
  term1 <- pmvnorm(upper=c(h1,k1),
                   corr=matrix(c(1,rho,rho,1),2))[1]
  
  # second phi2 (reflected drift)
  mu_p <- mu + 2*alpha*sigma_xi^2
  h2   <- (-mu_p*t_adj - alpha)/s
  k2   <-  mu_p / sigma_xi
  term2 <- pmvnorm(upper=c(h2,k2),
                   corr=matrix(c(1,-rho,-rho,1),2))[1]
  
  prob  <- (term1 + exp(2*alpha*mu + 2*alpha^2*sigma_xi^2)*term2) /
    pnorm(mu/sigma_xi)
  max(0,min(1,prob))
}


# numerical integration test
integrate(pdf_swtn, lower=theta, upper=0.5, alpha=alpha, theta=theta, mu_xi=mu_xi, sigma_xi=sigma_xi)
# analytic solution test
cdf_swtn(0.5, alpha, theta, mu_xi, sigma_xi)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. TESTING: Compare Analytic vs. Numeric
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Set Test Parameters ---
t <- 1.5      # Reaction time
alpha <- 0.8  # Boundary
mu <- 0.5     # Mean drift
sigma <- 0.3  # Drift variability (SD)
t0 <- 0.2     # Non-decision time

# Adjust time for t0
t_adj <- t - t0

# --- Calculate CDF using both methods ---
cat("--- Testing SWTN CDF ---\n")
cat(sprintf("Parameters: t=%.2f, alpha=%.2f, mu=%.2f, sigma=%.2f, t0=%.2f\n\n", t, alpha, mu, sigma, t0))

# Method 1: Analytic solution
cdf_analytic <- pswtn_analytic_R(t_adj, alpha, mu, sigma)
cat(sprintf("Analytic CDF Result:  \t%.10f\n", cdf_analytic))

# Method 2: Numerical integration (for verification)
cdf_numeric <- pswtn_numeric_R(t_adj, alpha, mu, sigma)
cat(sprintf("Numeric Integral Result:\t%.10f\n", cdf_numeric))
cat(sprintf("Difference:             \t%e\n\n", cdf_analytic - cdf_numeric))


# --- Test an edge case (low drift variability) ---
cat("--- Testing Edge Case (sigma -> 0) ---\n")
sigma_low <- 1e-8
cdf_analytic_low_sigma <- pswtn_analytic_R(t_adj, alpha, mu, sigma_low)
cdf_wald_classic <- pwald_classic_R(t_adj, alpha, mu)

cat(sprintf("Analytic (low sigma):   \t%.10f\n", cdf_analytic_low_sigma))
cat(sprintf("Standard Wald CDF:      \t%.10f\n", cdf_wald_classic))
cat(sprintf("Difference:             \t%e\n", cdf_analytic_low_sigma - cdf_wald_classic))


# Install necessary packages if you don't have them
# install.packages("mvtnorm")
# install.packages("statmod")
# install.packages("truncnorm")

# Load libraries
library(mvtnorm)
library(statmod)

#' Analytic CDF for the Drift-Mixture (SWTN+SPV) Model
#'
#' Calculates the cumulative distribution function (CDF) for a Wald process where
#' the drift rate is drawn from a truncated normal distribution and the start
#' point is drawn from a uniform distribution. This implementation uses the
#' analytic solution involving the bivariate normal CDF.
#'
#' @param t Numeric vector of response times.
#' @param alpha The base decision threshold (lower bound of the uniform SPV).
#' @param mu The mean of the (untruncated) normal drift rate distribution.
#' @param sigma The standard deviation of the normal drift rate distribution (`sv`).
#' @param A The width of the start-point variability (SPV) distribution. The
#'          threshold on any given trial is drawn from U(alpha, alpha + A).
#'          Defaults to 0 (no SPV).
#' @param n_gauss_nodes The number of nodes for Gauss-Legendre quadrature when
#'                      A > 0. 20 is typically very accurate.
#'
#' @return A numeric vector of CDF values F(t).
#' @examples
#' t_vals <- seq(0.1, 3, by = 0.05)
#'
#' # Case 1: No Start-Point Variability (A=0)
#' cdf_vals_no_spv <- prdmswtn_analytic(t_vals, alpha = 1.0, mu = 1.5, sigma = 0.4, A = 0)
#' plot(t_vals, cdf_vals_no_spv, type = 'l', main = "CDF (No SPV)", ylab = "F(t)")
#'
#' # Case 2: With Start-Point Variability (A > 0)
#' cdf_vals_spv <- prdmswtn_analytic(t_vals, alpha = 1.0, mu = 1.5, sigma = 0.4, A = 0.2)
#' lines(t_vals, cdf_vals_spv, col = "blue")
#' legend("bottomright", legend = c("A=0", "A=0.2"), col = c("black", "blue"), lty = 1)

prdmswtn_analytic_R <- function(t_adj, B, mu, sigma, A = 0.0, n_gauss_nodes = 20) {
  
  pswtn_R <- function(t_adj, alpha, mu, sigma){
    s      <- sqrt(t_adj*(1 + t_adj*sigma^2))
    rho    <-  (sigma*t_adj)/s
    
    # first phi2
    h1 <-  (mu*t_adj - alpha)/s
    k1 <-  mu / sigma
    term1 <- pmvnorm(upper=c(h1,k1),
                     corr=matrix(c(1,rho,rho,1),2))[1]
    
    # second phi2 (reflected drift)
    mu_p <- mu + 2*alpha*sigma^2
    h2   <- (-mu_p*t_adj - alpha)/s
    k2   <-  mu_p / sigma
    term2 <- pmvnorm(upper=c(h2,k2),
                     corr=matrix(c(1,-rho,-rho,1),2))[1]
    
    prob  <- (term1 + exp(2*alpha*mu + 2*alpha^2*sigma^2)*term2) /
      pnorm(mu/sigma)
    max(0,min(1,prob))
  }
  
  # --- Main logic ---
  if (A < 1e-7) {
    # Case 1: No start-point variability, just call the single-bound function
    return(pswtn_R(t_adj, B,mu,sigma))
  } else {
    # Case 2: Integrate over start-point variability using Gauss-Legendre
    
    # Get nodes and weights for [-1, 1] interval
    gl <- statmod::gauss.quad(n_gauss_nodes, "legendre")
    
    # Map nodes to the integration interval [alpha, alpha + A]
    k_nodes <- B + 0.5 * A * (gl$nodes + 1)
    
    # The integral is the sum of weights * f(nodes)
    # We need to compute this for each time point in 't'
    cdf_matrix <- vapply(k_nodes, function(k) pswtn_R(t_adj, k,mu,sigma), numeric(length(t)))
    
    # The integral is the weighted sum across the nodes (k) for each time (t)
    # The final scaling factor for Gauss-Legendre over [a,b] is (b-a)/2. Here, A/2.
    integral_val <- cdf_matrix %*% gl$weights * (A / 2)
    
    # The PDF of the uniform start point is 1/A.
    # The expected value E[CDF(t,k)] is integral CDF(t,k) * (1/A) dk
    return(as.vector(integral_val / A))
  }
}

#' Analytic PDF for the RDM-SWTN Model
#'
#' Calculates the probability density function (PDF) for a model incorporating
#' both RDM-style uniform start-point variability (SPV) and SWTN-style
#' truncated-normal drift-rate variability.
#'
#' @details
#' The function relies on the closed-form PDF for the SWTN model (drift-rate
#' variability only) derived by Steingroever et al. (2021). To account for
#' start-point variability, this function integrates the SWTN PDF over the
#' uniform distribution of start points, U(B, B+A), using Gauss-Legendre
#' quadrature for high accuracy and speed.
#'
#' @param t Numeric vector of response times. Must be positive.
#' @param B The base decision threshold (lower bound of the SPV uniform distribution).
#' @param mu The mean of the (untruncated) normal drift rate distribution.
#' @param A The width of the start-point variability (SPV) distribution. The
#'          threshold on any given trial is drawn from U(B, B+A).
#'          Defaults to 0 (no SPV).
#' @param sigma The standard deviation of the drift rate distribution (`sv`).
#' @param n_gauss_nodes The number of nodes for Gauss-Legendre quadrature.
#'                      Default of 20 is highly accurate for this problem.
#'
#' @return A numeric vector of probability density values f(t).

drdmswtn_analytic_R <- function(t, B, mu, A = 0.0, sigma = 0.0, n_gauss_nodes = 20) {
  
  # Helper function for the single-bound PDF (the closed-form SWTN formula)
  # This corresponds to dswtn in the C++ code and Eq. 3 in Steingroever et al. (2021)
  # This is vectorized over t.
  d_single_bound <- function(tt, alpha, mu, sigma) {
    tt[tt <= 0] <- NA # PDF is undefined for non-positive decision times
    
    # Handle the no-drift-variability case (reduces to standard Wald)
    if (sigma < 1e-7) {
      if (mu <= 0) return(rep(0, length(tt)))
      # Standard Wald PDF formula
      pdf_val <- alpha * exp(-((alpha - mu * tt)^2) / (2 * tt)) / sqrt(2 * pi * tt^3)
      pdf_val[is.nan(pdf_val)] <- 0
      return(pdf_val)
    }
    
    # --- Full SWTN PDF Formula ---
    v <- sigma^2 # The formula uses variance
    
    # Denominator of the normalization constant for the truncated normal drift
    # This is P(drift > 0), which is phi(mu/sigma)
    prob_xi_gt_0 <- pnorm(mu / sigma)
    # Avoid division by zero if drift is certainly negative
    prob_xi_gt_0[prob_xi_gt_0 < 1e-100] <- 1e-100 
    
    # Term 1: alpha
    term1 <- alpha
    
    # Term 2: sqrt part of the denominator
    term2 <- sqrt(2 * pi * tt^3 * (tt * v + 1))
    
    # Term 3: Normalization for truncated normal drift rate
    term3 <- 1 / prob_xi_gt_0
    
    # Term 4: The main exponential part
    exp_numerator <- -( (mu * tt - alpha)^2 )
    exp_denominator <- 2 * tt * (tt * v + 1)
    term4 <- exp(exp_numerator / exp_denominator)
    
    # Term 5: The normal CDF part
    pnorm_arg <- (alpha * v + mu) / sqrt(tt * v^2 + v)
    term5 <- pnorm(pnorm_arg)
    
    pdf_val <- (term1 / term2) * term3 * term4 * term5
    pdf_val[is.nan(pdf_val) | !is.finite(pdf_val) | pdf_val < 0] <- 0
    return(pdf_val)
  }
  
  # --- Main logic ---
  # Ensure t is a vector for vapply
  t_vec <- as.vector(t)
  
  if (A < 1e-7) {
    # Case 1: No start-point variability. Just call the single-bound function
    # with the fixed threshold B.
    return(d_single_bound(t_vec, alpha = B, mu = mu, sigma = sigma))
    
  } else {
    # Case 2: Integrate over start-point variability using Gauss-Legendre
    
    # Get nodes and weights for the [-1, 1] interval
    gl <- statmod::gauss.quad(n_gauss_nodes, "legendre")
    
    # Map nodes to the integration interval k ~ U(B, B+A)
    k_nodes <- B + 0.5 * A * (gl$nodes + 1)
    
    # Calculate the PDF at each time 't' for each start-point node 'k'
    pdf_matrix <- vapply(k_nodes,
                         function(k) d_single_bound(t_vec, alpha = k, mu = mu, sigma = sigma),
                         numeric(length(t_vec)))
    
    # The integral integral f(k) dk over [B, B+A] is approximated by:
    # (A/2) * sum(weights * f(k_nodes))
    # We want the *average* value, which is (1/A) * integral f(k) dk.
    # The 'A' terms cancel, leaving (1/2) * sum(...)
    integral_avg <- (pdf_matrix %*% gl$weights) * 0.5
    
    return(as.vector(integral_avg))
  }
}

# --- Example Usage and Sanity Check ---
# Check that the PDF integrates to ~1.0

# Parameters for the test
params <- list(B = 0.8, mu = 2.0, A = 0.3, sigma = 0.5)

# Create a fine grid of time points
t_grid <- seq(1e-4, 5, length.out = 5000)

# Calculate the PDF values
pdf_values <- drdmswtn_analytic(
  t = t_grid, B = params$B, mu = params$mu, A = params$A, sigma = params$sigma
)

# Use trapezoidal rule for numerical integration
# install.packages("pracma")
library(pracma)
area <- trapz(t_grid, pdf_values)

# Print the result
print(paste("Numerical integral of the PDF:", round(area, 5)))
# The result should be very close to 1.0, confirming correct normalization.

# Plot the PDF
plot(t_grid, pdf_values, type = 'l', lwd = 2,
     main = "RDM-SWTN Probability Density Function",
     xlab = "Response Time (s)", ylab = "Density f(t)",
     sub = paste0("Integral = ", round(area, 5)))


library(truncnorm) # For rtruncnorm

#' Random Number Generator for the Drift-Mixture (SWTN+SPV) Model
#'
#' Generates random response times from the specified model.
#'
#' @param n Number of samples to generate.
#' @param alpha, mu, sigma, A Model parameters, same as in prdmswtn_analytic.
#' @param t0 Non-decision time to be added to the generated RTs.
#' @return A numeric vector of n response times.
rrdmswtn <- function(n, alpha, mu, sigma, A = 0.0, t0 = 0.0) {
  
  # 1. Sample drift rates from N(mu, sigma^2) truncated at 0
  drifts <- rtruncnorm(n, a = 0, b = Inf, mean = mu, sd = sigma)
  
  # 2. Sample start points from U(alpha, alpha + A)
  if (A < 1e-7) {
    start_points <- rep(alpha, n)
  } else {
    start_points <- runif(n, min = alpha, max = alpha + A)
  }
  
  # 3. Sample from the Inverse Gaussian (Wald) distribution for each trial
  # The mean of the Wald is alpha/drift
  # The shape parameter 'lambda' is alpha^2. The 'statmod' package uses this
  # parameterization directly.
  # We need to handle cases where drift is near zero.
  drifts[drifts < 1e-9] <- 1e-9 # Avoid division by zero
  
  wald_mean <- start_points / drifts
  wald_shape <- start_points^2
  
  rt_samples <- statmod::rinvgauss(n, mean = wald_mean, shape = wald_shape)
  
  return(rt_samples + t0)
}

# --- Run the Simulation and Comparison Plot ---
set.seed(123)
params <- list(alpha = 1.0, mu = 1.5, sigma = 0.4, A = 0.2, t0 = 0.1)

# Generate a large number of random samples
simulated_rts <- rrdmswtn(
  n = 100000,
  alpha = params$alpha,
  mu = params$mu,
  sigma = params$sigma,
  A = params$A,
  t0 = params$t0
)

# Plot the empirical CDF from the simulation
plot(ecdf(simulated_rts),
     main = "Analytic CDF vs. Simulation",
     xlab = "Response Time (s)",
     ylab = "Cumulative Probability F(t)",
     xlim = c(0, 3.5),
     lwd = 2,
     col = "darkorange",
     pch = NA, # Don't plot points, just the step function line
     verticals = TRUE
)

# Overlay the analytic CDF
t_grid <- seq(0, 3.5, length.out = 200)
# We subtract t0 from the grid because the analytic function models the decision time part
analytic_cdf <- prdmswtn_analytic(
  t = t_grid - params$t0,
  alpha = params$alpha,
  mu = params$mu,
  sigma = params$sigma,
  A = params$A
)
lines(t_grid, analytic_cdf, col = "steelblue", lwd = 2.5)

legend("bottomright",
       legend = c("Empirical CDF (100k Simulations)", "Analytic CDF"),
       col = c("darkorange", "steelblue"),
       lwd = 3,
       bty = "n")
