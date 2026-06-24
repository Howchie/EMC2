# Tests that dexg's Mills-ratio tail branch (z < -8) agrees with R's pnorm
# reference implementation to within a tight tolerance.
#
# The tail branch activates when z = (x-mu)/sigma - sigma/tau < -8.
# It uses the Laplace continued fraction (mills_cf_denom) rather than a
# first-order asymptotic approximation, so relative error should be ~1e-12.

.dexg_r <- function(x, mu, sigma, tau, log_d = FALSE) {
  z <- (x - mu) / sigma - sigma / tau
  log_f <- -log(tau) + (mu - x) / tau + sigma^2 / (2 * tau^2) +
    pnorm(z, log.p = TRUE)
  if (log_d) log_f else exp(log_f)
}

# Parameter sets chosen so z = (x-mu)/sigma - sigma/tau is well below -8.
# Columns: mu, sigma, tau, x, expected_z (approx)
.tail_cases <- list(
  list(mu = 0.5,  sigma = 0.10, tau = 0.20, x = -0.30),  # z ≈ -8.5
  list(mu = 0.5,  sigma = 0.10, tau = 0.20, x = -0.50),  # z ≈ -10.5
  list(mu = 0.5,  sigma = 0.10, tau = 0.20, x = -1.00),  # z ≈ -15.5
  list(mu = 0.5,  sigma = 0.10, tau = 0.20, x = -2.00),  # z ≈ -25.5
  list(mu = 0.30, sigma = 0.05, tau = 0.10, x = -0.20),  # z ≈ -10.5
  list(mu = 0.30, sigma = 0.05, tau = 0.10, x = -0.50),  # z ≈ -16.5
  list(mu = 1.00, sigma = 0.20, tau = 0.05, x =  0.10),  # z ≈ -8.5 (large sigma/tau)
  list(mu = 1.00, sigma = 0.20, tau = 0.05, x = -0.20)   # z ≈ -12.5
)

test_that("dexg tail branch (Mills CF) matches R pnorm reference to 1e-8", {
  for (p in .tail_cases) {
    z <- (p$x - p$mu) / p$sigma - p$sigma / p$tau
    # Confirm this case actually hits the tail branch
    expect_lt(z, -8.0,
      label = sprintf("z check: mu=%.2f sig=%.2f tau=%.2f x=%.2f", p$mu, p$sigma, p$tau, p$x))

    cpp_val <- EMC2:::dexg_c(p$x, p$mu, p$sigma, p$tau, log_d = TRUE)
    r_val   <- .dexg_r(p$x, p$mu, p$sigma, p$tau, log_d = TRUE)

    expect_equal(cpp_val, r_val, tolerance = 1e-8,
      label = sprintf("dexg log-density: mu=%.2f sig=%.2f tau=%.2f x=%.2f (z=%.1f)",
                      p$mu, p$sigma, p$tau, p$x, z))
  }
})

test_that("dexg both branches agree with R pnorm reference at the z = -8 boundary", {
  # Tests that neither branch has a large error near the switch point.
  # (Comparing the two branches directly would just measure density change over
  # distance, not approximation error — so we compare each against R's pnorm.)
  mu <- 0.5; sigma <- 0.10; tau <- 0.20
  # x such that z = (x-mu)/sigma - sigma/tau = -8 exactly
  x_boundary <- mu + sigma * (sigma / tau - 8)
  x_inside  <- x_boundary - 1e-6   # z just below -8: tail branch
  x_outside <- x_boundary + 1e-6   # z just above -8: standard branch

  expect_equal(EMC2:::dexg_c(x_inside,  mu, sigma, tau, log_d = TRUE),
               .dexg_r(x_inside,  mu, sigma, tau, log_d = TRUE),
               tolerance = 1e-8, label = "tail branch near boundary")

  expect_equal(EMC2:::dexg_c(x_outside, mu, sigma, tau, log_d = TRUE),
               .dexg_r(x_outside, mu, sigma, tau, log_d = TRUE),
               tolerance = 1e-8, label = "standard branch near boundary")
})
