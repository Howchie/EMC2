library(EMC2)

# Test Wald sub-CDF with t0 > 0 and Erlang kill
# We compare the analytic t0-shifted solution with a very simple numerical check
# (PDF sum approximation or just checking for continuity/monotonicity)

test_wald_t0 <- function() {
  t <- 0.8
  mu <- 2.0
  b <- 1.0
  A <- 0.4
  sigma <- 1.0
  t0 <- 0.2
  k <- 1.5
  
  # 1. Monotonicity check
  ts <- seq(0, 1.5, length.out=20)
  # Signature: pwald(t, mu, b, A, sigma, t0, lambda_g, lambda_k, log_out, kill_shape, guess)
  p1 <- sapply(ts, function(ti) EMC2:::pwald(ti, mu, b, A, sigma, t0, 0.0, k, FALSE, 1, FALSE))
  p2 <- sapply(ts, function(ti) EMC2:::pwald(ti, mu, b, A, sigma, t0, 0.0, k, FALSE, 2, FALSE))
  
  if (any(is.na(p1)) || any(is.na(p2))) stop("NAs in Wald CDF")
  if (any(diff(p1) < -1e-10)) stop("Erlang-1 CDF not monotonic")
  if (any(diff(p2) < -1e-10)) stop("Erlang-2 CDF not monotonic")
  
  # 2. t < t0 check
  p_early <- EMC2:::pwald(t0 - 0.01, mu, b, A, sigma, t0, 0.0, k, FALSE, 1, FALSE)
  if (p_early != 0) stop("CDF should be 0 before t0 for non-guess path")
  
  # 3. Limit check (k=0)
  p_k0 <- EMC2:::pwald(t, mu, b, A, sigma, t0, 0.0, 0.0, FALSE, 1, FALSE)
  p_std <- EMC2:::pwald(t - t0, mu, b, A, sigma, 0.0, 0.0, 0.0, FALSE, 1, FALSE)
  if (abs(p_k0 - p_std) > 1e-10) stop("k=0 limit failed")
  
  message("Wald Erlang t0 tests passed!")
}

test_gbm_t0 <- function() {
  t <- 1.2
  mu <- 1.5
  b <- 2.0
  A <- 0.5
  sigma <- 0.5
  t0 <- 0.3
  k <- 1.0
  
  ts <- seq(0, 2, length.out=20)
  # Signature: pgbm(t, mu, b, A, sigma, t0, lambda_g, lambda_k, log_out, kill_shape, guess)
  p1 <- sapply(ts, function(ti) EMC2:::pgbm(ti, mu, b, A, sigma, t0, 0.0, k, FALSE, 1, FALSE))
  p2 <- sapply(ts, function(ti) EMC2:::pgbm(ti, mu, b, A, sigma, t0, 0.0, k, FALSE, 2, FALSE))
  
  if (any(is.na(p1)) || any(is.na(p2))) stop("NAs in GBM CDF")
  if (any(diff(p1) < -1e-10)) stop("GBM Erlang-1 CDF not monotonic")
  if (any(diff(p2) < -1e-10)) stop("GBM Erlang-2 CDF not monotonic")
  
  message("GBM Erlang t0 tests passed!")
}

test_lba_t0 <- function() {
  t <- 1.0
  v <- 2.5
  b <- 1.2
  A <- 0.5
  sv <- 0.3
  t0 <- 0.2
  k_leak <- 0.0 # non-leaky
  lambda <- 1.2
  
  ts <- seq(0, 2, length.out=20)
  # Signature: pkilledleakyba(t, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, log_out, kill_shape, guess)
  p1 <- sapply(ts, function(ti) EMC2:::pkilledleakyba(ti, v, b, A, sv, t0, k_leak, 0.0, lambda, TRUE, FALSE, 1, FALSE))
  p2 <- sapply(ts, function(ti) EMC2:::pkilledleakyba(ti, v, b, A, sv, t0, k_leak, 0.0, lambda, TRUE, FALSE, 2, FALSE))
  
  if (any(is.na(p1)) || any(is.na(p2))) stop("NAs in LBA CDF")
  # LBA quadrature can have slight numerical noise, allow 1e-7
  if (any(diff(p1) < -1e-7)) {
    print(p1)
    stop("LBA Erlang-1 CDF not monotonic")
  }
  if (any(diff(p2) < -1e-7)) stop("LBA Erlang-2 CDF not monotonic")
  
  message("LBA Erlang t0 tests passed!")
}

test_wald_t0()
test_gbm_t0()
test_lba_t0()
