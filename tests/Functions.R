# ===== Core pieces =====
rho_tau <- function(tau, tau_p) {
  s  <- sqrt(2*tau + 1); sp <- sqrt(2*tau_p + 1)
  (s - sp) / (s + sp)
}
C_tau <- function(tau, tau_p, b) {
  r <- rho_tau(tau, tau_p)
  (1 - 2*(b*b)*r) * exp(-(b*b)*r)
}
H_tau <- function(tau, tau_p, b, nu_tau, nu_tau_p) {
  if (tau_p >= tau) 0 else C_tau(tau, tau_p, b) * (nu_tau_p - nu_tau)
}

# exact panel weights for linear product with (τ-s)^(-3/2)
panel_weights_32 <- function(tau, a, b) {
  ra <- sqrt(max(tau - a, 0)); rb <- sqrt(max(tau - b, 0)); h <- b - a
  if (h <= 0) return(c(Wa=0, Wb=0))
  if (abs(rb) < .Machine$double.eps) return(c(Wa = 2*ra/h, Wb = 0))
  Wb <- (2/h) * ((ra - rb)^2) / rb
  Wa <- 2*(1/rb - 1/ra) - Wb
  c(Wa=Wa, Wb=Wb)
}

# your product rule with endpoint H (no spline, no last-panel patch)
integrate_pdf_product_R <- function(tau, tau_grid, nu_vals, b_scaled) {
  N <- length(tau_grid) - 1L
  if (N < 1) return(list(I=0, panels=numeric(0)))
  nu_tau <- tail(nu_vals, 1)
  contrib <- numeric(N)
  for (i in seq_len(N)) {
    a <- tau_grid[i]
    b <- if (i == N) tau else tau_grid[i+1]
    h <- b - a; if (h <= 0) { contrib[i] <- 0; next }
    Ha <- H_tau(tau, a, b_scaled, nu_tau, nu_vals[i])
    Hb <- if (i == N) 0 else H_tau(tau, b, b_scaled, nu_tau, nu_vals[i+1])
    w  <- panel_weights_32(tau, a, b)
    contrib[i] <- w["Wa"]*Ha + w["Wb"]*Hb
  }
  list(I = sum(contrib), panels = contrib)
}

# ===== Reference integrator on each panel via u = sqrt(τ - τ′) =====
# Uses linear interpolation for ν on [a,b]. No spline dependency.
panel_ref_u <- function(tau, a, b, b_scaled, nu_at_a, nu_at_b, nu_tau, K = 400L) {
  # map τ′ ∈ [a,b] to u ∈ [√(τ-b), √(τ-a)]
  u0 <- sqrt(max(tau - b, 0))
  u1 <- sqrt(max(tau - a, 0))
  if (u1 <= u0) return(0)
  du <- (u1 - u0) / K
  
  # linear ν(τ′) on the panel
  nu_lin <- function(tp) {
    if (b == a) return(nu_at_a)
    w <- (tp - a) / (b - a)
    (1 - w)*nu_at_a + w*nu_at_b
  }
  
  f <- function(u) {
    if (u <= 0) return(0)                 # finite limit; contribution ~ 0 at a single point
    tp <- tau - u*u
    C  <- C_tau(tau, tp, b_scaled)
    2 * C * (nu_lin(tp) - nu_tau) / (u*u)  # Jacobian already included
  }
  
  # trapezoid on [u0,u1]
  acc <- 0.5*f(u0)
  for (k in 1:(K-1)) acc <- acc + f(u0 + k*du)
  acc <- acc + 0.5*f(u1)
  acc * du
}

integrate_pdf_ref_R <- function(tau, tau_grid, nu_vals, b_scaled, K_panel = 400L) {
  N <- length(tau_grid) - 1L
  if (N < 1) return(list(I=0, panels=numeric(0)))
  nu_tau <- tail(nu_vals, 1)
  contrib <- numeric(N)
  for (i in seq_len(N)) {
    a <- tau_grid[i]
    b <- if (i == N) tau else tau_grid[i+1]
    contrib[i] <- panel_ref_u(tau, a, b, b_scaled,
                              nu_at_a = nu_vals[i],
                              nu_at_b = if (i == N) nu_tau else nu_vals[i+1],
                              nu_tau = nu_tau,
                              K = K_panel)
  }
  list(I = sum(contrib), panels = contrib)
}

# ===== Compare per-panel (reference vs product) =====
compare_panels <- function(tau, tau_grid, nu_vals, b_scaled, K_panel = 400L) {
  prod <- integrate_pdf_product_R(tau, tau_grid, nu_vals, b_scaled)
  ref  <- integrate_pdf_ref_R(tau, tau_grid, nu_vals, b_scaled, K_panel)
  
  N <- length(tau_grid) - 1L
  a_vec <- tau_grid[1:N]
  b_vec <- tau_grid[2:(N+1)]
  # ensure the last panel ends at tau (some grids already do; this just enforces it)
  b_vec[N] <- tau
  
  stopifnot(length(prod$panels) == N, length(ref$panels) == N)
  
  data.frame(
    panel = seq_len(N),
    a = a_vec,
    b = b_vec,
    I_prod = prod$panels,
    I_ref  = ref$panels,
    diff   = ref$panels - prod$panels,
    cum_prod_L = cumsum(prod$panels),
    cum_ref_L  = cumsum(ref$panels),
    cum_prod_R = rev(cumsum(rev(prod$panels))),
    cum_ref_R  = rev(cumsum(rev(ref$panels)))
  )
}

# Canonical A, B from t*, b, z (scaled)
AB_from_tstar <- function(t_scaled, b_scaled, z_scaled) {
  et <- exp(t_scaled)
  e2t <- exp(2*t_scaled)
  A   <- - (et*b_scaled - z_scaled) *
    exp(- (et*b_scaled - z_scaled)^2 / (e2t - 1) + 2*t_scaled) /
    (sqrt(pi) * (e2t - 1)^(3/2))
  Bc  <- et*b_scaled + e2t / (sqrt(pi) * sqrt(e2t - 1))
  list(A=A, Bcoeff=Bc, e2t=e2t, et=et)
}

# Needed integral from the ground-truth pdf (blue), using the SAME ν(τ) as your integral:
I_needed <- function(g_gt, lambda, A, Bcoeff, nu_tau, e2t) {
  ((g_gt / lambda) - A + Bcoeff * nu_tau) * sqrt(8 * pi) / e2t
}


# Needed integral from the ground-truth pdf (blue), using the SAME ν(τ) as your integral:
pf_needed <- function(g_gt, lambda, A, Bcoeff, nu_tau, integral) {
  ((g_gt / lambda) - A + Bcoeff * nu_tau ) 
}

# Assemble prediction from YOUR numerical integral, and compare to g_gt:
predict_pdf_from_I <- function(I_num, lambda, A, Bcoeff, nu_tau, e2t) {
  lambda * (A - Bcoeff * nu_tau + (e2t / sqrt(8 * pi)) * I_num)
}

trapz <- function(x, y) {
  n <- length(x); if (n < 2) return(0)
  sum(diff(x) * (head(y, -1) + tail(y, 1)) * 0.5)
}

rmse <- function(yhat, y) sqrt(mean((yhat - y)^2))

analyze_pdf_components <- function(t, A, B1, B2, C, g_true) {
  stopifnot(length(t)==length(A), length(A)==length(B1),
            length(B1)==length(B2), length(B2)==length(C),
            length(C)==length(g_true))

  # Canonical assembly
  S0 <- A
  S1 <- A - B1
  S2 <- A - B1 - B2
  S3 <- A - B1 - B2 + C   # your g(t)

  # Mass checks (should end up ~1)
  masses <- c(
    mass_A   = trapz(t, S0),
    mass_A_B1= trapz(t, S1),
    mass_A_B = trapz(t, S2),
    mass_all = trapz(t, S3),
    mass_C   = trapz(t, C)
  )

  # Fit diagnostics: how far are we from g_true at each stage?
  fits <- c(
    rmse_A    = rmse(S0, g_true),
    rmse_A_B1 = rmse(S1, g_true),
    rmse_A_B  = rmse(S2, g_true),
    rmse_all  = rmse(S3, g_true)
  )

  # Sign-combo sweep for B1/B2/C with A fixed +1
  combos <- expand.grid(sB1=c(-1,1), sB2=c(-1,1), sC=c(-1,1))
  sweep <- apply(combos, 1, function(s) {
    ghat <- A + s[1]*B1 + s[2]*B2 + s[3]*C
    c(rmse=rmse(ghat, g_true), mass=trapz(t, ghat),
      sB1=s[1], sB2=s[2], sC=s[3])
  })
  sweep <- as.data.frame(t(sweep))
  sweep <- sweep[order(sweep$rmse), ]

  # Constrained linear check (should return ~1s if the algebra is right)
  # g_true ≈ α*A + β*(-B1) + γ*(-B2) + δ*C, no intercept
  X <- cbind(A, -B1, -B2, C)
  coef <- as.numeric(solve(t(X)%*%X, t(X)%*%g_true))
  names(coef) <- c("alpha_A","beta_mB1","gamma_mB2","delta_C")

  list(masses=masses, fits=fits, best_signs=head(sweep, 6), coef=coef,
       S=list(A=S0, A_B1=S1, A_B=S2, all=S3))
}

# Inputs: vectors on the same t-grid
# A, B, C from your code; g_emp is empirical pdf (histogram/CDF' on t)
fit_prefactor_k <- function(t, A, B, C, g_emp) {
  # mass-matching k
  mass_AB  <- trapz(t, A - B)
  mass_C   <- trapz(t, C)
  k_mass   <- (1 - mass_AB) / mass_C

  # LS k with trapezoid weights
  w        <- c(diff(t), 0) + c(0, diff(t))   # trapezoid eqv. weights: Δt_left+Δt_right
  w        <- 0.5 * w
  y        <- g_emp - (A - B)
  num      <- sum(w * C * y)
  den      <- sum(w * C * C)
  k_ls     <- num / den

  list(k_mass = k_mass, k_ls = k_ls,
       masses = c(mAB = mass_AB, mC = mass_C))
}

# apply k and compute diagnostics
diagnose_prefactor <- function(t, A, B, C, g_emp, k) {
  g_hat <- (A - B) + k * C
  list(
    mass_hat = trapz(t, g_hat),
    rmse     = sqrt(mean( (g_hat - g_emp)^2 )),
    g_hat    = g_hat
  )
}

# b(t) = b_inf + (b0 - b_inf) * exp(-(t/tau)^p)
# - Vectorized over t
# - Stable for large t via log-power trick and clipping
decay_bound <- function(t, b0, b_inf = 0, tau, p = 1) {
  if (!is.numeric(t)) stop("t must be numeric (vector ok).")
  if (!(is.numeric(b0) && is.numeric(b_inf) && is.numeric(tau) && is.numeric(p))) {
    stop("b0, b_inf, tau, p must be numeric scalars.")
  }
  if (tau <= 0 || p <= 0) stop("tau and p must be > 0.")
  
  t <- pmax(t, 0)                 # guard negative times
  amp <- b0 - b_inf
  # y = p*(log t - log tau); handle t=0 separately; clip to avoid overflow
  y <- ifelse(t == 0, -Inf, p * (log(t) - log(tau)))
  y <- pmin(y, 700)               # prevents exp(y) overflow (~exp(709) is near double limit)
  s <- ifelse(is.infinite(y), 0, exp(y))  # s = (t/tau)^p; s=0 when t==0
  exp_term <- exp(-s)
  b_inf + amp * exp_term
}

# Example:
# decay_bound(c(0, 0.1, 0.5, 1, 2, 5), b0 = 1.0, b_inf = 0.2, tau = 1.5, p = 1.2)
