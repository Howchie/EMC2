# Compile (once per session)
# devtools::clean_dll()
# system("R CMD INSTALL .")
# devtools::build()
# devtools::load_all(reset=TRUE)
devtools::load_all()
# Parameters
lambda <- 1
theta  <- 0
sigma  <- 1
z0 <- 2
b=1
c = sqrt(lambda) / sigma
z_scaled = c * (z0 - theta)
b_scaled = c * (b - theta)

## First, replicate plots from the paper
cols = viridis::cividis(5)
bs = c(1, 0.5, 0.01, -0.5, -1)
for (i in seq(1, length(bs))) {
  test = ou_fht_pdf_forward_debug(2, lambda, theta, sigma, z0, bs[i], num_steps = 300)
  if (i == 1) {
#    par(pty="s")
    plot(test$nu_f_vals_theta_block ~ test$tau_vals_grid, type = "l", col=cols[i], ylab="nu_f", xlab="tau")
  } else {
    lines(test$nu_f_vals_theta_block ~ test$tau_vals_grid, type = "l", col=cols[i])
  }
}



t_vec = seq(0.01,3,length.out=1000)
sim_dat <- simulate_ou_hit_times_std(5e4, lambda, theta, sigma, z0, b, dt = 1e-3, t_max = 5)
misses = mean(is.na(sim_dat))
tmp=density(sim_dat[!is.na(sim_dat)]); f=approxfun(tmp)
df = data.frame(
  t = numeric(length(t_vec)), e2t_sqrtpi = numeric(length(t_vec)), et_sqrtpi = numeric(length(t_vec)), integral = numeric(length(t_vec)),
  et = numeric(length(t_vec)), e2t = numeric(length(t_vec)), scaling = numeric(length(t_vec)), ref = numeric(length(t_vec)), fit = numeric(length(t_vec)),
  fit_adj = numeric(length(t_vec)), A = numeric(length(t_vec)), B = numeric(length(t_vec)), nuf = numeric(length(t_vec)),
  b_scaled = numeric(length(t_vec)), z_scaled = numeric(length(t_vec)), lambda = numeric(length(t_vec)), theta = numeric(length(t_vec)), sigma = numeric(length(t_vec)),abel=numeric(length(t_vec)), prefactor_dif = numeric(length(t_vec))
)

for(i in 1:length(t_vec)) {

  test = ou_fht_pdf_forward_debug(t_vec[i], lambda, theta, sigma, z0, b, num_steps = 300)
  pdfs = c(
    "A_trap" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid_uniform, 1)) + test$prefactor * test$A_trap) * test$lambda,
    "B_trap" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid, 1)) + test$prefactor * test$B_trap) * test$lambda,
    "C_trap" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid_right, 1)) + test$prefactor * test$C_trap) * test$lambda,
    "Theta_grid" = (test$termA - (test$termB * tail(test$nu_f_vals_theta_grid, 1)) + test$prefactor_theta * test$Theta_grid) * test$lambda,
    "A_u" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid_uniform, 1)) + test$prefactor * test$A_u) * test$lambda,
    "B_u" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid, 1)) + test$prefactor * test$B_u) * test$lambda,
    "C_u" = (test$termA - (test$termB * tail(test$nu_f_vals_tau_grid_right, 1)) + test$prefactor * test$C_u) * test$lambda,
    "Theta_block" = (test$termA - (test$termB * tail(test$nu_f_vals_theta_block, 1)) + test$prefactor_theta * test$Theta_block) * test$lambda,
    "reference" = ou_fht_pdf_from_cdf_regularized(t_vec[i], lambda, theta, sigma, z0, b, num_steps = 300),
    "reference_dat" = f(t_vec[i]) * (1 - misses)
  )
  theta_grid <- test$theta_grid_block
  nu_vals <- test$nu_f_vals_theta_block
  theta_max = test$theta
  nu_theta   <- tail(nu_vals, 1)                     # consistent ν(τ)
  b_scaled <- test$b_scaled
  z_scaled <- test$z_scaled
  tau      <- test$tau
  e2t      <- test$e2t
  lambda   <- test$lambda
  g_gt <- pdfs["reference_dat"]
  t_scaled <- test$t_scaled
  termA = test$termA; termB = test$termB; prefactor = test$prefactor_theta
  integral = integrate_pdf_forward_theta_u(theta_max, theta_grid, nu_vals, b_scaled)
  # recompute A,B from raw scalars (ignore any termA/termB naming)
  AB  <- AB_from_tstar(t_scaled, b_scaled, z_scaled)
  A <- AB$A

  B  <- AB$Bcoeff*nu_theta
  e2t <- AB$e2t
  et <- AB$et
  B1 <- et * nu_theta
  B2 <- (e2t / sqrt(pi * (e2t - 1))) * nu_theta
  C <- integral
  df$scaling[i] = ((g_gt / lambda) - (A - B)) / C
  df$ref[i] = g_gt
  df$et[i] = et
  df$e2t[i] = e2t
  df$t[i] = t_vec[i]
  df$e2t_sqrtpi[i] = e2t/sqrt(pi)
  df$et_sqrtpi[i] = et/sqrt(pi)
  df$integral[i] = integral
  df$fit[i] = pdfs["Theta_block"]
  df$fit_adj[i] = lambda * (A - B + C * df$scaling[i])
  df$A[i] = A
  df$B[i] = AB$Bcoeff
  df$nuf[i] = nu_theta
  df$abel[i] = abel_approx_nu_f(theta_max, z_scaled, b_scaled)
}
  df$b_scaled = b_scaled
  df$z_scaled = z_scaled
  df$lambda = lambda
  df$theta = theta
  df$sigma = sigma
  df$c = sqrt(lambda) / sigma
  df$prefactor_dif=df$scaling/df$e2t_sqrtpi
plot(df$ref ~ df$t, type="l", ylab="PDF", xlab="t",ylim=c(-1,2))
# lines(df$A ~ df$t, lty="dotted", col='blue')
# lines(df$B*df$nuf ~ df$t, lty="dashed", col='red')
# lines((df$A - df$B * df$nuf) ~ df$t, lty="dashed",col='blue')
lines(((df$A - df$B * df$nuf) + df$et_sqrtpi * df$integral)*df$lambda ~ df$t, col = "red", lty = 3)
#lines(df$integral ~ df$t, col = "green", lty = 4)
#lines(df$nuf ~ df$t, col = "purple", lty = 5)
tt  <- seq(0.001, 5, length.out = 300)
pdf1 <- ou_fht_pdf_from_cdf_regularized_vec(tt, lambda, theta, sigma, z0, b, num_steps = 300)
pdf2 <- ou_fht_pdf_forward_vec(tt, lambda, theta, sigma, z0, b, num_steps = 600)
cdf <- ou_fht_cdf_vec(tt, lambda, theta, sigma, z0, b, num_steps = 300)

# Sanity checks
plot(tt, pdf1, type="l", main="OU FHT PDF")
lines(tt, pdf2, type="l", main="OU FHT PDF Forward", col='blue')
plot(tt, cdf, type="l", main="OU FHT CDF")
max(cdf)           # <= 1
trap <- sum(0.5*(head(pdf1,-1)+tail(pdf1,-1))*diff(tt))
trap                 # ~ probability of ever hitting by max(tt)

# Simulate hitting times
set.seed(1)
sim_dat <- simulate_ou_hit_times_std(5e4, lambda, theta, sigma, z0, b, dt = 1e-3, t_max = 5)
misses=mean(is.na(sim_dat))
hist(sim_dat, breaks=100, freq=FALSE, main="Sim vs. Numerically-Derived PDF")
lines(tt, pdf1/(1-misses), lwd=2, col='blue')
lines(tt, pdf2/(1-misses), lwd=2, col='red')
# Empirical CDF vs. model CDF
ec <- ecdf(sim_dat)

plot(tt, cdf/(1-misses), type="l", lwd=2, ylim=c(0,1), main="Sim vs. Backwards CDF")
lines(tt, ec(tt), col=2)

set.seed(1)
lambda <- 1.2; theta <- 0; sigma <- 0.8; z0 <- -0.5; b <- 0.5
tt  <- seq(0.01, 5, length.out = 400)

cdf <- ou_fht_cdf_vec(tt, lambda, theta, sigma, z0, b, num_steps = 400)  # even
pdf <- ou_fht_pdf_vec(tt, lambda, theta, sigma, z0, b, num_steps = 400)

# numerical derivative of CDF
eps <- 1e-4
cdfp <- ou_fht_cdf_vec(tt+eps, lambda, theta, sigma, z0, b, num_steps = 400)
cdfm <- ou_fht_cdf_vec(pmax(0,tt-eps), lambda, theta, sigma, z0, b, num_steps = 400)
pdf_from_cdf <- (cdfp - cdfm) / (2*eps)

plot(tt, pdf_from_cdf, type="l")
max(abs(pdf - pdf_from_cdf), na.rm=TRUE)
