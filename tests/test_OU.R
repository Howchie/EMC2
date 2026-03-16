# Compile (once per session)
# devtools::clean_dll()
# system("R CMD INSTALL .")
# devtools::build()
# devtools::load_all(reset=TRUE)
devtools::load_all()

# Helper: windowed PDF estimator with error bars
mc_pdf_window <- function(sim, t0, h) {
  sim <- sim[is.finite(sim)]
  N <- length(sim)
  if (N == 0) return(list(f=NA_real_, se=NA_real_, k=0L, N=0L))
  k <- sum(abs(sim - t0) <= h)
  fhat <- k / (N * 2*h)
  p <- max(k / N, .Machine$double.eps)
  se <- sqrt(p*(1-p)/N) / (2*h)
  list(f=fhat, se=se, k=k, N=N)
}

# Helper: KS distance between empirical CDF and model CDF evaluated on grid
ks_against_model <- function(sim, tgrid, F_model) {
  sim <- sim[is.finite(sim)]
  if (!length(sim)) return(NA_real_)
  ec <- ecdf(sim)
  max(abs(ec(tgrid) - F_model))
}

# Parameters
tt = seq(0.1,10,b=.01)
par(mfrow=c(2,2))
for (i in 1:4) {
    lambda <- runif(1,0.1,1.5)
    theta  <- runif(1,0.5,3)
    sigma  <- runif(1,0.1,1)
    z0 <- 0 # no start point variability
    b = theta #runif(1, 1, 3)
    b_inf = b # no decay (fixed_b)
    pow=3
    tau_exp=4

    c = sqrt(lambda) / sigma
    z_scaled = c * (z0 - theta)
    b_scaled = c * (b - theta)

    print(paste0("z_scaled = ", round(z_scaled,3), "; b_scaled = ", round(b_scaled,3), "; theta = ", round(theta,3)))
    
    library(microbenchmark)
    bench = microbenchmark(
        grid = ou_fht_pdf_vec_grid(tt, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow, 0.001, 300),
        chunked = ou_fht_pdf_vec_grid_chunked(
            tt, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow,
            0.001, 300, 1.5, 500, 12
        ),
        times = 5
    )

    sim_dat <- simulate_ou_hit_times_bb(5e6, lambda, theta, sigma, z0, b, z0, tau_exp, pow,
                                        dt = 1e-2, t_max = 10, adapt_factor = 1e6)
    misses = mean(is.na(sim_dat))
    # check at t=5 seconds via windowed estimator
    t = round(runif(1, 1, length(tt) - 1))
    t0 <- tt[t]
    h <- 0.1
    w <- mc_pdf_window(sim_dat, t0, h)
    fhat <- w$f; se <- w$se
    if (b==theta) {
        pdf_analytic = ou_fht_pdf_vec_closed_form(tt[t], lambda, theta, sigma, b, b, b_inf)
        cdf_analytic = ou_fht_cdf_vec_closed_form(tt[t], lambda, theta, sigma, z0, b, b_inf)
    }
    # analytic, on-grid pdf (no chunking or splines)
    pdf_on_grid = ou_fht_pdf_vec(tt[t], lambda, theta, sigma, z0, b, b, tau_exp, pow, 0.01, 0)
    cdf_on_grid = ou_fht_cdf_vec(tt[t], lambda, theta, sigma, z0, b, b, tau_exp, pow, 0.01, 0)
    # Solve kernel once at max_t and use the same kernel for all t
    pdf_grid = ou_fht_pdf_vec_grid(tt, lambda, theta, sigma, z0, b, b, tau_exp, pow, 0.0001,0, 0.001)[t]
    # as above but also smartly chunk the time grid used for kernel estimation using a minimum level for t_min, then gemoetrically widening the grid width across chunks
    pdf_chunked = ou_fht_pdf_vec_grid_chunked(tt, lambda, theta, sigma, z0, b, b, tau_exp, pow,
        0.01,0,.001, 5,100,30)
    pdf_chunked2 = ou_fht_pdf_vec_grid_chunked(tt, 2.75, theta, sigma, z0, b, b, tau_exp, pow,
        0.01,0,.001, 5,100,30)
    if (b==theta) {
        print(paste0("Closed-form PDF at t=", round(t0, 2), "s: ", round(pdf_analytic,6)))
    }
    print(paste0(
        "PDF at t=", round(t0, 2), "s: volterra=", round(pdf_on_grid, 6),
        "; spline=", round(pdf_spline, 6),
        "; grid=", round(pdf_grid, 6),
        "; chunked=", round(pdf_chunked, 6),
        "; mc_window=", round(fhat, 6),
        " (±", round(3*se,6), ")"
    ))

    cdf_chunked = ou_fht_cdf_vec_grid_chunked(tt, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow, 0.001, 300, 2, 300, 12)[tt[t]]
    cdf_backwards = ou_fht_cdf_vec(tt[t], lambda, theta, sigma, z0, b, b_inf, tau_exp, pow, 0.001, 300)
    ec <- ecdf(sim_dat)
    print(paste0(
        "CDF at t=", round(t0, 2), "s: volterra=", round(cdf_backwards, 6),
        "; chunked=", round(cdf_chunked, 6),
        "; ecdf=", round(ec(t0), 6)
    ))
}



    # Plot pdf against data
    hist(sim_dat2, breaks=100, freq=FALSE, main="TV b, no start point variability")
    lines(tt, pdf_on_grid/(1-misses), lwd=2, col='red')
    lines(tt, pdf_chunked2 / (1 - misses), lwd = 2, col = "red", lty="dashed")
    abline(v=t0, col="gray70", lty=3)
    segments(x0=t0-h, y0=fhat, x1=t0+h, y1=fhat, col="black", lwd=2)
    segments(x0=t0, y0=fhat-3*se, x1=t0, y1=fhat+3*se, col="black", lwd=2)
    # Empirical CDF vs. model CDF
    ec <- ecdf(sim_dat)
    plot(tt, ec(tt), type="l", lwd=2, ylim=c(0,1), main="CDF")
    lines(tt, cdf/(1-misses), col="red")
    lines(tt, cdf_chunked/(1-misses), col="blue", lty="dashed")
    ks <- ks_against_model(sim_dat, tt, cdf)
    cat("KS distance:", ks, "\n")
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


## Test GBM
mu = 2
sigma = 1
z0 = 0.5
b = 3
tt = seq(0.01, 5, by = 0.02)
gbm = simulate_gbm_hit_times_bb(5e5, mu + 0.5 * sigma * sigma, sigma, z0 * exp(b), exp(b), exp(b), adaptive = TRUE)
bm = simulate_bm_hit_times_bb(5e5, mu, sigma, z0, b, b, adaptive = TRUE)
pdf_bm = numeric(length(tt))
pdf_gbm = numeric(length(tt))
for (i in seq_along(tt)) {
  pdf_bm[i] = drdmswtn(tt[i], b, mu, 0, 0)
  pdf_gbm[i] = drdmswtn(tt[i], b, mu - 0.5, 0, 0)
}

## We can do a GBM with uniform start in natural units, 
## If we take the spv as a proportion of the natural-scale boundary

# Given t, vprime, sigma, rho
u0 <- -log(rho)
f  <- 0
for (i in 1:m) {
  v   <- nodes[i]         # standard Laguerre nodes
  w   <- weights[i]       # standard Laguerre weights
  u   <- u0 + v
  muIG <- u / vprime
  lam  <- (u*u) / (sigma*sigma)
  f   <- f + w * dinvgauss(t, mean = muIG, shape = lam)  # log-sum-exp in practice
}
f <- exp(-u0) * f / rho

## Test SWTN
b=3; v=2; sv=0.5; s=.8; c=0.5
tmp = rSWTN(1e5, b/s, v/s, 0/s, sv/s)
misses = mean(tmp>5)
tmp = tmp[tmp<5]
hist(tmp, breaks = 100, freq = FALSE)

misses = mean(out>5)
out = out[out<5]
hist(out, breaks=100, freq=FALSE)

pdf = dSWTNspv(seq(0.01,5,by=0.02), v, b, 0, 0, sv, s, c)
lines(seq(0.01,5,by=0.02), pdf/(1-misses), col='blue', lty='dashed',lwd=2)

b=3; v=2; sv=0.5; s=.8
time=1.65
pdf = dGBMspv(seq(0.01,time,by=0.02), b, v, s, 0, 0)
cdf_numeric = sum(pracma::trapz(seq(0.01, time, by = 0.02), pdf))
cdf_analytic = pgbm(time, b, v, s, 0, 0)

pdf = dWALDspv(seq(0.01,time,by=0.02), b, v, s, 0, 0)
cdf_numeric = sum(pracma::trapz(seq(0.01, time, by = 0.02), pdf))
cdf_analytic = pwald(time, b, v, s, 0, 0)
