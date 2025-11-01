## This script will loop through many parameter settings for the
## Ornstein-Uhlenbeck first hitting time distribution, computing the analytic solution where b_scaled == 0,
## and comparing several numerical solutions (grid-based, chunked grid-based, spline-based) 
devtools::load_all()
iterations = 10000
min_t = 0.1
max_t = 10
tt = seq(min_t,max_t,b=.01)
params = data.frame(
    t = runif(iterations, 1, length(tt)-1),
    lambda = runif(iterations, 0.01, 3),
    theta  = runif(iterations, 0.5, 3),
    sigma  = runif(iterations, 0.1, 1),
    z0     = runif(iterations,0,.49),
    pdf_analytic = rep(NA, iterations),
    pdf_single_volterra = rep(NA, iterations),
    pdf_chunked_max_t = rep(NA, iterations),
)
# Parameters
tt = seq(0.1,10,b=.01)
par(mfrow=c(2,2))
for (i in 1:iterations) {
    lambda params$lambda[i]
    theta params$theta[i]
    sigma params$sigma[i]
    z0 params$z0[i]
    b = theta # fixed to theta for analytic solution
    b_inf = b # no decay (fixed_b)
    pow=3
    tau_exp=4

    params$pdf_analytic[i] = ou_fht_pdf_vec_closed_form(tt[t], lambda, theta, sigma, z0, b, b_inf)
    # analytic, on-grid pdf (no chunking or splines)
    params$pdf_single_volterra[i] = ou_fht_pdf_vec(tt[t], lambda, theta, sigma, z0, b, z0, tau_exp, pow, 0.0005, 0)
    # Solve kernel once at max_t and smartly chunk the time grid used for kernel estimation using a minimum level for t_min, then gemoetrically widening the grid width across chunks
    pdf$pdf_chunked_max_t[i] = ou_fht_pdf_vec_grid_chunked(tt, lambda, theta, sigma, z0, b, b, tau_exp, pow,
        0.0005,0,.001, 2,200,15)[t]

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
    hist(sim_dat, breaks=1000, freq=FALSE, main="Fixed-b, no start point variability")
    lines(tt, pdf_on_grid/(1-misses), lwd=2, col='red')
    lines(tt, pdf_chunked / (1 - misses), lwd = 2, col = "green", lty="dashed")
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
