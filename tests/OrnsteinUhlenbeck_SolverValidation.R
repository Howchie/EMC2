## This script will loop through many parameter settings for the
## Ornstein-Uhlenbeck first hitting time distribution, computing the analytic solution where b_scaled == 0,
## and comparing several numerical solutions (grid-based, chunked grid-based, spline-based) 
devtools::load_all()
library(logspline)
ou_debug_set(FALSE, 1)
iterations = 1000
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
    ecdf = rep(NA, iterations),
    cdf_analytic = rep(NA, iterations),
    cdf_single_volterra = rep(NA, iterations),
    cdf_chunked_max_t = rep(NA, iterations),
    pdf_logspline = rep(NA, iterations)
)
# Parameters
for (i in 1:iterations) {
    lambda = params$lambda[i]
    theta = params$theta[i]
    sigma = params$sigma[i]
    z0 = params$z0[i]
    b = theta # fixed to theta for analytic solution
    b_inf = b # no decay (fixed_b)
    pow=3
    tau_exp = 4
    t = round(params$t[i])
    sim_dat <- simulate_ou_hit_times_bb(1e6, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow,
        dt = 1e-2, t_max = 10, adapt_factor = 1e6)
    params$pdf_analytic[i] = ou_fht_pdf_vec_closed_form(tt[t], lambda, theta, sigma, z0, b, b_inf)
    # analytic, on-grid pdf (no chunking or splines)
    params$pdf_single_volterra[i] = ou_fht_pdf_vec(tt[t], lambda, theta, sigma, z0, b, b_inf, tau_exp, pow, 0.01, 0)
    # Solve kernel once at max_t and smartly chunk the time grid used for kernel estimation
    params$pdf_chunked_max_t[i] = ou_fht_pdf_vec_grid_chunked(tt, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow,
        0.01, 0, .01, 5, 100, 30)[t]
    #params$cdf_analytic[i] = ou_fht_cdf_vec_closed_form(tt[t], lambda, theta, sigma, z0, b, b_inf)
    params$cdf_chunked_max_t[i] = ou_fht_cdf_vec_grid_chunked(tt, lambda, theta, sigma, z0, b, b_inf, tau_exp, pow, 0.01, 0, .01, 5, 100, 30)[t]
    params$cdf_single_volterra[i] = ou_fht_cdf_vec(tt[t], lambda, theta, sigma, z0, b, b_inf, tau_exp, pow)
    ec <- ecdf(sim_dat)
    params$ecdf[i] = ec(tt[t])
    
    miss_rate = mean(is.na(sim_dat))
    valid_rts <- sim_dat[!is.na(sim_dat)]
    if (length(valid_rts) > 1) {
        fit <- logspline::logspline(valid_rts)
        params$pdf_logspline[i] <- logspline::dlogspline(tt[t], fit) * (1 - miss_rate)
    }
}
write.csv(params, "ou_solver_validation_results.csv", row.names = FALSE)
par(mfrow=c(2,2))
plot(params$pdf_analytic, params$pdf_single_volterra, xlab="Analytic", ylab="Single Volterra"); abline(0,1,col="red")
plot(params$pdf_analytic, params$pdf_chunked_max_t, xlab="Analytic", ylab="Chunked Max t"); abline(0,1,col="red")
plot(params$pdf_analytic, params$pdf_logspline, xlab="Analytic", ylab="Logspline from Sim"); abline(0,1,col="red")
