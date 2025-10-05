# Compile (once per session)
# devtools::clean_dll()
# system("R CMD INSTALL .")
# devtools::build()
# devtools::load_all(reset=TRUE)
devtools::load_all()
# Parameters
lambda <- 0.5
theta  <- 2.5
sigma  <- 1
z0 <- 0
b=2
pow=3
tau_exp=4

c = sqrt(lambda) / sigma
z_scaled = c * (z0 - theta)
b_scaled = c * (b - theta)

tt = seq(0.01,10,length.out=1000)
sim_dat <- simulate_ou_hit_times_std(5e4, lambda, theta, sigma, z0, b, z0, tau_exp,pow, dt = 1e-3, t_max = 10)
misses = mean(is.na(sim_dat))
tmp=density(sim_dat[!is.na(sim_dat)]); f=approxfun(tmp)

pdf <- ou_fht_pdf_forward_vec(tt, lambda, theta, sigma, z0, b, z0, tau_exp, pow, num_steps = 300)
tictoc::tic()
cdf <- ou_fht_cdf_vec(tt, lambda, theta, sigma, z0, b, z0, tau_exp, pow, num_steps = 300)
tictoc::toc()
# Sanity checks
plot(tt, pdf, type="l", main="OU FHT PDF")
lines(tt, f(tt) * (1 - misses), col = "purple", lwd = 2)
plot(tt, cdf, type="l", main="OU FHT CDF")
max(cdf)           # <= 1
trap <- sum(0.5*(head(pdf1,-1)+tail(pdf1,-1))*diff(tt))
trap                 # ~ probability of ever hitting by max(tt)

# Simulate hitting times
hist(sim_dat, breaks=100, freq=FALSE, main="TV b")
lines(tt, pdf/(1-misses), lwd=2, col='red')
# Empirical CDF vs. model CDF
ec <- ecdf(sim_dat)
plot(tt, cdf/(1-misses), type="l", lwd=2, ylim=c(0,1), main="Our version")
lines(tt, ec(tt), col=2)


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
