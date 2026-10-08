source("benchmarks/subject-hmc/replay/setup11.R")
cat(sprintf("p %d, floor %.0f, ll at stuck state %.2f\n", p, floor_s, ll_stuck))
# Force stability across h at several states: stuck state, Laplace-conditional
# draws, and a state near the A lower region.
set.seed(2)
lap <- ns$.emc_laplace_conditional(structure(numeric(p), lik = ch$lik_laplace[[s]]), list(mu = mu, var = V))
states <- list(stuck = a_stuck, lapmean = lap$mu)
Ls <- t(chol(lap$var))
for (k in 1:3) states[[paste0("lapdraw", k)]] <- lap$mu + drop(Ls %*% rnorm(p))
far <- lap$mu + drop(Ls %*% rnorm(p)) * 2.5; states$far2.5sd <- far
hs <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
for (nm in names(states)) {
  a <- states[[nm]]; z <- drop(fr$R %*% (a - fr$mu))
  G <- sapply(hs, function(h) { e <- .emc_hmc_eval(z, fr, ll_fn, h, floor_s); if (e$ok) e$g else rep(NA, p) })
  ref <- G[, 3]
  rel <- apply(G, 2, function(g) sqrt(sum((g - ref)^2)) / sqrt(sum(ref^2)))
  adj <- sapply(1:(length(hs) - 1), function(k) sqrt(sum((G[, k] - G[, k + 1])^2)) / sqrt(sum(G[, k + 1]^2)))
  cat(sprintf("%-9s |g| %8.2f | rel diff vs h=1e-4: %s | adjacent: %s\n", nm, sqrt(sum(ref^2)),
      paste(sprintf("%.1e", rel), collapse = " "), paste(sprintf("%.1e", adj), collapse = " ")))
}
