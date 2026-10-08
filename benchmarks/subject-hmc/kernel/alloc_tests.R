oracle <- function(rhoP, rhoH, cost, C_fixed, rho_G, grid = seq(0, 1, by = 0.005)) {
  n <- nrow(rhoP); usable <- !apply(is.na(rhoH), 1, all)
  best <- Inf; arg <- NULL
  G <- if (any(usable)) as.matrix(expand.grid(rep(list(grid), sum(usable)))) else matrix(0, 1, 0)
  for (k in seq_len(nrow(G))) {
    r <- numeric(n); r[usable] <- G[k, ]
    rho <- max(rho_G, sapply(1:n, function(s) if (usable[s]) max(rhoP[s, ] * (1 - r[s] + r[s] * rhoH[s, ])) else max(rhoP[s, ])))
    obj <- .emc_hmc_tau(rho) * (C_fixed + sum(r * cost))
    if (obj < best - 1e-12) { best <- obj; arg <- r }
  }
  list(objective = best, rate = arg)
}
set.seed(7); worst <- 0
for (trial in 1:60) {
  n <- sample(1:2, 1); p <- sample(1:3, 1)
  rhoP <- matrix(runif(n * p, 0, 0.995), n); rhoH <- matrix(runif(n * p, 0, 0.99), n)
  if (runif(1) < .2) rhoH[1, ] <- NA
  if (runif(1) < .2) rhoP[n, 1] <- 0
  if (runif(1) < .2) rhoH[n, ] <- pmin(rhoH[n, ], 1)  # keep
  cost <- runif(n, 50, 500); Cf <- runif(1, 100, 2000); rG <- runif(1, 0, .9)
  a <- .emc_hmc_allocate(rhoP, rhoH, cost, Cf, rG)
  o <- oracle(rhoP, rhoH, cost, Cf, rG)
  gap <- if (is.infinite(o$objective) && is.infinite(a$objective)) 0 else (a$objective - o$objective) / o$objective
  worst <- max(worst, gap)
  if (gap > 1e-3) { cat("trial", trial, "solver", a$objective, "oracle", o$objective, "\n"); print(list(rhoP, rhoH, cost, Cf, rG, a$rate, o$rate)) }
}
cat(sprintf("  60 random problems: worst solver excess over grid oracle %.2e (negative = solver below grid)\n", worst))
# Named edge cases.
show <- function(lbl, ...) { a <- .emc_hmc_allocate(...); cat(sprintf("  %-34s rate %s tau %.3g unresolved %s\n", lbl, paste(signif(a$rate, 3), collapse = ","), a$tau, paste(a$unresolved, collapse = ","))) }
show("rhoP=0", matrix(0, 1, 2), matrix(.5, 1, 2), 100, 1000, 0)
show("rhoH=1 (no contraction)", matrix(.99, 1, 2), matrix(1, 1, 2), 100, 1000, 0)
show("rhoH=0 trapped subject", matrix(c(.999, .5), 1), matrix(0, 1, 2), 100, 1000, 0)
show("rhoH>=rhoP (0.2 vs 0.5)", matrix(.2, 1, 1), matrix(.5, 1, 1), 1, 1000, 0)
show("group-limited (rhoG .99)", matrix(.9, 2, 1), matrix(.1, 2, 1), c(100, 100), 1000, .99)
show("rhoP=1 unaddressable", matrix(c(1, .5), 2, 1), matrix(c(NA, .1), 2, 1), c(100, 100), 1000, 0)
show("rhoP=1 addressable", matrix(c(1, .5), 2, 1), matrix(c(.3, .1), 2, 1), c(100, 100), 1000, 0)
show("tie: identical subjects", matrix(.95, 2, 1), matrix(.2, 2, 1), c(100, 100), 1000, 0)
show("heterogeneous costs", matrix(c(.97, .97), 2, 1), matrix(c(.2, .2), 2, 1), c(10, 1000), 1000, 0)
show("all missing evidence (NA)", matrix(.97, 2, 1), matrix(NA_real_, 2, 1), c(10, 1000), 1000, 0)
