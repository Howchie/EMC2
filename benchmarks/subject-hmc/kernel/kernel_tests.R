# Scratch correctness tests of the subject HMC kernel, outside the package.
source("R/subject_hmc.R")
set.seed(1)
cat("== 1. Gaussian force: sign, scale and finite-difference error ==\n")
p <- 5
A <- matrix(rnorm(p * p), p); SigL <- crossprod(A) / p + diag(0.3, p); PL <- solve(SigL)
m <- rnorm(p); mu <- rnorm(p); V <- diag(runif(p, 0.5, 2))
ll_fn <- function(X) { D <- sweep(X, 2, m); -0.5 * rowSums((D %*% PL) * D) }
fr <- .emc_hmc_frame(PL, mu, V)
z <- rnorm(p)
a <- fr$mu + drop(fr$S %*% z)
g_exact <- drop(crossprod(fr$S, -PL %*% (a - m) - solve(V, a - mu)))
for (h in c(1e-2, 1e-4, 1e-6, 1e-8)) {
  e <- .emc_hmc_eval(z, fr, ll_fn, h, -Inf)
  cat(sprintf("  h %.0e: max |g - exact| %.2e (|g| %.2f), rows %d\n", h,
              max(abs(e$g - g_exact)), sqrt(sum(g_exact^2)), e$rows))
}
# The whitened posterior is exactly N(z*, I): gradient = -(z - z*).
zstar <- drop(fr$R %*% (solve(PL + solve(V), PL %*% m + solve(V, mu)) - mu))
cat(sprintf("  whitened Gaussian check: max |g + (z - z*)| = %.2e\n",
            max(abs(.emc_hmc_eval(z, fr, ll_fn, 1e-4, -Inf)$g + (z - zstar)))))

cat("== 2. Energy error at fixed integration time ==\n")
banana_b <- 1.5
llb <- function(X) dnorm(X[, 2], banana_b * (X[, 1]^2 - 1), 0.4, log = TRUE)
mu2 <- c(0, 0); V2 <- diag(c(1, 4))
H2 <- matrix(c(4 * banana_b^2 / 0.4^2 * 0.2, 0, 0, 1 / 0.4^2), 2)  # rough metric
fr2 <- .emc_hmc_frame(H2, mu2, V2)
energies <- sapply(c(0.4, 0.2, 0.1, 0.05, 0.025), function(eps) {
  set.seed(3); e <- replicate(200, {
    a0 <- c(rnorm(1, 0, .7), rnorm(1)); l0 <- llb(matrix(a0, 1))
    tr <- .emc_hmc_transition(a0, l0, fr2, llb, eps / 1.0, 1e-5, -Inf)
    abs(tr$energy) })
  c(eps = eps, med = median(e, na.rm = TRUE), mean = mean(e, na.rm = TRUE))
})
print(round(t(energies), 5))
sl <- coef(lm(log(energies["med", ]) ~ log(energies["eps", ])))[2]
cat(sprintf("  log-log slope of median |energy error| vs eps: %.2f (expect ~2)\n", sl))

cat("== 3. Reversibility ==\n")
leap <- function(z0, v0, eps, L, fr, ll_fn, h) {
  cur <- .emc_hmc_eval(z0, fr, ll_fn, h, -Inf); v <- v0 + .5 * eps * cur$g; z <- z0
  for (l in 1:L) { z <- z + eps * v; new <- .emc_hmc_eval(z, fr, ll_fn, h, -Inf)
    v <- v + (if (l < L) eps else .5 * eps) * new$g }
  list(z = z, v = v)
}
z0 <- c(0.3, -1.2); v0 <- c(0.7, 0.4)
f <- leap(z0, v0, 0.13, 12, fr2, llb, 1e-5); b <- leap(f$z, -f$v, 0.13, 12, fr2, llb, 1e-5)
cat(sprintf("  forward-then-reverse: |z - z0| %.2e, |v + v0| %.2e\n", max(abs(b$z - z0)), max(abs(b$v + v0))))

cat("== 4. Invariance on the curved target from exact draws ==\n")
# Target: N(a; mu2, V2) * N(a2; b (a1^2 - 1), 0.4^2). Exact: a1 by rejection
# from N(0, 1) with acceptance N(b (a1^2 - 1); 0, 4 + 0.16) / max, then a2 | a1.
exact <- function(n) {
  out <- numeric(0)
  while (length(out) < n) { x <- rnorm(4 * n); k <- dnorm(banana_b * (x^2 - 1), 0, sqrt(4.16)) / dnorm(0, 0, sqrt(4.16))
    out <- c(out, x[runif(length(x)) < k]) }
  x <- out[1:n]; pm <- (banana_b * (x^2 - 1) / 0.16) / (1 / 4 + 1 / 0.16); ps <- sqrt(1 / (1 / 4 + 1 / 0.16))
  cbind(x, rnorm(n, pm, ps))
}
kernel_once <- function(X, eps, broken = FALSE, steps = 1) {
  t(apply(X, 1, function(a0) {
    for (k in seq_len(steps)) {
      l0 <- llb(matrix(a0, 1))
      if (!broken) { a0 <- .emc_hmc_transition(a0, l0, fr2, llb, eps, 1e-5, -Inf)$alpha } else {
        # Broken control: acceptance ignores the kinetic energy.
        U <- runif(1, .9, 1.1); v0 <- rnorm(2); u <- runif(1); e <- eps * U; L <- ceiling((pi/2)/e)
        cur <- .emc_hmc_eval(drop(fr2$R %*% (a0 - fr2$mu)), fr2, llb, 1e-5, -Inf, l0)
        v <- v0 + .5 * e * cur$g; z <- cur$z
        for (l in 1:L) { z <- z + e * v; new <- .emc_hmc_eval(z, fr2, llb, 1e-5, -Inf); v <- v + (if (l < L) e else .5 * e) * new$g }
        if (u < min(1, exp(new$lp - cur$lp))) a0 <- new$a
      }
    }
    a0 }))
}
stat <- function(X) c(m1 = mean(X[, 1]), m2 = mean(X[, 2]), v1 = mean(X[, 1]^2), q2 = mean(X[, 2] > 1), c12 = mean(X[, 1]^2 * X[, 2]))
check <- function(label, eps, broken, steps = 3, n = 6000) {
  set.seed(11); X0 <- exact(n); X1 <- kernel_once(X0, eps, broken, steps); R <- exact(n)
  ks1 <- suppressWarnings(ks.test(X1[, 1], R[, 1])$p.value); ks2 <- suppressWarnings(ks.test(X1[, 2], R[, 2])$p.value)
  # Paired comparison against the start draws: the kernel must not change the law.
  s0 <- stat(X0); s1 <- stat(X1)
  se <- apply(rbind(sapply(1:200, function(b) { i <- sample(n, n, TRUE); stat(X1[i, ]) - stat(X0[i, ]) })), 1, sd)
  zz <- (s1 - s0) / se
  cat(sprintf("  %-28s KS p (a1, a2) %.3f %.3f | max |z| moment shift %.2f\n", label, ks1, ks2, max(abs(zz))))
}
check("valid, eps 0.5", 0.5, FALSE)
check("valid, eps 1.2 (rough)", 1.2, FALSE)
check("BROKEN (no kinetic term)", 1.2, TRUE)

cat("== 5. Support boundary ==\n")
# Truncated target: a1 > 0 required (ll at the floor below). Exact draws by rejection.
fl <- -1e6
llt <- function(X) ifelse(X[, 1] > 0, dnorm(X[, 2], banana_b * (X[, 1]^2 - 1), 0.4, log = TRUE), fl - 1)
kernel_t <- function(X, eps) t(apply(X, 1, function(a0) { for (k in 1:3) { l0 <- llt(matrix(a0, 1))
  a0 <- .emc_hmc_transition(a0, l0, fr2, llt, eps, 1e-5, fl)$alpha }; a0 }))
set.seed(5); n <- 6000
ex_t <- function(n) { X <- exact(4 * n); X <- X[X[, 1] > 0, ]; X[1:n, ] }
X0 <- ex_t(n); X1 <- kernel_t(X0, 0.5); R <- ex_t(n)
cat(sprintf("  truncated: any a1 <= 0 after kernel: %s | KS p (a1, a2) %.3f %.3f\n", any(X1[, 1] <= 0),
            suppressWarnings(ks.test(X1[, 1], R[, 1])$p.value), suppressWarnings(ks.test(X1[, 2], R[, 2])$p.value)))
# Abort symmetry: an abort at a position is a function of the position, so a
# start whose probes cross the bound must abort without moving.
a_edge <- c(1e-6, 0); tr <- .emc_hmc_transition(a_edge, llt(matrix(a_edge, 1)), fr2, llt, 0.5, 1e-5, fl)
cat(sprintf("  start with a probe outside: status %s, moved %s, rows %d\n", tr$status, any(tr$alpha != a_edge), tr$rows))
cat("== 6. Allocator vs dense oracle ==\n")
source("benchmarks/subject-hmc/kernel/alloc_tests.R")
