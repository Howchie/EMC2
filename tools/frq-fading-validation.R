library(EMC2)

conditional_density <- function(t, alpha, beta, h, kappa) {
  q_inf <- qbeta(h, alpha, beta)
  lambda <- -log1p(-q_inf) * kappa
  EMC2:::dfrqfade(t, alpha, beta, lambda, kappa, 0, 0) / h
}

conditional_kl <- function(alpha, beta, h1 = 0.95, h2 = 0.60,
                           free_second_shapes = FALSE) {
  f1 <- function(t) conditional_density(t, alpha, beta, h1, 1)
  kl <- function(theta) {
    if (free_second_shapes) {
      alpha2 <- 1 + exp(theta[1])
      beta2 <- 1 + exp(theta[2])
      kappa2 <- exp(theta[3])
    } else {
      alpha2 <- alpha
      beta2 <- beta
      kappa2 <- exp(theta[1])
    }
    f2 <- function(t) conditional_density(t, alpha2, beta2, h2, kappa2)
    integrate(function(t) {
      x <- f1(t)
      y <- f2(t)
      out <- numeric(length(t))
      keep <- x > 0 & y > 0
      out[keep] <- x[keep] * log(x[keep] / y[keep])
      out
    }, 0, 50, rel.tol = 1e-8)$value
  }

  start <- if (free_second_shapes) c(log(alpha - 1), log(beta - 1), log(1.8))
           else log(1.8)
  fit <- if (free_second_shapes) {
    optim(start, kl, method = "Nelder-Mead",
          control = list(maxit = 400, reltol = 1e-9))
  } else {
    optimize(kl, interval = c(-3, 3))
  }
  as.numeric(if (!is.null(fit$value)) fit$value else fit$objective)
}

kl_11 <- conditional_kl(1, 1)
kl_23 <- conditional_kl(2, 3)
kl_23_free <- conditional_kl(2, 3, free_second_shapes = TRUE)
stopifnot(abs(kl_11 - 0.0111) < 0.001,
          abs(kl_23 - 0.0097) < 0.001,
          abs(kl_23_free - 0.0077) < 0.001)
cat(sprintf("Conditional RT KL (h .95 to .60): (1,1) %.5f; (2,3) %.5f; (2,3), second shapes free %.5f\n",
            kl_11, kl_23, kl_23_free))

set.seed(20260923)
options(emc2.cpp_rfun = TRUE)
n <- 8000L
truth <- c(alpha = 1.4, beta = 2.2, lambda = 1.3)
kappa_truth <- c(0.55, 0.8, 1.1, 1.5)
sim <- lapply(kappa_truth, function(kappa) {
  pars <- cbind(alpha = truth["alpha"], beta = truth["beta"],
    lambda = truth["lambda"], kappa = kappa, t0 = 0.1, delta = 0, cv_u = 0)
  pars <- pars[rep(1, n), , drop = FALSE]
  EMC2:::.rfun_FRQfade(factor(rep("go", n), levels = "go"), pars)
})

loglik <- function(theta) {
  alpha <- 1 + exp(theta[1])
  beta <- 1 + exp(theta[2])
  lambda <- exp(theta[3])
  kappas <- exp(theta[4:7])
  total <- 0
  for (j in seq_along(sim)) {
    h <- EMC2:::frq_fade_summary(alpha, beta, lambda, kappas[j], 0, 0)[1, "h"]
    u <- sim[[j]]$rt[is.finite(sim[[j]]$rt)] - 0.1
    dens <- EMC2:::dfrqfade(u, alpha, beta, lambda, kappas[j], 0, 0)
    if (!is.finite(h) || h <= 0 || h >= 1 || any(!is.finite(dens)) ||
        any(dens <= 0)) return(1e100)
    total <- total + sum(log(dens)) + (n - length(u)) * log1p(-h)
  }
  if (is.finite(total)) -total else 1e100
}

start <- c(log(1.35 - 1), log(2.5 - 1), log(1.2),
           log(0.6), log(0.85), log(1), log(1.4))
fit <- optim(start, loglik, method = "BFGS",
             control = list(maxit = 400, reltol = 1e-10))
recovered <- c(alpha = 1 + exp(fit$par[1]), beta = 1 + exp(fit$par[2]),
  lambda = exp(fit$par[3]), kappa = exp(fit$par[4:7]))
h_truth <- vapply(kappa_truth, function(kappa)
  EMC2:::frq_fade_summary(truth["alpha"], truth["beta"], truth["lambda"],
                          kappa, 0, 0)[1, "h"], numeric(1))
completion <- vapply(sim, function(x) mean(is.finite(x$rt)), numeric(1))
median_rt <- vapply(sim, function(x) median(x$rt[is.finite(x$rt)]), numeric(1))
stopifnot(fit$convergence == 0,
  abs(recovered["alpha"] / truth["alpha"] - 1) < 0.15,
  abs(recovered["beta"] / truth["beta"] - 1) < 0.25,
  abs(recovered["lambda"] / truth["lambda"] - 1) < 0.25,
  all(abs(recovered[4:7] / kappa_truth - 1) < 0.10),
  max(abs(completion - h_truth)) < 0.012,
  all(diff(completion) < 0), median_rt[1] - median_rt[4] > 0.02)

cat("Multi-condition single-accumulator recovery (8,000 trials/condition):\n")
print(data.frame(kappa = kappa_truth, h = as.numeric(h_truth),
  completed = completion, median_rt = median_rt))
print(recovered)
