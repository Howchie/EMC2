# Surrogate group marginal log p(mu, V | y) (Laplace fits as exact likelihoods) along two scale coordinates of
# parameter i: (a) sd_i with correlations fixed (current shear), (b) residual sd of i given the others.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2"); for (f in ls(ns, all.names = TRUE)) if (grepl("^\\.emc_", f)) assign(f, get(f, ns))
a <- commandArgs(TRUE); e <- readRDS(a[1]); I <- as.integer(strsplit(a[2], ",")[[1]])
ch <- e[[1]]; p <- ch$n_pars; n <- ch$n_subjects; pn <- ch$par_names
fits <- .emc_shear_fits(ch); spec <- .emc_shear_v_spec(ch)
logmarg <- function(mu, V) {
  lp <- .emc_shear_log_prior_mu(mu, ch$prior) + .emc_shear_log_prior_v(V, spec)
  Vi <- chol2inv(chol(V)); ldV <- 2 * sum(log(diag(chol(V))))
  for (f in fits) { Q <- f$prec + Vi; RQ <- chol(Q); b <- f$lin + drop(Vi %*% mu)
    z <- forwardsolve(t(RQ), b)
    lp <- lp - 0.5 * ldV - sum(log(diag(RQ))) + 0.5 * sum(z^2) - 0.5 * sum(mu * (Vi %*% mu)) }
  lp }
tg <- seq(-2, 2, by = 0.05)
sdw <- function(lf) { lf[!is.finite(lf)] <- -Inf; wt <- exp(lf - max(lf)); wt <- wt / sum(wt); m <- sum(wt * tg); sqrt(sum(wt * (tg - m)^2)) }
# chain draws of log sd_i and log residual sd_i
w <- which(ch$samples$stage == "sample")
draws <- lapply(e, function(c) { ww <- which(c$samples$stage == "sample"); t(sapply(ww, function(t) { V <- as.matrix(c$samples$theta_var[, , t])
  c(0.5 * log(diag(V)), -0.5 * log(diag(solve(V)))) })) })
D <- do.call(rbind, draws)
lag1 <- function(k) { x <- lapply(draws, function(d) d[, k]); 1 - mean(unlist(lapply(x, function(z) diff(z)^2))) / (2 * var(unlist(x))) }
suppressPackageStartupMessages(library(posterior))
essk <- function(k) { m <- min(sapply(draws, nrow)); ess_bulk(sapply(draws, function(d) d[1:m, k])) }
for (i in I) {
  cat(sprintf("%s: marginal sd  log sd %.3f  log resid sd %.3f | lag1 rho %.2f / %.2f | ESS %.0f / %.0f\n", pn[i], sd(D[, i]), sd(D[, p + i]),
    lag1(i), lag1(p + i), essk(i), essk(p + i)))
  for (it in w[round(c(0.2, 0.5, 0.8) * length(w))]) {
    mu <- ch$samples$theta_mu[, it]; V <- as.matrix(ch$samples$theta_var[, , it])
    la <- vapply(tg, function(t) { d <- replace(rep(1, p), i, exp(t)); logmarg(mu, V * outer(d, d)) + (spec$size[i] + 1) * t }, 0)
    s2 <- 1 / solve(V)[i, i]
    lb <- vapply(tg, function(t) { V2 <- V; V2[i, i] <- V[i, i] + s2 * (exp(2 * t) - 1); logmarg(mu, V2) + 2 * t }, 0)
    cat(sprintf("   iter %4d: lsd %.2f, R2_i %.2f | orbit sd (a) sd_i,R fixed %.3f  (b) residual sd %.3f\n", it, 0.5 * log(V[i, i]), 1 - s2 / V[i, i], sdw(la), sdw(lb)))
  }
}
