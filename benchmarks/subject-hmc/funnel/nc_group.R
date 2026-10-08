# Full non-centred group update: given z_s = L^-1 (a_s - mu), propose (mu, L) jointly from the Gaussian
# conditional under the subjects' quadratic likelihood fits (prior on L by a quadratic expansion at the
# block-average V), accept by MH with the true likelihood. Reports acceptance and how far V moves.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2"); for (f in ls(ns, all.names = TRUE)) if (grepl("^\\.emc_", f)) assign(f, get(f, ns))
a <- commandArgs(TRUE); e <- readRDS(a[1]); K <- 40
ch <- e[[1]]; p <- ch$n_pars; n <- ch$n_subjects; pn <- ch$par_names
fits <- .emc_shear_fits(ch); spec <- .emc_shear_v_spec(ch)
w <- which(ch$samples$stage == "sample")
lt <- which(lower.tri(diag(p), diag = TRUE)); dl <- which(diag(p) %in% 1)  # positions
idx_l <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)
nL <- nrow(idx_l)
L_of <- function(v) { L <- matrix(0, p, p); L[idx_l] <- v; L }
logJ <- function(L) sum((p - seq_len(p) + 1) * log(diag(L)))   # |dV/dL| up to a constant
lpV <- function(L) { if (any(diag(L) <= 0)) return(-Inf); .emc_shear_log_prior_v(L %*% t(L), spec) + logJ(L) }
lpmu <- function(mu) .emc_shear_log_prior_mu(mu, ch$prior)
# posterior spread of lsd and correlations, for scale
post <- do.call(rbind, lapply(e, function(c) t(sapply(which(c$samples$stage == "sample"), function(t) {
  V <- as.matrix(c$samples$theta_var[, , t]); R <- cov2cor(V); c(0.5 * log(diag(V)), atanh(R[lower.tri(R)])) }))))
psd <- apply(post, 2, sd)
Vbar <- apply(ch$samples$theta_var[, , w], 1:2, mean); Lbar <- t(chol(Vbar))
# quadratic expansion of the L prior at Lbar (numerical)
f0 <- function(v) lpV(L_of(v)); v0 <- Lbar[idx_l]
gr <- numeric(nL); Hp <- matrix(0, nL, nL); h <- 1e-4
for (i in 1:nL) { ei <- replace(numeric(nL), i, h); gr[i] <- (f0(v0 + ei) - f0(v0 - ei)) / (2 * h) }
for (i in 1:nL) for (j in i:nL) { ei <- replace(numeric(nL), i, h); ej <- replace(numeric(nL), j, h)
  Hp[i, j] <- Hp[j, i] <- (f0(v0 + ei + ej) - f0(v0 + ei - ej) - f0(v0 - ei + ej) + f0(v0 - ei - ej)) / (4 * h^2) }
ev <- eigen(-Hp, symmetric = TRUE); Pp <- ev$vectors %*% (pmax(ev$values, 1e-3) * t(ev$vectors))
bp <- gr + Pp %*% v0      # linear term of the quadratic: -1/2 v'Pp v + bp'v
llf <- function(A) matrix(.emc_wpool_subject_ll(NULL, lapply(1:n, function(s) matrix(A[[s]], ncol = p, byrow = FALSE,
  dimnames = list(NULL, pn))), ch)$ll, ncol = n)
res <- NULL
for (it in w[round(seq(0.1, 1, length.out = 5) * length(w))]) {
  alpha <- ch$samples$alpha[, , it]; mu <- ch$samples$theta_mu[, it]; V <- ch$samples$theta_var[, , it]
  L <- t(chol(V)); z <- forwardsolve(L, alpha - mu)
  # a_s = mu + L z_s = A_s theta, theta = (mu, vech L)
  P <- matrix(0, p + nL, p + nL); b <- numeric(p + nL)
  for (s in 1:n) { Zs <- matrix(0, p, nL); for (k in 1:nL) Zs[idx_l[k, 1], k] <- z[idx_l[k, 2], s]
    A <- cbind(diag(p), Zs); P <- P + t(A) %*% fits[[s]]$prec %*% A; b <- b + drop(t(A) %*% fits[[s]]$lin) }
  P[1:p, 1:p] <- P[1:p, 1:p] + ch$prior$theta_mu_invar; b[1:p] <- b[1:p] + drop(ch$prior$theta_mu_invar %*% ch$prior$theta_mu_mean)
  P[p + 1:nL, p + 1:nL] <- P[p + 1:nL, p + 1:nL] + Pp; b[p + 1:nL] <- b[p + 1:nL] + bp
  R <- chol((P + t(P)) / 2); m <- backsolve(R, forwardsolve(t(R), b))
  lq <- function(th) -0.5 * sum((R %*% (th - m))^2)
  th0 <- c(mu, L[idx_l])
  props <- lapply(1:K, function(k) m + backsolve(R, rnorm(p + nL)))
  ok <- vapply(props, function(th) all(diag(L_of(th[p + 1:nL])) > 0), logical(1))
  # true target pieces
  amat <- function(th) { mu2 <- th[1:p]; L2 <- L_of(th[p + 1:nL]); mu2 + L2 %*% z }
  Acur <- amat(th0)
  Aall <- lapply(1:n, function(s) rbind(Acur[, s], t(sapply(props, function(th) amat(th)[, s]))))
  ll <- llf(lapply(Aall, function(M) t(M)))          # (K+1) x n
  tgt <- function(k, th) sum(ll[k, ]) + lpmu(th[1:p]) + lpV(L_of(th[p + 1:nL]))
  t0 <- tgt(1, th0)
  la <- vapply(1:K, function(k) if (!ok[k]) -Inf else (tgt(k + 1, props[[k]]) - t0) - (lq(props[[k]]) - lq(th0)), 0)
  acc <- mean(pmin(1, exp(la)))
  # how far: lsd and correlation change of proposals, in posterior sds
  feat <- function(th) { V2 <- L_of(th[p + 1:nL]) %*% t(L_of(th[p + 1:nL])); R2 <- cov2cor(V2); c(0.5 * log(diag(V2)), atanh(R2[lower.tri(R2)])) }
  f0v <- feat(th0); dist <- if (any(ok)) vapply(props[ok], function(th) sqrt(mean(((feat(th) - f0v) / psd)^2)), 0) else NA; cat("valid proposals", sum(ok), "of", K, "\n")
  cat(sprintf("iter %4d: acceptance %.3f | median log ratio %.1f | rms move of lsd+cor in posterior sds %.2f (median over proposals)\n",
    it, acc, median(la[is.finite(la)]), median(dist)))
}
