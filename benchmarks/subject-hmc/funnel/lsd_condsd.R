# Conditional sd of each log group sd given the other group coordinates, from the posterior draws.
e <- readRDS(commandArgs(TRUE)[1]); p <- e[[1]]$n_pars; pn <- e[[1]]$par_names
X <- do.call(rbind, lapply(e, function(ch) { w <- which(ch$samples$stage == "sample")
  t(sapply(w, function(t) { V <- as.matrix(ch$samples$theta_var[, , t]); R <- cov2cor(V)
    c(ch$samples$theta_mu[, t], 0.5 * log(diag(V)), atanh(R[upper.tri(R)])) })) }))
im <- 1:p; il <- p + 1:p; ic <- 2 * p + seq_len(p * (p - 1) / 2)
res <- function(y, Z) sd(lm.fit(cbind(1, Z), y)$residuals)
out <- t(sapply(1:p, function(i) { y <- X[, il[i]]
  c(marg = sd(y), given_mu = res(y, X[, im]), given_mu_lsd = res(y, X[, c(im, il[-i])]),
    given_all = res(y, X[, c(im, il[-i], ic)])) }))
rownames(out) <- pn; print(round(out, 3))
