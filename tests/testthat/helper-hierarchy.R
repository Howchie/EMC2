# Exact draws from the diagonal hierarchy the group-move invariance tests
# use: y_s | a_s ~ N(a_s, S), a_s | mu, V ~ N(mu, V), mu ~ N(0, 1/P0 I),
# gamma(shape, rate) on each group precision, with the first parameter
# bounded below at b for every subject (the likelihood is zero beyond it),
# so a_1s is a truncated normal and the (mu_1, log sd_1) marginal carries
# the subjects' normalisers. (mu_k, log sd_k) come from a fine grid of their
# marginal posterior, then a_s from its conditional.
diag_hierarchy_draws <- function(N, b = -Inf, p = 2, n = 5, S = c(2, 1.5), m0 = 0,
                                 P0 = 0.5, shape = 2, rate = 0.3, seed = 21) {
  set.seed(seed)
  y <- matrix(rnorm(p * n, 0.5, 1.5), p)
  grid_draws <- function(k) {
    mu_g <- seq(-6, 6, length.out = 301); ls_g <- seq(-5, 2.5, length.out = 301)
    lp <- outer(mu_g, ls_g, Vectorize(function(m, l) {
      tau <- exp(-2 * l); v <- exp(2 * l)
      pv <- 1 / (1 / S[k] + 1 / v)
      trunc <- if (k == 1 && is.finite(b))
        sum(pnorm((pv * (y[k, ] / S[k] + m / v) - b) / sqrt(pv), log.p = TRUE)) else 0
      dnorm(m, m0, sqrt(1 / P0), log = TRUE) +
        sum(dnorm(y[k, ], m, sqrt(S[k] + v), log = TRUE)) + trunc +
        dgamma(tau, shape, rate, log = TRUE) + log(2 * tau)
    }))
    w <- exp(lp - max(lp))
    cell <- sample(length(w), N, TRUE, prob = w)
    hm <- diff(mu_g[1:2]); hl <- diff(ls_g[1:2])
    cbind(mu = mu_g[row(w)[cell]] + runif(N, -hm / 2, hm / 2),
          lsd = ls_g[col(w)[cell]] + runif(N, -hl / 2, hl / 2))
  }
  g <- lapply(seq_len(p), grid_draws)
  mu <- sapply(g, function(x) x[, "mu"]); lsd <- sapply(g, function(x) x[, "lsd"])
  a <- array(NA_real_, c(N, p, n))
  for (k in seq_len(p)) {
    v <- exp(2 * lsd[, k])
    pv <- 1 / (1 / S[k] + 1 / v)
    for (s in seq_len(n)) {
      pm <- pv * (y[k, s] / S[k] + mu[, k] / v)
      a[, k, s] <- if (k == 1 && is.finite(b)) {
        lo <- pnorm((b - pm) / sqrt(pv))
        pm + sqrt(pv) * qnorm(lo + runif(N) * (1 - lo))
      } else rnorm(N, pm, sqrt(pv))
    }
  }
  Si <- diag(1 / S)
  list(p = p, n = n, S = S, Si = Si, y = y, b = b, mu = mu, lsd = lsd, a = a,
       before = cbind(mu, lsd, matrix(a, N)),
       prior = list(theta_mu_mean = rep(m0, p), theta_mu_invar = diag(P0, p),
                    shape = rep(shape, p), rate = rep(rate, p)),
       # the likelihood as the mocked calc_ll_manager returns it
       ll = function(proposals, dadm) {
         dd <- sweep(proposals, 2, dadm)
         ifelse(proposals[, 1] >= b, -0.5 * rowSums((dd %*% Si) * dd), -1e4)
       },
       walls = list(U = matrix(c(1, numeric(p - 1)), p, n), b = rep(b, n),
                    Uperp = rep(list(diag(p)[, -1L, drop = FALSE]), n), has = rep(TRUE, n)))
}
