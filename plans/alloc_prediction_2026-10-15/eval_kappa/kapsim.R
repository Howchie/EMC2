KNOWN <- TRUE
set.seed(2)
ar <- function(n, a) { x <- numeric(n); x[1] <- rnorm(1); for (t in 2:n) x[t] <- a * x[t-1] + sqrt(1 - a^2) * rnorm(1); x }
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-12); (1 + r) / (1 - r) }
est <- function(X, M, mono) { n <- length(X[[1]]); pv <- if (KNOWN) 1 else var(unlist(X))
  r <- sapply(1:M, function(k) 1 - mean(sapply(X, function(x) mean((x[-(1:k)] - x[1:(n - k)])^2))) / (2 * pv))
  P <- r[seq(1, M, 2)] + r[seq(2, M, 2)]; s <- 1; prev <- Inf
  for (m in seq_along(P)) { pr <- P[m]; if (pr <= 0) break; if (mono) pr <- min(pr, prev); prev <- pr; s <- s + 2 * pr }
  c(k = max(1, s / tau1(r[1])), tau = s) }
cases <- list(ar0.5 = function() ar(100, 0.5), ar0.7 = function() ar(100, 0.7), ar0.9 = function() ar(100, 0.9),
  mix = function() sqrt(0.3) * ar(100, 0) + sqrt(0.7) * ar(100, 0.9), slow = function() sqrt(0.2) * ar(100, 0.3) + sqrt(0.8) * ar(100, 0.97))
truth <- c(ar0.5 = 1, ar0.7 = 1, ar0.9 = 1, mix = (0.3 + 0.7 * 19) / tau1(0.63), slow = (0.2 * tau1(.3) + 0.8 * tau1(.97)) / tau1(0.2 * .3 + 0.8 * .97))
for (nm in names(cases)) { R <- replicate(400, { X <- lapply(1:3, function(i) cases[[nm]]())
  c(est(X, 50, FALSE)[1], est(X, 50, TRUE)[1], est(X, 24, TRUE)[1], est(X, 12, TRUE)[1]) })
  cat(sprintf("%-6s true kappa %.2f | IPS50 %s | IMS50 %s | IMS24 %s | IMS12 %s\n", nm, truth[nm],
    paste(round(quantile(R[1, ], c(.1, .5, .9)), 2), collapse = "/"), paste(round(quantile(R[2, ], c(.1, .5, .9)), 2), collapse = "/"),
    paste(round(quantile(R[3, ], c(.1, .5, .9)), 2), collapse = "/"), paste(round(quantile(R[4, ], c(.1, .5, .9)), 2), collapse = "/"))) }
