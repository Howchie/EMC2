# Offline emulation of skew-normal likelihood frames for shear transport.
# Surrogate per subject: lin'a - a'Ha/2 + log Phi(c + lambda'a). The latent
# w | a ~ N(c + lambda'a, 1) truncated to w > 0 makes (a, w) jointly Gaussian
# under any group state; the map holds the truncated-normal CDF of w and the
# whitened a | w fixed. Walls are not applied (every method alike).
# Usage: Rscript skew_frames.R fixture.rds output.rds coord [chains] [draws] [reps]
#   coord: a parameter name (its log group sd) or "mean:<k>" (k-th mean direction)
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
coord <- args[3]
chains <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else seq_along(emc)
ch1 <- emc[[1]]; p <- ch1$n_pars; n <- ch1$n_subjects; pn <- ch1$par_names
last <- ch1$samples$idx
draws <- if (length(args) >= 5 && nzchar(args[5])) as.integer(strsplit(args[5], ",")[[1]]) else
  round(last - c(401, 301, 201, 101, 51, 21, 1, 0))
REPS <- if (length(args) >= 6) as.integer(args[6]) else 16L
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample")
k <- if (startsWith(coord, "mean:")) as.integer(sub("mean:", "", coord)) else p + match(coord, pn)
stopifnot(!is.na(k))
train_idx <- setdiff(idx, unlist(lapply(draws, function(d) (d - 5):(d + 5))))
grid <- seq(-2, 2, by = 0.5); g0 <- which(grid == 0)
spec <- ns$.emc_shear_v_spec(ch1)
ctx <- ns$.emc_shear_context(ch1)
lap <- ctx$packed
ll_rows <- function(X, s, ch = ch1) ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, pn)
set.seed(sum(utf8ToInt(paste("skew-frames", coord))))

# ---- Training points: window draws (stored likelihoods) and points the
# Laplace transport reaches along the coordinate from probe draws.
cpu_fit <- proc.time()
probe <- unlist(lapply(chains, function(j) sample(train_idx, 2)))
probe_chain <- rep(chains, each = 2)
path_pts <- vector("list", n)
path_rows <- 0
for (i in seq_along(probe)) {
  c <- emc[[probe_chain[i]]]; d <- probe[i]
  mu <- c$samples$theta_mu[, d]; V <- c$samples$theta_var[, , d]
  aa <- matrix(c$samples$alpha[, , d], p)
  path <- ns$.emc_shear_coord_path(c, idx, k, mu, V, c(-2, -1, 1, 2, 0))
  st <- lapply(path, function(g) ns$.emc_shear_state(g$mu, g$V, lap, c$prior, spec))
  tr <- lapply(st[1:4], function(g) ns$.emc_shear_transport(aa, st[[5]], g, NULL)$alpha)
  for (s in seq_len(n)) {
    X <- t(vapply(tr, function(a) a[, s], numeric(p)))
    path_pts[[s]] <- rbind(path_pts[[s]], cbind(X, ll_rows(X, s)))
    path_rows <- path_rows + nrow(X)
  }
}
train_of <- function(s) {
  A <- do.call(rbind, lapply(chains, function(j) t(matrix(emc[[j]]$samples$alpha[, s, train_idx], p))))
  y <- unlist(lapply(chains, function(j) emc[[j]]$samples$subj_ll[s, train_idx]))
  A <- rbind(A, path_pts[[s]][, 1:p, drop = FALSE]); y <- c(y, path_pts[[s]][, p + 1])
  ok <- is.finite(y) & complete.cases(A)
  list(A = A[ok, , drop = FALSE], y = y[ok], is_path = c(rep(FALSE, sum(ok) - nrow(path_pts[[s]])), rep(TRUE, nrow(path_pts[[s]])))[seq_len(sum(ok))])
}

huber_w <- function(r, delta = 3) ifelse(abs(r) <= delta, 1, delta / abs(r))
# Full quadratic by IRLS Huber regression in whitened coordinates.
fit_quad <- function(tr) {
  A <- tr$A; y <- tr$y
  ac <- colMeans(A); R <- chol(cov(A) + diag(1e-8, p)); W <- backsolve(R, diag(p))
  Z <- sweep(A, 2, ac) %*% W
  ij <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  X <- cbind(1, Z, Z[, ij[, 1], drop = FALSE] * Z[, ij[, 2], drop = FALSE])
  w <- rep(1, length(y)); pen <- c(0, rep(1e-3, ncol(X) - 1))
  for (it in 1:8) {
    b <- solve(crossprod(X * w, X) + diag(pen), crossprod(X * w, y))[, 1]
    w <- huber_w(y - drop(X %*% b))
  }
  cq <- b[(p + 2):length(b)]
  B <- matrix(0, p, p); B[ij] <- -cq; B <- B + t(B); diag(B) <- -2 * cq[ij[, 1] == ij[, 2]]
  ev <- eigen(B, symmetric = TRUE); B <- ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
  H <- W %*% B %*% t(W); H <- (H + t(H)) / 2
  # Refit the linear part with the projected curvature held.
  off <- 0.5 * rowSums((A %*% H) * A)
  X1 <- cbind(1, A); w <- rep(1, length(y))
  for (it in 1:8) {
    b1 <- solve(crossprod(X1 * w, X1) + diag(c(0, rep(1e-8, p))), crossprod(X1 * w, y + off))[, 1]
    w <- huber_w(y + off - drop(X1 %*% b1))
  }
  list(H = H, lin = b1[-1])
}
# Skew factor on top of a Gaussian part (lin, H): fit (c, lambda) and a constant.
fit_skew <- function(tr, lin, H) {
  base <- tr$y - drop(tr$A %*% lin) + 0.5 * rowSums((tr$A %*% H) * tr$A)
  am <- colMeans(tr$A); Ac <- sweep(tr$A, 2, am)
  sc <- apply(Ac, 2, sd); sc[sc <= 0] <- 1
  Zs <- sweep(Ac, 2, sc, "/")
  loss <- function(th) {
    eta <- th[1] + drop(Zs %*% th[-1])
    r <- base - pnorm(eta, log.p = TRUE)
    r <- r - median(r)
    sum(ifelse(abs(r) <= 3, 0.5 * r^2, 3 * abs(r) - 4.5)) + 1e-3 * sum(th[-1]^2)
  }
  r0 <- base - median(base)
  best <- list(value = sum(ifelse(abs(r0) <= 3, 0.5 * r0^2, 3 * abs(r0) - 4.5)), par = c(10, numeric(p)))
  starts <- list(c(0, numeric(p)))
  for (j in seq_len(p)) for (sg in c(-1, 1)) starts[[length(starts) + 1]] <- c(1, replace(numeric(p), j, sg * 2))
  for (st in starts) {
    o <- tryCatch(optim(st, loss, method = "BFGS", control = list(maxit = 200)), error = function(e) NULL)
    if (!is.null(o) && o$value < best$value) best <- o
  }
  lam <- best$par[-1] / sc
  list(c = best$par[1] - sum(lam * am), lambda = lam, gain = loss(c(10, numeric(p))) - best$value)
}
# Held-out (every 5th point) check of each surrogate against plain Laplace.
surr <- function(f, A) {
  v <- drop(A %*% f$lin) - 0.5 * rowSums((A %*% f$H) * A)
  if (!is.null(f$lambda)) v <- v + pnorm(f$c + drop(A %*% f$lambda), log.p = TRUE)
  v
}
fits_file <- Sys.getenv("SKEW_FITS")
fits <- if (nzchar(fits_file) && file.exists(fits_file)) readRDS(fits_file) else lapply(seq_len(n), function(s) {
  tr <- train_of(s)
  L <- list(lin = lap$lin[, s], H = lap$cube[, , s])
  ho <- seq_along(tr$y) %% 5 == 0
  trn <- list(A = tr$A[!ho, , drop = FALSE], y = tr$y[!ho])
  Q <- tryCatch(fit_quad(trn), error = function(e) L)
  sk_L <- fit_skew(trn, L$lin, L$H); sk_Q <- fit_skew(trn, Q$lin, Q$H)
  cand <- list(laplace = L, quad = Q, laplace_skew = c(L, sk_L[c("c", "lambda")]),
               quad_skew = c(Q, sk_Q[c("c", "lambda")]))
  ho_rmse <- vapply(cand, function(f) {
    r <- tr$y[ho] - surr(f, tr$A[ho, , drop = FALSE]); sqrt(mean((r - median(r))^2))
  }, 1)
  # Gate: a richer surrogate replaces Laplace only if it predicts held-out points better.
  gated <- cand[[names(which.min(ho_rmse))]]
  # Refit on all points with the chosen form.
  list(cand = cand, gated = gated, gated_name = names(which.min(ho_rmse)), ho_rmse = ho_rmse)
})
cpu_fit <- sum((proc.time() - cpu_fit)[1:2])
if (nzchar(fits_file) && !file.exists(fits_file)) saveRDS(fits, fits_file)
cat("fit cpu", round(cpu_fit, 1), "s; path rows", path_rows, "\n")
print(table(vapply(fits, `[[`, "", "gated_name")))
print(round(t(apply(vapply(fits, `[[`, numeric(4), "ho_rmse"), 1, quantile, c(.1, .5, .9, 1))), 3))

# ---- Frames and transport -----------------------------------------------
# State of one subject's surrogate conditional at group state (mu, V).
frame <- function(f, mu, Vi) {
  P0 <- f$H + Vi; R0 <- chol((P0 + t(P0)) / 2)
  m0 <- drop(backsolve(R0, forwardsolve(t(R0), f$lin + drop(Vi %*% mu))))
  if (is.null(f$lambda) || all(f$lambda == 0))
    return(list(m0 = m0, R = R0, skew = FALSE))
  Sl <- drop(backsolve(R0, forwardsolve(t(R0), f$lambda)))
  vw <- 1 + sum(f$lambda * Sl)
  Q <- P0 + tcrossprod(f$lambda); R <- chol((Q + t(Q)) / 2)
  list(m0 = m0, R = R, skew = TRUE, Sl = Sl, vw = vw, mw = f$c + sum(f$lambda * m0))
}
# log density of the truncated (w > 0) normal N(m, v)
ltn <- function(w, m, v) dnorm(w, m, sqrt(v), log = TRUE) - pnorm(m / sqrt(v), log.p = TRUE)
cond_mean <- function(fr, w) fr$m0 + fr$Sl * (w - fr$mw) / fr$vw
transport <- function(a, w, f0, f1) {
  if (!f0$skew) return(list(a = f1$m0 + drop(backsolve(f1$R, f0$R %*% (a - f0$m0))), w = w,
                            jac = sum(log(diag(f0$R))) - sum(log(diag(f1$R)))))
  s0 <- sqrt(f0$vw); s1 <- sqrt(f1$vw)
  lu <- pnorm(w, f0$mw, s0, lower.tail = FALSE, log.p = TRUE) - pnorm(f0$mw / s0, log.p = TRUE)
  w1 <- qnorm(lu + pnorm(f1$mw / s1, log.p = TRUE), f1$mw, s1, lower.tail = FALSE, log.p = TRUE)
  a1 <- cond_mean(f1, w1) + drop(backsolve(f1$R, f0$R %*% (a - cond_mean(f0, w))))
  list(a = a1, w = w1, jac = ltn(w, f0$mw, f0$vw) - ltn(w1, f1$mw, f1$vw) +
         sum(log(diag(f0$R))) - sum(log(diag(f1$R))))
}
log_q <- function(f, a, w) if (is.null(f$lambda) || all(f$lambda == 0)) 0 else
  ltn(w, f$c + sum(f$lambda * a), 1)
draw_w <- function(f, a) {
  if (is.null(f$lambda) || all(f$lambda == 0)) return(0)
  m <- f$c + sum(f$lambda * a)
  qnorm(runif(1, pnorm(0, m, 1), 1), m, 1)
}

# ---- Evaluation at held-out draws -----------------------------------------
methods <- c("laplace", "quad", "laplace_skew", "quad_skew", "gated")
out <- list(); persub <- list(); cpu_eval <- proc.time()
for (j in chains) for (d in draws) {
  c <- emc[[j]]
  mu <- c$samples$theta_mu[, d]; V <- c$samples$theta_var[, , d]
  alpha <- matrix(c$samples$alpha[, , d], p)
  path <- ns$.emc_shear_coord_path(c, idx, k, mu, V, grid)
  lp <- vapply(path, function(g) ns$.emc_shear_log_prior_v(g$V, spec) + ns$.emc_shear_log_prior_mu(g$mu, c$prior) +
                 sum((spec$size + 1) * g$u_sd), 1)
  cV <- lapply(path, function(g) chol(g$V)); Vi <- lapply(cV, chol2inv)
  hier <- function(a, g) { z <- forwardsolve(t(cV[[g]]), a - path[[g]]$mu); -0.5 * sum(z^2) - sum(log(diag(cV[[g]]))) }
  for (m in methods) {
    f_of <- function(s) if (m == "gated") fits[[s]]$gated else fits[[s]]$cand[[m]]
    reps <- if (m %in% c("laplace", "quad")) 1L else REPS
    per <- array(0, c(length(grid), n, reps))
    for (s in seq_len(n)) {
      f <- f_of(s)
      frs <- lapply(seq_along(grid), function(g) frame(f, path[[g]]$mu, Vi[[g]]))
      X <- NULL; extra <- NULL
      for (r in seq_len(reps)) {
        w0 <- draw_w(f, alpha[, s])
        for (g in seq_along(grid)) {
          tr <- if (g == g0) list(a = alpha[, s], w = w0, jac = 0) else transport(alpha[, s], w0, frs[[g0]], frs[[g]])
          X <- rbind(X, tr$a)
          extra <- c(extra, hier(tr$a, g) + log_q(f, tr$a, tr$w) + tr$jac)
        }
      }
      ll <- ll_rows(X, s, c)
      per[, s, ] <- matrix(ll + extra, length(grid))
    }
    base <- per[g0, , , drop = FALSE]
    per <- per - base[rep(1L, length(grid)), , , drop = FALSE]
    tot <- apply(per, c(1, 3), sum) + (lp - lp[g0])
    tot <- matrix(tot, length(grid))
    acc <- rowMeans(pmin(exp(tot), 1))
    spread <- mean(apply(per, 3, function(P) sum(apply(P, 2, var))))
    row <- data.frame(chain = j, draw = d, method = m, spread = spread,
                      am1 = acc[grid == -1], amh = acc[grid == -0.5], aph = acc[grid == 0.5], ap1 = acc[grid == 1],
                      am2 = acc[grid == -2], ap2 = acc[grid == 2])
    print(row, digits = 3); out[[length(out) + 1]] <- row
    persub[[paste(j, d, m)]] <- apply(per, c(1, 2), mean)
  }
  saveRDS(list(results = do.call(rbind, out), persub = persub, fits = fits, coord = coord, k = k,
               draws = draws, chains = chains, cpu_fit = cpu_fit, path_rows = path_rows,
               head = system("git rev-parse HEAD", intern = TRUE)), args[2])
}
res <- do.call(rbind, out)
cat("\nmeans over draws\n")
print(aggregate(cbind(spread, am2, am1, amh, aph, ap1, ap2) ~ method, res, mean), digits = 3)
cat("\nmedian spread\n"); print(aggregate(spread ~ method, res, median), digits = 3)
cat("eval cpu", round(sum((proc.time() - cpu_eval)[1:2]), 1), "s\n")
