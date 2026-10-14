# Offline full-dimensional transport with a non-Gaussian first marginal.
# Usage: Rscript ridge_first.R fixture.rds output.rds [chains] [draws]
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2")
args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
# Reuse the analytically integrated piecewise exponential/Gaussian CDFs.
helper <- file.path("plans", "shear_lsd_transport_2026-10-14", "line_collapse.R")
needed <- c("build_pieces", "ldphi", "piece_logmass", "lse", "cdf_side", "inv_cdf", "log_ltilde")
for (ex in parse(helper)) if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.symbol(ex[[2]]) &&
                            as.character(ex[[2]]) %in% needed) eval(ex)
chains <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else 1L
draws <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else c(1001, 1000, 980, 950, 900, 800, 700, 600)
ch <- emc[[1]]; p <- ch$n_pars; n <- ch$n_subjects; pn <- ch$par_names
l <- match("B_Espeed", pn); stopifnot(!is.na(l))
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
train_idx <- setdiff(idx, draws)
ridge_solve <- function(X, y, penalty = 0.01) {
  pen <- rep(penalty, ncol(X)); pen[1] <- 0
  solve(crossprod(X) + diag(pen), crossprod(X, y))[, 1]
}
hinges <- function(x, knots) outer(x, knots, function(a, b) pmax(a - b, 0))
directions <- list()
for (nm in intersect(c("B", "B_Espeed", "B_lRright", "A", "v", "t0"), pn))
  directions[[nm]] <- replace(numeric(p), match(nm, pn), 1)
if (all(c("B", "B_Espeed", "B_lRright") %in% pn)) {
  for (speed in c(-1, 1)) for (resp in c(-1, 1)) {
    v <- numeric(p); v[match(c("B", "B_Espeed", "B_lRright"), pn)] <- c(1, speed, resp)
    directions[[paste0("threshold", speed, ":", resp)]] <- v
  }
}
fit_one <- function(s) {
  A <- do.call(rbind, lapply(emc, function(c) t(matrix(c$samples$alpha[, s, train_idx], p))))
  y <- unlist(lapply(emc, function(c) c$samples$subj_ll[s, train_idx]))
  keep <- complete.cases(A) & is.finite(y); A <- A[keep, ]; y <- y[keep]
  ac <- colMeans(A); R <- chol(cov(A) + diag(1e-10, p)); W <- solve(R)
  Z <- sweep(A, 2, ac) %*% W
  ij <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  X <- cbind(1, Z, Z[, ij[, 1], drop = FALSE] * Z[, ij[, 2], drop = FALSE])
  # Interleaved held-out rows select a ridge; replay draws are never fitted.
  va <- seq_len(nrow(A)) %% 5 == 0
  scores <- vapply(directions, function(v) {
    x <- drop(A %*% v); x <- (x - mean(x)) / sd(x)
    knots <- as.numeric(quantile(x[!va], seq(0.1, 0.9, length.out = 7)))
    XX <- cbind(X, hinges(x, knots)); b <- ridge_solve(XX[!va, ], y[!va])
    mean((y[va] - XX[va, ] %*% b)^2)
  }, 1)
  name <- names(which.min(scores)); v <- directions[[name]]
  raw <- drop(A %*% v); xc <- mean(raw); xs <- sd(raw); x <- (raw - xc) / xs
  knots <- as.numeric(quantile(x, seq(0.1, 0.9, length.out = 7)))
  F <- hinges(x, knots); b <- ridge_solve(cbind(X, F), y)
  cq <- b[(p + 2):(p + 1 + nrow(ij))]
  B <- matrix(0, p, p); B[ij] <- -cq; B <- B + t(B)
  diag(B) <- -2 * cq[ij[, 1] == ij[, 2]]
  ev <- eigen(B, symmetric = TRUE); B <- ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
  # Z = (A-ac) W, hence precision on the sampled scale is W B W'.
  H <- W %*% B %*% t(W)
  bb <- ridge_solve(cbind(1, Z, F), y + 0.5 * rowSums((Z %*% B) * Z))
  lin <- drop(H %*% ac + W %*% bb[1 + seq_len(p)])
  hb <- tail(bb, ncol(F)); hn <- xc + xs * knots
  hy <- drop(hinges(knots, knots) %*% hb)
  slopes <- diff(hy) / diff(hn)
  P <- rbind(c(-Inf, hn[1], hy[1], 0),
             cbind(head(hn, -1), tail(hn, -1), head(hy, -1) - slopes * head(hn, -1), slopes),
             c(tail(hn, 1), Inf, tail(hy, 1) - sum(hb) / xs * tail(hn, 1), sum(hb) / xs))
  # The first hinge is zero to the left of its first knot.
  list(prec = H, lin = lin, v = v, P = P, name = name,
       validation_mse = min(scores), train_rmse = sqrt(mean((y - (cbind(1, Z, F) %*% bb - 0.5 * rowSums((Z %*% B) * Z)))^2)))
}
fit_work <- function() {
  c <- emc[[1]]; ct <- ns$.emc_shear_context(c)
  probe_draws <- setdiff(c(600, 700, 800, 900), draws)
  if (length(probe_draws) < 3) probe_draws <- train_idx[unique(round(seq(1, length(train_idx), length.out = 4)))]
  gg <- c(-1, -0.5, 0, 0.5, 1); base <- which(gg == 0)
  props <- lapply(probe_draws, function(d) {
    mu <- c$samples$theta_mu[, d]; V <- c$samples$theta_var[, , d]
    aa <- matrix(c$samples$alpha[, , d], p)
    path <- ns$.emc_shear_coord_path(c, idx, p + l, mu, V, gg)
    st <- lapply(path, function(g) ns$.emc_shear_state(g$mu, g$V, ct$packed, c$prior, ct$spec))
    lapply(st, function(g) ns$.emc_shear_transport(aa, st[[base]], g, ct$walls)$alpha)
  })
  lapply(seq_len(n), function(s) {
    H <- ct$packed$cube[, , s]; lin <- ct$packed$lin[, s]
    A <- do.call(rbind, lapply(props, function(pr) t(vapply(pr, function(aa) aa[, s], numeric(p)))))
    yy <- ns$.emc_burn_ll_batch(A, c$data[[s]], c$model, pn) - drop(A %*% lin) + 0.5 * rowSums((A %*% H) * A)
    origin <- rep(seq(0, nrow(A) - 1, by = length(gg)) + base, each = length(gg))
    dy <- yy - yy[origin]
    va <- seq_len(nrow(A)) > nrow(A) - length(gg)
    candidates <- lapply(directions, function(v) {
      raw <- drop(A %*% v); xc <- mean(raw); xs <- sd(raw); x <- (raw - xc) / xs
      knots <- unique(as.numeric(quantile(x, seq(0.05, 0.95, length.out = 9))))
      F <- cbind(x, hinges(x, knots)); DX <- F - F[origin, ]
      b <- solve(crossprod(DX[!va, ]) + diag(0.1, ncol(DX)), crossprod(DX[!va, ], dy[!va]))[, 1]
      score <- mean((dy[va] - DX[va, ] %*% b)^2)
      b <- solve(crossprod(DX) + diag(0.1, ncol(DX)), crossprod(DX, dy))[, 1]
      hn <- xc + xs * knots; hy <- drop(cbind(knots, hinges(knots, knots)) %*% b)
      slopes <- diff(hy) / diff(hn); bl <- b[1] / xs; br <- sum(b) / xs
      P <- rbind(c(-Inf, hn[1], hy[1] - bl * hn[1], bl),
                 cbind(head(hn, -1), tail(hn, -1), head(hy, -1) - slopes * head(hn, -1), slopes),
                 c(tail(hn, 1), Inf, tail(hy, 1) - br * tail(hn, 1), br))
      list(prec = H, lin = lin, v = v, P = P, validation_mse = score,
           baseline_mse = mean(dy[va]^2), train_rmse = sqrt(mean((dy - DX %*% b)^2)))
    })
    name <- names(which.min(vapply(candidates, `[[`, 1, "validation_mse")))
    f <- candidates[[name]]; f$name <- name; f
  })
}
fit_time <- system.time(fits <- if (Sys.getenv("RIDGE_FIT") == "work") fit_work() else lapply(seq_len(n), fit_one))
cat("fit elapsed", fit_time[[3]], "selection\n"); print(table(vapply(fits, `[[`, "", "name")))
spec <- ns$.emc_shear_v_spec(ch)
make_state <- function(mu, V, s, gaussian = FALSE) {
  f <- fits[[s]]; Vi <- solve(V); S <- solve(f$prec + Vi)
  m <- drop(S %*% (f$lin + Vi %*% mu)); v <- f$v
  mx <- sum(v * m); sx2 <- sum(v * (S %*% v)); beta <- drop(S %*% v) / sx2
  U <- qr.Q(qr(cbind(v, diag(p))))[, -1, drop = FALSE]
  K <- crossprod(U, (S - sx2 * tcrossprod(beta)) %*% U)
  P <- if (gaussian) matrix(c(-Inf, Inf, 0, 0), 1) else f$P
  list(m = m, S = S, mx = mx, beta = beta, U = U, L = t(chol(K)), P = P,
       q = c(1 / sx2, mx / sx2, -0.5 * mx^2 / sx2 - 0.5 * log(2 * pi * sx2)))
}
map_one <- function(a, from, to, s) {
  v <- fits[[s]]$v; x <- sum(v * a)
  cf <- cdf_side(from$P, from$q, x); xx <- inv_cdf(to$P, to$q, cf$side, cf$lf)
  r <- forwardsolve(from$L, drop(crossprod(from$U, a - from$m - from$beta * (x - from$mx))))
  aa <- to$m + to$beta * (xx - to$mx) + drop(to$U %*% (to$L %*% r))
  logpx <- function(st, z) log_ltilde(st$P, z) + st$q[3] - 0.5 * st$q[1] * z^2 + st$q[2] * z - lse(piece_logmass(st$P, st$q))
  list(a = aa, jac = logpx(from, x) - logpx(to, xx) + sum(log(diag(to$L))) - sum(log(diag(from$L))))
}
grid <- seq(-2, 2, by = 0.5); g0 <- which(grid == 0)
out <- list(); cpu0 <- proc.time()
for (j in chains) for (d in draws) {
  c <- emc[[j]]; c$samples$idx <- d
  mu <- c$samples$theta_mu[, d]; V <- c$samples$theta_var[, , d]; alpha <- matrix(c$samples$alpha[, , d], p)
  path <- ns$.emc_shear_coord_path(c, idx, p + l, mu, V, grid)
  fixed <- ns$.emc_shear_orbit_terms(c, idx, p + l, grid, function(st) ns$.emc_shear_context(c)$packed, 1L)
  for (method in c("ridge", "gaussian")) {
    started <- proc.time(); states <- lapply(path, function(st) lapply(seq_len(n), function(s) make_state(st$mu, st$V, s, method == "gaussian")))
    mapped <- lapply(seq_along(grid), function(g) lapply(seq_len(n), function(s) map_one(alpha[, s], states[[g0]][[s]], states[[g]][[s]], s)))
    per <- vapply(seq_len(n), function(s) {
      X <- t(vapply(mapped, function(tr) tr[[s]]$a, numeric(p)))
      ll <- ns$.emc_burn_ll_batch(X, c$data[[s]], c$model, pn)
      vapply(seq_along(grid), function(g) {
        z <- forwardsolve(t(chol(path[[g]]$V)), X[g, ] - mu)
        ll[g] - 0.5 * sum(z^2) - sum(log(diag(chol(path[[g]]$V)))) + mapped[[g]][[s]]$jac
      }, 1)
    }, numeric(length(grid)))
    lp <- vapply(path, function(st) ns$.emc_shear_log_prior_v(st$V, spec) + sum((spec$size + 1) * st$u_sd), 1)
    delta <- rowSums(per) + lp; delta <- delta - delta[g0]
    back <- vapply(seq_len(n), function(s) max(abs(map_one(mapped[[length(grid)]][[s]]$a, states[[length(grid)]][[s]], states[[g0]][[s]], s)$a - alpha[, s])), 1)
    row <- data.frame(chain = j, draw = d, method = method, spread = sum(apply(per, 2, var)),
                      am1 = min(1, exp(delta[grid == -1])), amh = min(1, exp(delta[grid == -0.5])),
                      aph = min(1, exp(delta[grid == 0.5])), ap1 = min(1, exp(delta[grid == 1])),
                      elapsed = (proc.time() - started)[[3]], cpu = sum((proc.time() - started)[1:2]), inverse = max(back))
    print(row); out[[length(out) + 1]] <- row
  }
  dl <- rowSums(fixed$per) + fixed$lp; dl <- dl - dl[g0]
  out[[length(out) + 1]] <- data.frame(chain = j, draw = d, method = "fixed", spread = sum(apply(fixed$per, 2, var)),
    am1 = min(1, exp(dl[grid == -1])), amh = min(1, exp(dl[grid == -0.5])), aph = min(1, exp(dl[grid == 0.5])), ap1 = min(1, exp(dl[grid == 1])), elapsed = NA, cpu = NA, inverse = NA)
  saveRDS(list(results = do.call(rbind, out), fits = fits, fit_time = fit_time, grid = grid,
               train_idx = train_idx, head = system("git rev-parse HEAD", intern = TRUE)), args[2])
}
print(aggregate(cbind(spread, am1, amh, aph, ap1) ~ method, do.call(rbind, out), mean))
