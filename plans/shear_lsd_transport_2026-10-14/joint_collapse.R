# Joint group update over line-collapsed subjects. Each subject is collapsed
# along its fixed line for the target coordinate (as line_collapse.R); the
# group update then moves a set S of group coordinates jointly (mean offsets
# and log-sd offsets with correlations held) by random-scan slice steps on
#   pi~(phi_S) = p(mu) p(V) J prod_s Z_s(phi),
# subjects follow by their line CDFs, and the exact correction is the
# interpolation residual. For each S it also reports the same marginal under
# every subject integrated out under its Gaussian frame, to separate subject
# pinning from coupling between group coordinates.
#
# Rscript joint_collapse.R <emc.rds> <chain> lsd:<par> <draws|auto> [ncorr=4] [iter=300] [reps=8] [K=16] [NC=8]
suppressPackageStartupMessages(library(EMC2, lib.loc = if (nzchar(Sys.getenv("EMC_LIB"))) Sys.getenv("EMC_LIB") else NULL))
ns <- asNamespace("EMC2")
a <- commandArgs(TRUE)
kv <- function(key, def) { m <- grep(paste0("^", key, "="), a, value = TRUE); if (length(m)) sub(".*=", "", m[1]) else def }
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(a[1])), quiet = TRUE)
j <- as.integer(a[2]); what <- sub(".*:", "", a[3])
K <- as.integer(kv("K", 16)); SPAN <- 5; NCORR <- as.integer(kv("ncorr", 4)); ITER <- as.integer(kv("iter", 300))
REPS <- as.integer(kv("reps", 8)); NC <- as.integer(kv("NC", 8)); M_SWEEP <- as.integer(kv("m", 5))
set.seed(as.integer(kv("seed", 1)))
ch0 <- emc[[j]]; p <- ch0$n_pars; n <- ch0$n_subjects; pn <- ch0$par_names; last <- ch0$samples$idx
draws <- if (a[4] == "auto") round(seq(last - 7 * 60, last, by = 60)) else as.integer(strsplit(a[4], ",")[[1]])
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample"); idx <- idx[idx <= last]
ctx <- ns$.emc_shear_context(ch0); H <- ctx$packed$cube; LIN <- ctx$packed$lin; spec <- ctx$spec; prior <- ch0$prior
l <- if (what %in% pn) match(what, pn) else as.integer(what)
# group coordinates 1..p = means, p+1..2p = log sds
gw <- rbind(matrix(ch0$samples$theta_mu[, idx], p), 0.5 * log(matrix(apply(ch0$samples$theta_var[, , idx, drop = FALSE], 3, diag), p)))
gsd <- apply(gw, 1, sd); gname <- c(paste0("mu ", pn), paste0("lsd ", pn)); tgt <- p + l
cr <- cor(t(gw))[tgt, ]; cr[tgt] <- 0
top <- order(abs(cr), decreasing = TRUE)[seq_len(NCORR)]
sets <- list(lsd = tgt, `lsd+mu` = c(tgt, l), corr = unique(c(tgt, l, top)), all = seq_len(2 * p))
cat(sprintf("%s chain %d: target log sd(%s), window sd %.3f; most correlated group coords: %s\n", basename(a[1]), j, pn[l], gsd[tgt],
            paste(sprintf("%s (%.2f)", gname[top], cr[top]), collapse = ", ")))
state_at <- function(mu, V, u) { um <- u[seq_len(p)]; d <- exp(u[p + seq_len(p)]); list(mu = mu + um, V = V * outer(d, d)) }
lp_u <- function(st, u) { g <- tryCatch(ns$.emc_shear_log_prior_mu(st$mu, prior) + ns$.emc_shear_log_prior_v(st$V, spec), error = function(e) -Inf)
  g + sum((spec$size + 1) * u[p + seq_len(p)]) }
# lines, fixed per block
mu_ref <- rowMeans(matrix(ch0$samples$theta_mu[, idx], p))
V_ref <- apply(ch0$samples$theta_var[, , idx, drop = FALSE], 1:2, mean); Vi_ref <- solve(V_ref); e <- replace(numeric(p), l, 1)
Cl <- matrix(0, p, n); Ell <- matrix(0, p, n); tc <- numeric(n); tsd <- numeric(n)
for (s in seq_len(n)) { Q <- H[, , s] + Vi_ref; c_s <- drop(solve(Q, Vi_ref %*% e)); q <- sum(c_s * (Q %*% c_s))
  Cl[, s] <- c_s; Ell[, s] <- drop(Q %*% c_s) / q; tc[s] <- sum(Ell[, s] * solve(Q, LIN[, s] + Vi_ref %*% mu_ref)); tsd[s] <- 1 / sqrt(q) }
support <- ns$.emc_shear_support(ch0)
feasible <- function(s, w) { lo <- -Inf; hi <- Inf
  for (h in support[[s]]) { uc <- sum(h$u * Cl[, s]); r <- h$b - sum(h$u * w)
    if (abs(uc) < 1e-14) next; if (uc > 0) lo <- max(lo, r / uc) else hi <- min(hi, r / uc) }
  c(lo, hi) }
llb <- function(X, s) ns$.emc_burn_ll_batch(X, ch0$data[[s]], ch0$model, pn)
ll_rows <- function(rows) unlist(parallel::mclapply(seq_len(n), function(s) list(llb(rows[[s]], s)), mc.cores = NC), recursive = FALSE)
one_row <- function(A) unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(A[, s], 1); colnames(X) <- pn; X })))
ldphi <- function(zr, zl) { out <- array(-Inf, dim(as.matrix(zr))); zr <- as.matrix(zr); zl <- as.matrix(zl); up <- zl >= 0
  i <- up & zr > zl; if (any(i)) { s1 <- pnorm(zl[i], lower.tail = FALSE, log.p = TRUE); s2 <- pnorm(zr[i], lower.tail = FALSE, log.p = TRUE); out[i] <- s1 + log1p(-exp(s2 - s1)) }
  i <- !up & zr > zl; if (any(i)) { f2 <- pnorm(zr[i], log.p = TRUE); f1 <- pnorm(zl[i], log.p = TRUE); out[i] <- f2 + log1p(-exp(f1 - f2)) }
  out }
lse <- function(x) { mx <- max(x); if (!is.finite(mx)) return(mx); mx + log(sum(exp(x - mx))) }
rlse <- function(X) { m <- apply(X, 1, max); m + log(rowSums(exp(X - m))) }
build_pieces <- function(tn, y, lo, hi, ylo, yhi) {
  Kn <- length(tn); b <- diff(y) / diff(tn); a0 <- y[-Kn] - b * tn[-Kn]
  left <- if (is.finite(lo) && abs(tn[1] - lo) < 1e-12) c(-Inf, tn[1], ylo, 0) else { bl <- max(b[1], 0); c(-Inf, tn[1], y[1] - bl * tn[1], bl) }
  right <- if (is.finite(hi) && abs(tn[Kn] - hi) < 1e-12) c(tn[Kn], Inf, yhi, 0) else { br <- min(b[Kn - 1], 0); c(tn[Kn], Inf, y[Kn] - br * tn[Kn], br) }
  rbind(left, cbind(tn[-Kn], tn[-1], a0, b), right) }
piece_lm <- function(P, q) { m <- (q[2] + P[, 4]) / q[1]; sA <- sqrt(q[1])
  q[3] + P[, 3] + (q[2] + P[, 4])^2 / (2 * q[1]) + 0.5 * log(2 * pi / q[1]) + drop(ldphi((P[, 2] - m) * sA, (P[, 1] - m) * sA)) }
piece_val <- function(P, t) { i <- which(t >= P[, 1] & t <= P[, 2])[1]; P[i, 3] + P[i, 4] * t }
cdf_side <- function(P, q, t) { lm <- piece_lm(P, q); lz <- lse(lm); i <- which(t >= P[, 1] & t <= P[, 2])[1]
  m <- (q[2] + P[i, 4]) / q[1]; sA <- sqrt(q[1]); zl <- (P[i, 1] - m) * sA; zr <- (P[i, 2] - m) * sA; zt <- (t - m) * sA
  base <- lm[i] - drop(ldphi(zr, zl)); L <- lse(c(lm[seq_len(i - 1)], base + drop(ldphi(zt, zl)))); R <- lse(c(lm[-seq_len(i)], base + drop(ldphi(zr, zt))))
  if (L <= R) list(side = "L", lf = L - lz) else list(side = "R", lf = R - lz) }
inv_cdf <- function(P, q, side, lf) { lm <- piece_lm(P, q); target <- lf + lse(lm)
  ord <- if (side == "L") seq_len(nrow(P)) else rev(seq_len(nrow(P))); acc <- -Inf
  for (i in ord) { nxt <- lse(c(acc, lm[i]))
    if (nxt >= target || i == ord[length(ord)]) {
      rem <- if (is.finite(acc)) target + log1p(-min(1, exp(acc - target))) else target
      frac <- min(rem - lm[i], 0); m <- (q[2] + P[i, 4]) / q[1]; sA <- sqrt(q[1])
      zl <- (P[i, 1] - m) * sA; zr <- (P[i, 2] - m) * sA; tot <- drop(ldphi(zr, zl))
      g <- if (side == "L") function(z) drop(ldphi(z, zl)) - tot - frac else function(z) frac - (drop(ldphi(zr, z)) - tot)
      return(m + uniroot(g, c(max(zl, -38), min(zr, 38)), tol = 1e-14)$root / sA) }
    acc <- nxt } }
# random-scan coordinate slice on f over coordinates S, widths w; nstep steps
slice_rs <- function(f, u, lu, S, nstep) {
  for (it in seq_len(nstep)) { i <- S[sample.int(length(S), 1)]; w0 <- gsd[i]
    y <- lu - rexp(1); L <- u; R <- u; L[i] <- u[i] - runif(1) * w0; R[i] <- L[i] + w0
    for (z in 1:10) { if (f(L) < y) break; L[i] <- L[i] - w0 }; for (z in 1:10) { if (f(R) < y) break; R[i] <- R[i] + w0 }
    repeat { v <- u; v[i] <- runif(1, L[i], R[i]); lv <- f(v); if (lv >= y) { u <- v; lu <- lv; break }; if (v[i] < u[i]) L[i] <- v[i] else R[i] <- v[i] } }
  list(u = u, lu = lu) }

for (d in draws) {
  t0 <- proc.time()[3]
  alpha <- matrix(ch0$samples$alpha[, , d], p); mu <- ch0$samples$theta_mu[, d]; V <- ch0$samples$theta_var[, , d]
  t_cur <- colSums(Ell * alpha); W <- alpha - sweep(Cl, 2, t_cur, `*`)
  fe <- t(vapply(seq_len(n), function(s) feasible(s, W[, s]), numeric(2)))
  nodes <- lapply(seq_len(n), function(s) { lo <- max(tc[s] - SPAN * tsd[s], fe[s, 1]); hi <- min(tc[s] + SPAN * tsd[s], fe[s, 2])
    list(tn = seq(lo, hi, length.out = K), extra = c(if (is.finite(fe[s, 1])) fe[s, 1] - tsd[s], if (is.finite(fe[s, 2])) fe[s, 2] + tsd[s]), lo = fe[s, 1], hi = fe[s, 2]) })
  Y <- ll_rows(lapply(seq_len(n), function(s) { X <- t(W[, s] + outer(Cl[, s], c(nodes[[s]]$tn, nodes[[s]]$extra))); colnames(X) <- pn; X }))
  Pl <- lapply(seq_len(n), function(s) { nd <- nodes[[s]]; y <- Y[[s]]; Kn <- K
    ylo <- if (is.finite(nd$lo)) y[Kn + 1] else NA; yhi <- if (is.finite(nd$hi)) y[Kn + 1 + is.finite(nd$lo)] else NA
    y0 <- y[seq_len(Kn)]; y0[!is.finite(y0)] <- -1e6; build_pieces(nd$tn, y0, nd$lo, nd$hi, ylo, yhi) })
  # padded piece arrays for the vectorised log Z
  Pm <- max(vapply(Pl, nrow, 1L)); PA <- array(0, c(n, Pm, 4))
  for (s in seq_len(n)) { P <- Pl[[s]]; PA[s, seq_len(nrow(P)), ] <- P; if (nrow(P) < Pm) PA[s, (nrow(P) + 1):Pm, 1:2] <- 0 }
  quads <- function(st) { Vi <- solve(st$V); ldV <- determinant(st$V)$modulus; r <- W - st$mu
    cbind(colSums(Cl * (Vi %*% Cl)), -colSums(Cl * (Vi %*% r)), -0.5 * colSums(r * (Vi %*% r)) - 0.5 * ldV - 0.5 * p * log(2 * pi)) }
  logZ <- function(q) { A <- q[, 1]; B <- q[, 2]; m <- (B + PA[, , 4]) / A; sA <- sqrt(A)
    lm <- q[, 3] + PA[, , 3] + (B + PA[, , 4])^2 / (2 * A) + 0.5 * log(2 * pi / A) + ldphi((PA[, , 2] - m) * sA, (PA[, , 1] - m) * sA)
    rlse(lm) }
  f_line <- function(u) { st <- state_at(mu, V, u); g <- lp_u(st, u); if (!is.finite(g)) return(-Inf); g + sum(logZ(quads(st))) }
  f_gauss <- function(u) { st <- state_at(mu, V, u); g <- lp_u(st, u); if (!is.finite(g)) return(-Inf)
    Vi <- solve(st$V); ldV <- determinant(st$V)$modulus; Vm <- drop(Vi %*% st$mu)
    g + sum(vapply(seq_len(n), function(s) { Q <- H[, , s] + Vi; h <- LIN[, s] + Vm
      -0.5 * ldV - 0.5 * determinant(Q)$modulus + 0.5 * sum(h * solve(Q, h)) - 0.5 * sum(st$mu * Vm) }, 1)) }
  f_held <- function(u) { st <- state_at(mu, V, u); g <- lp_u(st, u); if (!is.finite(g)) return(-Inf)
    R <- chol(st$V); z <- backsolve(R, alpha - st$mu, transpose = TRUE); g - sum(z^2) / 2 - n * sum(log(diag(R))) }
  cur_ll <- one_row(alpha); r0 <- cur_ll - vapply(seq_len(n), function(s) piece_val(Pl[[s]], t_cur[s]), 1)
  q0 <- quads(state_at(mu, V, numeric(2 * p)))
  out <- sprintf("draw %4d:", d)
  for (nm in names(sets)) { S <- sets[[nm]]
    # marginal sd of the target under each collapsed conditional (long slice chains)
    msd <- vapply(list(held = f_held, line = f_line, gauss = f_gauss), function(f) {
      u <- numeric(2 * p); lu <- f(u); tr <- numeric(ITER)
      for (it in seq_len(ITER)) { o <- slice_rs(f, u, lu, S, length(S)); u <- o$u; lu <- o$lu; tr[it] <- u[tgt] }
      sd(tr[-seq_len(ITER %/% 5)]) / gsd[tgt] }, 1)
    # the move: M_SWEEP * |S| random-scan steps on pi~, transport, exact correction
    acc <- numeric(REPS); dj <- numeric(REPS)
    for (r in seq_len(REPS)) {
      o <- slice_rs(f_line, numeric(2 * p), f_line(numeric(2 * p)), S, M_SWEEP * length(S))
      q1 <- quads(state_at(mu, V, o$u))
      t1 <- vapply(seq_len(n), function(s) { cs <- cdf_side(Pl[[s]], q0[s, ], t_cur[s]); inv_cdf(Pl[[s]], q1[s, ], cs$side, cs$lf) }, 1)
      ll1 <- one_row(W + sweep(Cl, 2, t1, `*`)); r1 <- ll1 - vapply(seq_len(n), function(s) piece_val(Pl[[s]], t1[s]), 1)
      acc[r] <- min(1, exp(sum(r1 - r0))); dj[r] <- o$u[tgt]
    }
    out <- paste0(out, sprintf(" [%s, |S|=%d] target sd held %.2f line %.2f gauss %.2f; move acc %.2f ESJD %.3f |",
      nm, length(S), msd[1], msd[2], msd[3], mean(acc), mean(acc * dj^2) / gsd[tgt]^2))
  }
  cat(out, sprintf("(%.0fs)\n", proc.time()[3] - t0))
}
