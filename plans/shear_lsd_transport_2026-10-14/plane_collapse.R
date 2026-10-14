# Emulation of plane-collapsed shear rounds: each subject moves in a fixed
# plane a = w + C t (w held), its likelihood known on a K x K node grid.
#
# Surrogate per subject: log Ltilde(t1, t2) is piecewise linear in t2 on each
# node column and linear in t1 between columns, so along t2 at any t1 the
# surrogate times the exact conditional group Gaussian G2(t2 | t1) integrates
# in closed form to Z2(t1). The t1 marginal is G1(t1) Zhat2(t1), with
# log Zhat2 interpolated on a fine likelihood-free grid. The move draws the
# group coordinate from pi~ = prior * prod_s int G1 Zhat2, carries every
# subject by the Knothe-Rosenblatt map of its surrogate (t1 marginal, then t2
# given t1), and corrects by prod_s [L / Ltilde * Z2 / Zhat2] at the new over
# the old point.
#
# Rscript plane_collapse.R <emc.rds> <chain> lsd:<par>|mean:<k>|mean:slow <draws|auto>
#   [K=12] [span=4] [c2=auto|axis|ridge|e<j>|f<j>] [reps=8] [NC=8] [predict=1]
suppressPackageStartupMessages(library(EMC2, lib.loc = if (nzchar(Sys.getenv("EMC_LIB"))) Sys.getenv("EMC_LIB") else NULL))
ns <- asNamespace("EMC2")
a <- commandArgs(TRUE)
kv <- function(key, def) { m <- grep(paste0("^", key, "="), a, value = TRUE); if (length(m)) sub(".*=", "", m[1]) else def }
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(a[1])), quiet = TRUE)
j <- as.integer(a[2]); coord <- a[3]
K <- as.integer(kv("K", 12)); SPAN <- as.numeric(kv("span", 4)); C2 <- kv("c2", "auto")
REPS <- as.integer(kv("reps", 8)); NC <- as.integer(kv("NC", 8)); M_SLICE <- as.integer(kv("m", 5))
FINE <- as.numeric(kv("fine", 0.2)); PRED_ONLY <- kv("predict", "0") == "1"
set.seed(as.integer(kv("seed", 1)))
ch0 <- emc[[j]]; p <- ch0$n_pars; n <- ch0$n_subjects; pn <- ch0$par_names
last <- ch0$samples$idx
draws <- if (a[4] == "auto") round(seq(last - 7 * 60, last, by = 60)) else as.integer(strsplit(a[4], ",")[[1]])
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample"); idx <- idx[idx <= last]
ctx <- ns$.emc_shear_context(ch0)
H <- ctx$packed$cube; LIN <- ctx$packed$lin; spec <- ctx$spec; prior <- ch0$prior
mu_ref <- rowMeans(matrix(ch0$samples$theta_mu[, idx], p))
V_ref <- apply(ch0$samples$theta_var[, , idx, drop = FALSE], 1:2, mean); Vi_ref <- solve(V_ref)
plan <- ch0$shear$plan
if (is.null(plan$dirs)) plan <- ns$.emc_shear_plan(ns$.emc_shear_fits(ch0), prior, V_ref)
kind <- sub(":.*", "", coord); what <- sub(".*:", "", coord)
if (kind == "lsd") {
  l <- if (what %in% pn) match(what, pn) else as.integer(what)
  e <- replace(numeric(p), l, 1); wsd <- sd(0.5 * log(ch0$samples$theta_var[l, l, idx]))
  state_at <- function(mu, V, u) { d <- rep(1, p); d[l] <- exp(u); list(mu = mu, V = V * outer(d, d)) }
  jac_u <- function(u) (spec$size[l] + 1) * u; label <- sprintf("log sd(%s)", pn[l])
} else {
  kk <- if (what == "slow") which.max(plan$lam) else as.integer(what)
  dvec <- plan$dirs[, kk]; e <- dvec; l <- NA
  wsd <- sd(drop(crossprod(dvec, matrix(ch0$samples$theta_mu[, idx], p))) / sum(dvec^2))
  state_at <- function(mu, V, u) list(mu = mu + u * dvec, V = V); jac_u <- function(u) 0
  label <- sprintf("mean direction %d (lambda %.3f)", kk, plan$lam[kk])
}
lp_group <- function(st) ns$.emc_shear_log_prior_mu(st$mu, prior) + ns$.emc_shear_log_prior_v(st$V, spec)
Qref <- lapply(seq_len(n), function(s) H[, , s] + Vi_ref)
mref <- lapply(seq_len(n), function(s) drop(solve(Qref[[s]], LIN[, s] + Vi_ref %*% mu_ref)))
dir_frame <- function(s, ev) drop(solve(Qref[[s]], Vi_ref %*% ev))
c1_of <- function(s) dir_frame(s, e)

# ---- Gaussian-frame predictor of pi~ width for any collapsed subspace --------
# Subject s collapsed over span(C) with t = Lproj' a, w = a - C t; the frame
# likelihood exp(lin'a - a'Ha/2) integrated against N(a | mu, V) over t.
ug <- seq(-4, 4, by = 0.25) * wsd
sd_of <- function(lg) { w <- exp(lg - max(lg)); w <- w / sum(w); sqrt(sum(w * ug^2) - sum(w * ug)^2) / wsd }
predict_width <- function(d, Cs) {   # Cs: list of p x r matrices (NULL = held, "full" = everything)
  alpha <- matrix(ch0$samples$alpha[, , d], p); mu <- ch0$samples$theta_mu[, d]; V <- ch0$samples$theta_var[, , d]
  pre <- lapply(seq_len(n), function(s) { C <- Cs[[s]]; if (is.null(C)) return(NULL)
    Lp <- Qref[[s]] %*% C %*% solve(crossprod(C, Qref[[s]] %*% C)); w <- alpha[, s] - drop(C %*% crossprod(Lp, alpha[, s])); list(C = C, w = w) })
  sd_of(vapply(ug, function(u) { st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g)) return(-Inf); Vi <- solve(st$V); ldV <- determinant(st$V)$modulus
    g + jac_u(u) + sum(vapply(seq_len(n), function(s) { x <- pre[[s]]
      if (is.null(x)) { r <- alpha[, s] - st$mu; return(-0.5 * sum(r * (Vi %*% r)) - 0.5 * ldV) }
      C <- x$C; w <- x$w; r <- w - st$mu
      A <- crossprod(C, (H[, , s] + Vi) %*% C); b <- crossprod(C, LIN[, s] - H[, , s] %*% w - Vi %*% r)
      -0.5 * sum(r * (Vi %*% r)) - 0.5 * ldV + 0.5 * sum(b * solve(A, b)) - 0.5 * determinant(A)$modulus }, 1)) }, 1))
}
cands <- c("axis", "ridge", paste0("f", setdiff(seq_len(p), if (is.na(l)) integer(0) else l)), paste0("e", setdiff(seq_len(p), if (is.na(l)) integer(0) else l)))
c2_of <- function(s, nm) {
  if (nm == "axis") return(e / sqrt(sum(e^2)))
  if (nm == "ridge") { ev <- eigen(H[, , s] + 1e-8 * diag(p), symmetric = TRUE); return(ev$vectors[, p]) }
  jj <- as.integer(substring(nm, 2)); ej <- replace(numeric(p), jj, 1)
  if (substr(nm, 1, 1) == "f") dir_frame(s, ej) else ej
}
plane_of <- function(nm) lapply(seq_len(n), function(s) cbind(c1_of(s), c2_of(s, nm)))
line_of <- lapply(seq_len(n), function(s) matrix(c1_of(s), p))
pd <- draws[unique(round(seq(1, length(draws), length.out = 3)))]
cat(sprintf("%s: chain %d, %s; window sd %.4f; n=%d p=%d\n", basename(a[1]), j, label, wsd, n, p))
cat("Gaussian-predicted pi~ sd (window sds) over draws", paste(pd, collapse = ","), "\n")
pw <- function(Cs) mean(vapply(pd, function(d) predict_width(d, Cs), 1))
cat(sprintf("  held %.2f | 1-D frame line %.2f | full collapse %.2f\n", pw(vector("list", n)), pw(line_of), pw(lapply(seq_len(n), function(s) diag(p)))))
if (C2 == "auto") {
  sc <- vapply(cands, function(nm) pw(plane_of(nm)), 1)
  print(round(sort(sc, decreasing = TRUE), 2))
  C2 <- names(which.max(sc))
} else cat(sprintf("  plane (frame, %s) %.2f\n", C2, pw(plane_of(C2))))
cat("chosen second direction:", C2, "\n")
if (PRED_ONLY) quit(save = "no")

# ---- plane surrogate --------------------------------------------------------
Cp <- list(); Lp <- list(); tc <- matrix(0, 2, n)
for (s in seq_len(n)) {   # whiten so t ~ N(tc, I) under the reference frame
  C <- plane_of(C2)[[s]]; R <- chol(crossprod(C, Qref[[s]] %*% C)); C <- C %*% solve(R)
  Cp[[s]] <- C; Lp[[s]] <- Qref[[s]] %*% C; tc[, s] <- drop(crossprod(Lp[[s]], mref[[s]]))
}
llb <- function(X, s) ns$.emc_burn_ll_batch(X, ch0$data[[s]], ch0$model, pn)
ll_rows <- function(rows) unlist(parallel::mclapply(seq_len(n), function(s) list(llb(rows[[s]], s)), mc.cores = NC), recursive = FALSE)
one_row <- function(A) unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(A[, s], 1); colnames(X) <- pn; X })))
ldphi <- function(zr, zl) {
  out <- rep(-Inf, length(zr)); up <- zl >= 0
  i <- up & zr > zl
  if (any(i)) { s1 <- pnorm(zl[i], lower.tail = FALSE, log.p = TRUE); s2 <- pnorm(zr[i], lower.tail = FALSE, log.p = TRUE); out[i] <- s1 + log1p(-exp(s2 - s1)) }
  i <- !up & zr > zl
  if (any(i)) { f2 <- pnorm(zr[i], log.p = TRUE); f1 <- pnorm(zl[i], log.p = TRUE); out[i] <- f2 + log1p(-exp(f1 - f2)) }
  out
}
lse <- function(x) { mx <- max(x); if (!is.finite(mx)) return(mx); mx + log(sum(exp(x - mx))) }
# pieces [L, R, a, b]: log value a + b t; tails extend the end slopes outward-decreasing
pieces <- function(tn, y) {
  Kn <- length(tn); b <- diff(y) / diff(tn); a0 <- y[-Kn] - b * tn[-Kn]
  bl <- max(b[1], 0); br <- min(b[Kn - 1], 0)
  rbind(c(-Inf, tn[1], y[1] - bl * tn[1], bl), cbind(tn[-Kn], tn[-1], a0, b), c(tn[Kn], Inf, y[Kn] - br * tn[Kn], br))
}
# q = (A, B, C): log Gaussian factor C - A t^2 / 2 + B t
piece_lm <- function(P, q) { m <- (q[2] + P[, 4]) / q[1]; sA <- sqrt(q[1])
  q[3] + P[, 3] + (q[2] + P[, 4])^2 / (2 * q[1]) + 0.5 * log(2 * pi / q[1]) + ldphi((P[, 2] - m) * sA, (P[, 1] - m) * sA) }
piece_val <- function(P, t) { i <- which(t >= P[, 1] & t <= P[, 2])[1]; P[i, 3] + P[i, 4] * t }
cdf_side <- function(P, q, t) {
  lm <- piece_lm(P, q); lz <- lse(lm); i <- which(t >= P[, 1] & t <= P[, 2])[1]
  m <- (q[2] + P[i, 4]) / q[1]; sA <- sqrt(q[1]); zl <- (P[i, 1] - m) * sA; zr <- (P[i, 2] - m) * sA; zt <- (t - m) * sA
  base <- lm[i] - ldphi(zr, zl)
  ll <- lse(c(lm[seq_len(i - 1)], base + ldphi(zt, zl))); lr <- lse(c(lm[-seq_len(i)], base + ldphi(zr, zt)))
  if (ll <= lr) list(side = "L", lf = ll - lz) else list(side = "R", lf = lr - lz)
}
inv_cdf <- function(P, q, side, lf) {
  lm <- piece_lm(P, q); target <- lf + lse(lm)
  ord <- if (side == "L") seq_len(nrow(P)) else rev(seq_len(nrow(P))); acc <- -Inf
  for (i in ord) {
    nxt <- lse(c(acc, lm[i]))
    if (nxt >= target || i == ord[length(ord)]) {
      rem <- if (is.finite(acc)) target + log1p(-min(1, exp(acc - target))) else target
      frac <- min(rem - lm[i], 0); m <- (q[2] + P[i, 4]) / q[1]; sA <- sqrt(q[1])
      zl <- (P[i, 1] - m) * sA; zr <- (P[i, 2] - m) * sA; tot <- ldphi(zr, zl)
      g <- if (side == "L") function(z) ldphi(z, zl) - tot - frac else function(z) frac - (ldphi(zr, z) - tot)
      return(m + uniroot(g, c(max(zl, -38), min(zr, 38)), tol = 1e-14)$root / sA)
    }
    acc <- nxt
  }
}
logq <- function(P, q, t) piece_val(P, t) + q[3] - 0.5 * q[1] * t^2 + q[2] * t - lse(piece_lm(P, q))

tn <- seq(-SPAN, SPAN, length.out = K)          # node grid in whitened t (offset by tc)
T1 <- seq(-SPAN - 3, SPAN + 3, by = FINE)       # fine likelihood-free grid for log Zhat2
# column pieces of subject s at t1 (linear in t1 between columns, flat beyond)
col_at <- function(cols, t1) {
  x <- (t1 - tn[1]) / (tn[2] - tn[1]); i <- min(max(floor(x), 0), K - 2); lam <- min(max(x - i, 0), 1)
  P <- cols[[i + 1]]; P[, 3:4] <- (1 - lam) * cols[[i + 1]][, 3:4] + lam * cols[[i + 2]][, 3:4]; P
}
# group Gaussian on the plane: log G(t) = C0 - t'At/2 + b't, in subject coords t (not offset)
plane_quad <- function(s, w, st, Vi, ldV) {
  C <- Cp[[s]]; r <- w - st$mu
  list(A = crossprod(C, Vi %*% C), b = -drop(crossprod(C, Vi %*% r)), C0 = -0.5 * sum(r * (Vi %*% r)) - 0.5 * ldV - 0.5 * p * log(2 * pi))
}
# everything for one subject at one group state: t1 pieces of log Zhat2 with G1 quad
subj_state <- function(s, g, cols) {
  A <- g$A; b <- g$b
  q2 <- function(t1) { m2 <- (b[2] - A[2, 1] * t1) / A[2, 2]; c(A[2, 2], A[2, 2] * m2, -0.5 * A[2, 2] * m2^2 + 0.5 * log(A[2, 2] / (2 * pi))) }
  # t2 pieces live in coordinates t2 (absolute): shift node grid by tc
  lz2 <- vapply(T1 + tc[1, s], function(t1) lse(piece_lm(col_at(cols, t1 - tc[1, s]), q2(t1))), 1)
  al <- A[1, 1] - A[1, 2]^2 / A[2, 2]; be <- b[1] - A[1, 2] * b[2] / A[2, 2]
  q1 <- c(al, be, g$C0 + 0.5 * b[2]^2 / A[2, 2] + 0.5 * log(2 * pi / A[2, 2]))
  list(P1 = pieces(T1 + tc[1, s], lz2), q1 = q1, q2 = q2)
}

run_draw <- function(d) {
  alpha <- matrix(ch0$samples$alpha[, , d], p); mu <- ch0$samples$theta_mu[, d]; V <- ch0$samples$theta_var[, , d]
  tcur <- vapply(seq_len(n), function(s) drop(crossprod(Lp[[s]], alpha[, s])), numeric(2))
  W <- alpha - vapply(seq_len(n), function(s) drop(Cp[[s]] %*% tcur[, s]), numeric(p))
  rows <- lapply(seq_len(n), function(s) { G <- as.matrix(expand.grid(t2 = tn + tc[2, s], t1 = tn + tc[1, s]))[, 2:1]
    X <- t(W[, s] + Cp[[s]] %*% t(G)); colnames(X) <- pn; X })
  el0 <- proc.time()[3]; Y <- ll_rows(rows); el <- proc.time()[3] - el0
  cols <- lapply(seq_len(n), function(s) { y <- matrix(Y[[s]], K, byrow = TRUE)   # row i = column t1_i
    y[!is.finite(y)] <- -1e6; lapply(seq_len(K), function(i) pieces(tn + tc[2, s], y[i, ])) })
  cur_ll <- one_row(alpha)
  states <- function(u) { st <- state_at(mu, V, u); g0 <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g0)) return(NULL); Vi <- solve(st$V); ldV <- determinant(st$V)$modulus
    list(lp = g0 + jac_u(u), ss = lapply(seq_len(n), function(s) subj_state(s, plane_quad(s, W[, s], st, Vi, ldV), cols[[s]]))) }
  lpt_of <- function(S) if (is.null(S)) -Inf else S$lp + sum(vapply(S$ss, function(x) lse(piece_lm(x$P1, x$q1)), 1))
  memo <- new.env(); getS <- function(u) { k <- sprintf("%.12g", u); if (is.null(memo[[k]])) memo[[k]] <- states(u); memo[[k]] }
  lpt <- function(u) lpt_of(getS(u))
  # correction term log(L / Ltilde) + log(Z2 / Zhat2) of subject s at t under state x
  corr <- function(s, x, t, ll) {
    t1 <- t[1]; P2 <- col_at(cols[[s]], t1 - tc[1, s]); lz2 <- lse(piece_lm(P2, x$q2(t1)))
    ll - piece_val(P2, t[2]) + lz2 - piece_val(x$P1, t1)
  }
  S0 <- getS(0); r0 <- vapply(seq_len(n), function(s) corr(s, S0$ss[[s]], tcur[, s], cur_ll[s]), 1)
  transport <- function(Sf, St, tf) vapply(seq_len(n), function(s) {
    xf <- Sf$ss[[s]]; xt <- St$ss[[s]]
    c1 <- cdf_side(xf$P1, xf$q1, tf[1, s]); t1 <- inv_cdf(xt$P1, xt$q1, c1$side, c1$lf)
    Pf <- col_at(cols[[s]], tf[1, s] - tc[1, s]); Pt <- col_at(cols[[s]], t1 - tc[1, s])
    c2 <- cdf_side(Pf, xf$q2(tf[1, s]), tf[2, s]); c(t1, inv_cdf(Pt, xt$q2(t1), c2$side, c2$lf)) }, numeric(2))
  move <- function(u) {
    S1 <- getS(u); if (is.null(S1)) return(list(lr = -Inf))
    t1 <- transport(S0, S1, tcur)
    A1 <- W + vapply(seq_len(n), function(s) drop(Cp[[s]] %*% t1[, s]), numeric(p))
    ll1 <- one_row(A1); r1 <- vapply(seq_len(n), function(s) corr(s, S1$ss[[s]], t1[, s], ll1[s]), 1)
    list(t1 = t1, rr = r1 - r0, lr = lpt(u) - lpt(0) + sum(r1 - r0))
  }
  # checks: inverse map
  S1 <- getS(wsd); tf <- transport(S0, S1, tcur); tb <- transport(S1, S0, tf); inv_err <- max(abs(tb - tcur))
  lg <- vapply(ug, lpt, 1); mob <- sd_of(lg)
  jumps <- c(-1, -0.5, 0.5, 1)
  jacc <- vapply(jumps, function(z) min(1, exp(move(z * wsd)$lr)), 1)
  slice <- function(u0, lu0, w0) { y <- lu0 - rexp(1); L <- u0 - runif(1) * w0; R <- L + w0
    for (i in 1:10) { if (lpt(L) < y) break; L <- L - w0 }; for (i in 1:10) { if (lpt(R) < y) break; R <- R + w0 }
    repeat { u1 <- runif(1, L, R); l1 <- lpt(u1); if (l1 >= y) return(c(u1, l1)); if (u1 < u0) L <- u1 else R <- u1 } }
  props <- t(vapply(seq_len(REPS), function(r) { u <- 0; lu <- lpt(0)
    for (m in seq_len(M_SLICE)) { o <- slice(u, lu, wsd); u <- o[1]; lu <- o[2] }
    mv <- move(u); rs <- sum(mv$rr); c(u = u, acc = min(1, exp(rs)), rs = rs, rsd = sd(mv$rr), rmax = max(abs(mv$rr))) }, numeric(5)))
  list(mob = mob, jacc = jacc, props = props, el = el, inv = inv_err, r0sd = sd(r0), r0max = max(abs(r0)))
}

cat(sprintf("plane surrogate: K=%d (%d rows/subject), span %g, fine grid %d\n", K, K * K, SPAN, length(T1)))
tot <- NULL
for (d in draws) {
  t0 <- proc.time()[3]; o <- run_draw(d); pr <- o$props
  esjd <- mean(pr[, "acc"] * pr[, "u"]^2) / wsd^2
  cat(sprintf("draw %4d: pi~ sd %.2f | jump acc -1 %.2f -.5 %.2f +.5 %.2f +1 %.2f | slice props acc mean %.2f med %.2f, resid sum sd %.2f per-subj sd %.3f max %.2f | cur corr sd %.3f max %.2f | ESJD/wsd^2 %.3f | node call %.2fs | inv err %.1e (%.0fs)\n",
    d, o$mob, o$jacc[1], o$jacc[2], o$jacc[3], o$jacc[4], mean(pr[, "acc"]), median(pr[, "acc"]), sd(pr[, "rs"]), mean(pr[, "rsd"]), max(pr[, "rmax"]), o$r0sd, o$r0max, esjd, o$el, o$inv, proc.time()[3] - t0))
  tot <- rbind(tot, c(o$mob, o$jacc, mean(pr[, "acc"]), esjd))
}
colnames(tot) <- c("pi_sd", "j-1", "j-.5", "j+.5", "j+1", "acc", "esjd")
cat("\nmeans over draws:\n"); print(round(colMeans(tot), 3))
