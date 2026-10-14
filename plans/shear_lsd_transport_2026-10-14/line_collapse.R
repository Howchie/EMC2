# Emulation of line-collapsed shear rounds on a replay or checkpoint.
#
# Each subject is moved along one line a = w + t c (w held), its likelihood
# known at K nodes on the line and interpolated log-linearly; the group
# coordinate is drawn from the collapsed target pi~(u) = prior * prod_s Z_s(u)
# and every subject is carried by its own line CDF. The exact correction is
# the interpolation residual at the old and new points. Compared against
# fixed-frame shear rounds (the current move) on the same draws.
#
# Rscript line_collapse.R <emc.rds> <chain> <coord: lsd:<par>|mean:<k>|mean:slow>
#   <draws a,b,..> [K=16] [span=5] [dir=frame|axis] [reps=12] [NC=8]
suppressPackageStartupMessages(library(EMC2, lib.loc = if (nzchar(Sys.getenv("EMC_LIB"))) Sys.getenv("EMC_LIB") else NULL))
ns <- asNamespace("EMC2")
a <- commandArgs(TRUE)
kv <- function(key, def) { m <- grep(paste0("^", key, "="), a, value = TRUE); if (length(m)) sub(".*=", "", m[1]) else def }
emc <- readRDS(a[1])
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
j <- as.integer(a[2]); coord <- a[3]
K <- as.integer(kv("K", 16)); SPAN <- as.numeric(kv("span", 5)); DIR <- kv("dir", "frame")
REPS <- as.integer(kv("reps", 12)); NC <- as.integer(kv("NC", 8)); M_SLICE <- as.integer(kv("m", 5))
SPACING <- kv("spacing", "uniform")
set.seed(as.integer(kv("seed", 1)))
ch0 <- emc[[j]]; p <- ch0$n_pars; n <- ch0$n_subjects; pn <- ch0$par_names
last <- ch0$samples$idx
draws <- if (a[4] == "auto") round(seq(last - 7 * 60, last, by = 60)) else as.integer(strsplit(a[4], ",")[[1]])
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample")
idx <- idx[idx <= last]
ctx <- ns$.emc_shear_context(ch0)
H <- ctx$packed$cube; LIN <- ctx$packed$lin
spec <- ctx$spec; prior <- ch0$prior
mu_ref <- rowMeans(matrix(ch0$samples$theta_mu[, idx], p))
V_ref <- apply(ch0$samples$theta_var[, , idx, drop = FALSE], 1:2, mean)
Vi_ref <- solve(V_ref)
plan <- ch0$shear$plan
if (is.null(plan$dirs)) plan <- ns$.emc_shear_plan(ns$.emc_shear_fits(ch0), prior, V_ref)

# ---- coordinate ----------------------------------------------------------
kind <- sub(":.*", "", coord); what <- sub(".*:", "", coord)
if (kind == "lsd") {
  l <- if (what %in% pn) match(what, pn) else as.integer(what)
  kk <- p + l
  e <- replace(numeric(p), l, 1)
  win <- 0.5 * log(ch0$samples$theta_var[l, l, idx])
  wsd <- sd(win)
  state_at <- function(mu, V, u) { d <- rep(1, p); d[l] <- exp(u); list(mu = mu, V = V * outer(d, d)) }
  jac_u <- function(u) (spec$size[l] + 1) * u
  label <- sprintf("log sd(%s)", pn[l])
} else {
  kk <- if (what == "slow") which.max(plan$lam) else as.integer(what)
  dvec <- plan$dirs[, kk]
  e <- dvec
  win <- drop(crossprod(dvec, matrix(ch0$samples$theta_mu[, idx], p))) / sum(dvec^2)
  wsd <- sd(win)
  state_at <- function(mu, V, u) list(mu = mu + u * dvec, V = V)
  jac_u <- function(u) 0
  label <- sprintf("mean direction %d (lambda %.3f)", kk, plan$lam[kk])
}
lp_group <- function(st) ns$.emc_shear_log_prior_mu(st$mu, prior) + ns$.emc_shear_log_prior_v(st$V, spec)

# ---- lines, fixed per block -----------------------------------------------
Cl <- matrix(0, p, n); Ell <- matrix(0, p, n); tc <- numeric(n); tsd <- numeric(n)
for (s in seq_len(n)) {
  Q <- H[, , s] + Vi_ref
  c_s <- if (DIR == "axis") e / sqrt(sum(e^2)) else drop(solve(Q, Vi_ref %*% e))
  q <- drop(crossprod(c_s, Q %*% c_s))
  Cl[, s] <- c_s; Ell[, s] <- drop(Q %*% c_s) / q
  m_ref <- drop(solve(Q, LIN[, s] + Vi_ref %*% mu_ref))
  tc[s] <- sum(Ell[, s] * m_ref); tsd[s] <- 1 / sqrt(q)
}
support <- ns$.emc_shear_support(ch0)
# Feasible t interval of line w + t c under all half-spaces u'a >= b.
feasible <- function(s, w) {
  lo <- -Inf; hi <- Inf
  for (h in support[[s]]) {
    uc <- sum(h$u * Cl[, s]); r <- h$b - sum(h$u * w)
    if (abs(uc) < 1e-14) { if (r > 0) return(c(NA, NA)); next }
    if (uc > 0) lo <- max(lo, r / uc) else hi <- min(hi, r / uc)
  }
  c(lo, hi)
}

llb <- function(X, s) ns$.emc_burn_ll_batch(X, ch0$data[[s]], ch0$model, pn)
ll_rows <- function(rows) {  # rows: list of per-subject matrices
  unlist(parallel::mclapply(seq_len(n), function(s) list(llb(rows[[s]], s)), mc.cores = NC), recursive = FALSE)
}

# ---- piecewise line density -------------------------------------------------
# Pieces per subject: [L, R] with log Ltilde = a + b t. Built from node
# values; floor pieces beyond bounds take the floor value measured there.
build_pieces <- function(tn, y, lo, hi, yfl_lo, yfl_hi) {
  Kn <- length(tn)
  b <- diff(y) / diff(tn); a0 <- y[-Kn] - b * tn[-Kn]
  L <- tn[-Kn]; R <- tn[-1]
  bl <- if (is.finite(lo) && abs(tn[1] - lo) < 1e-12) NA else max(b[1], 0)
  br <- if (is.finite(hi) && abs(tn[Kn] - hi) < 1e-12) NA else min(b[Kn - 1], 0)
  left <- if (is.na(bl)) c(-Inf, tn[1], yfl_lo, 0) else c(-Inf, tn[1], y[1] - bl * tn[1], bl)
  right <- if (is.na(br)) c(tn[Kn], Inf, yfl_hi, 0) else c(tn[Kn], Inf, y[Kn] - br * tn[Kn], br)
  P <- rbind(left, cbind(L, R, a0, b), right)
  P[is.finite(P[, 3]) | TRUE, , drop = FALSE]
}
ldphi <- function(zr, zl) {  # log(Phi(zr) - Phi(zl)), zr >= zl
  out <- rep(-Inf, length(zr))
  up <- zl >= 0
  i <- up & zr > zl
  if (any(i)) { s1 <- pnorm(zl[i], lower.tail = FALSE, log.p = TRUE); s2 <- pnorm(zr[i], lower.tail = FALSE, log.p = TRUE)
    out[i] <- s1 + log1p(-exp(s2 - s1)) }
  i <- !up & zr > zl
  if (any(i)) { f2 <- pnorm(zr[i], log.p = TRUE); f1 <- pnorm(zl[i], log.p = TRUE)
    out[i] <- f2 + log1p(-exp(f1 - f2)) }
  out
}
# Line quadratic of the group density: log N(w + t c) = C - A t^2 / 2 + B t.
line_quad <- function(w, c_s, st, Vi, ldV) {
  r <- w - st$mu
  A <- sum(c_s * (Vi %*% c_s)); B <- -sum(c_s * (Vi %*% r))
  C <- -0.5 * sum(r * (Vi %*% r)) - 0.5 * ldV - 0.5 * p * log(2 * pi)
  c(A, B, C)
}
piece_logmass <- function(P, q) {
  A <- q[1]; B <- q[2]; C <- q[3]
  m <- (B + P[, 4]) / A; sA <- sqrt(A)
  C + P[, 3] + (B + P[, 4])^2 / (2 * A) + 0.5 * log(2 * pi / A) + ldphi((P[, 2] - m) * sA, (P[, 1] - m) * sA)
}
lse <- function(x) { mx <- max(x); if (!is.finite(mx)) return(mx); mx + log(sum(exp(x - mx))) }
log_ltilde <- function(P, t) { i <- which(t >= P[, 1] & t <= P[, 2])[1]; P[i, 3] + P[i, 4] * t }
# Forward CDF as (side, log mass on that side) and the inverse.
cdf_side <- function(P, q, t) {
  lm <- piece_logmass(P, q); lz <- lse(lm)
  i <- which(t >= P[, 1] & t <= P[, 2])[1]
  A <- q[1]; m <- (q[2] + P[i, 4]) / A; sA <- sqrt(A)
  lleft <- lse(c(lm[seq_len(i - 1)], piece_logmass(P[i, , drop = FALSE], q) - ldphi((P[i, 2] - m) * sA, (P[i, 1] - m) * sA) + ldphi((t - m) * sA, (P[i, 1] - m) * sA)))
  lright <- lse(c(lm[-seq_len(i)], piece_logmass(P[i, , drop = FALSE], q) - ldphi((P[i, 2] - m) * sA, (P[i, 1] - m) * sA) + ldphi((P[i, 2] - m) * sA, (t - m) * sA)))
  if (lleft <= lright) list(side = "L", lf = lleft - lz) else list(side = "R", lf = lright - lz)
}
inv_cdf <- function(P, q, side, lf) {
  lm <- piece_logmass(P, q); lz <- lse(lm); target <- lf + lz
  ord <- if (side == "L") seq_len(nrow(P)) else rev(seq_len(nrow(P)))
  acc <- -Inf
  for (i in ord) {
    nxt <- lse(c(acc, lm[i]))
    if (nxt >= target || i == ord[length(ord)]) {
      rem <- if (is.finite(acc)) target + log1p(-min(1, exp(acc - target))) else target
      frac <- min(rem - lm[i], 0)
      A <- q[1]; m <- (q[2] + P[i, 4]) / A; sA <- sqrt(A)
      zl <- (P[i, 1] - m) * sA; zr <- (P[i, 2] - m) * sA; tot <- ldphi(zr, zl)
      g <- if (side == "L") function(z) ldphi(z, zl) - tot - frac else function(z) frac - (ldphi(zr, z) - tot)
      z <- uniroot(g, c(max(zl, -38), min(zr, 38)), tol = 1e-14)$root
      return(m + z / sA)
    }
    acc <- nxt
  }
}
log_q <- function(P, q, t) {
  lz <- lse(piece_logmass(P, q))
  log_ltilde(P, t) + q[3] - 0.5 * q[1] * t^2 + q[2] * t - lz
}

# ---- one draw -------------------------------------------------------------
run_draw <- function(d) {
  alpha <- matrix(ch0$samples$alpha[, , d], p)
  mu <- as.numeric(ch0$samples$theta_mu[, d]); V <- as.matrix(ch0$samples$theta_var[, , d])
  t_cur <- colSums(Ell * alpha); W <- alpha - sweep(Cl, 2, t_cur, `*`)
  # nodes
  fe <- t(vapply(seq_len(n), function(s) feasible(s, W[, s]), numeric(2)))
  nodes <- lapply(seq_len(n), function(s) {
    lo <- tc[s] - SPAN * tsd[s]; hi <- tc[s] + SPAN * tsd[s]
    lo2 <- max(lo, fe[s, 1]); hi2 <- min(hi, fe[s, 2])
    if (hi2 <= lo2) { lo2 <- fe[s, 1]; hi2 <- fe[s, 1] + 2 * SPAN * tsd[s] }
    tn <- if (SPACING == "uniform") seq(lo2, hi2, length.out = K) else {
      g <- qnorm(seq(0.5 / K, 1 - 0.5 / K, length.out = K)); g <- (g - g[1]) / (g[K] - g[1]); lo2 + g * (hi2 - lo2) }
    extra <- c(if (is.finite(fe[s, 1])) fe[s, 1] - tsd[s], if (is.finite(fe[s, 2])) fe[s, 2] + tsd[s])
    list(tn = tn, extra = extra, lo = fe[s, 1], hi = fe[s, 2])
  })
  rows <- lapply(seq_len(n), function(s) { tt <- c(nodes[[s]]$tn, nodes[[s]]$extra)
    X <- t(W[, s] + outer(Cl[, s], tt)); colnames(X) <- pn; X })
  t0 <- proc.time()[3]; ylist <- ll_rows(rows); el_nodes <- proc.time()[3] - t0
  cur_ll <- ll_rows(lapply(seq_len(n), function(s) { X <- matrix(alpha[, s], 1); colnames(X) <- pn; X }))
  cur_ll <- unlist(cur_ll)
  Pl <- lapply(seq_len(n), function(s) {
    nd <- nodes[[s]]; y <- ylist[[s]]; Kn <- length(nd$tn)
    yfl_lo <- if (is.finite(nd$lo)) y[Kn + 1] else NA
    yfl_hi <- if (is.finite(nd$hi)) y[Kn + 1 + is.finite(nd$lo)] else NA
    yy <- y[seq_len(Kn)]; yy[!is.finite(yy)] <- min(yy[is.finite(yy)], -1e6)
    build_pieces(nd$tn, yy, nd$lo, nd$hi, yfl_lo, yfl_hi)
  })
  quads <- function(st) { Vi <- solve(st$V); ldV <- determinant(st$V)$modulus
    lapply(seq_len(n), function(s) line_quad(W[, s], Cl[, s], st, Vi, ldV)) }
  lpt <- function(u) {
    st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g)) return(-Inf)
    qs <- quads(st)
    g + jac_u(u) + sum(vapply(seq_len(n), function(s) lse(piece_logmass(Pl[[s]], qs[[s]])), 1))
  }
  # residual of the current point
  res_cur <- cur_ll - vapply(seq_len(n), function(s) log_ltilde(Pl[[s]], t_cur[s]), 1)
  # mobility: pi~ on a grid in window sds
  ug <- seq(-4, 4, by = 0.25) * wsd
  lg <- vapply(ug, lpt, 1); wg <- exp(lg - max(lg)); wg <- wg / sum(wg)
  mob_mean <- sum(wg * ug) / wsd; mob_sd <- sqrt(sum(wg * ug^2) - sum(wg * ug)^2) / wsd
  # the same coordinate's conditional with every subject held (centred Gibbs)
  lfix <- vapply(ug, function(u) { st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g)) return(-Inf); ch <- chol(st$V); g + jac_u(u) + ns$.emc_shear_log_hier(alpha, list(mu = st$mu, chol_V = ch)) }, 1)
  wf <- exp(lfix - max(lfix)); wf <- wf / sum(wf); fix_sd <- sqrt(sum(wf * ug^2) - sum(wf * ug)^2) / wsd
  # every subject integrated out under its Gaussian frame (the bound for any subject-collapsing move)
  lgau <- vapply(ug, function(u) { st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g)) return(-Inf); Vi <- solve(st$V); ldV <- determinant(st$V)$modulus; Vm <- drop(Vi %*% st$mu)
    g + jac_u(u) + sum(vapply(seq_len(n), function(s) { Q <- H[, , s] + Vi; h <- LIN[, s] + Vm
      -0.5 * ldV - 0.5 * determinant(Q)$modulus + 0.5 * sum(h * solve(Q, h)) - 0.5 * sum(st$mu * Vm) }, 1)) }, 1)
  wq <- exp(lgau - max(lgau)); wq <- wq / sum(wq); gau_sd <- sqrt(sum(wq * ug^2) - sum(wq * ug)^2) / wsd
  # transport to u, residual
  q0 <- quads(state_at(mu, V, 0))
  transport <- function(u, from = 0, tf = t_cur) {
    q1 <- quads(state_at(mu, V, u)); qf <- if (from == 0) q0 else quads(state_at(mu, V, from))
    t1 <- numeric(n)
    for (s in seq_len(n)) {
      cs <- cdf_side(Pl[[s]], qf[[s]], tf[s])
      t1[s] <- inv_cdf(Pl[[s]], q1[[s]], cs$side, cs$lf)
    }
    t1
  }
  resid_at <- function(t1s) {
    A1 <- W + sweep(Cl, 2, t1s, `*`)
    ll1 <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(A1[, s], 1); colnames(X) <- pn; X })))
    r1 <- ll1 - vapply(seq_len(n), function(s) log_ltilde(Pl[[s]], t1s[s]), 1)
    r1 - res_cur
  }
  l0 <- lpt(0)
  tb <- transport(0, wsd, transport(wsd)); inv_err <- max(abs(tb - t_cur) / tsd)
  jumps <- c(-1, -0.5, 0.5, 1)
  jacc <- vapply(jumps, function(z) {
    u <- z * wsd; t1 <- transport(u); r <- resid_at(t1)
    min(1, exp(lpt(u) - l0 + sum(r))) }, 1)
  # slice proposals from pi~
  slice <- function(u0, lu0, w0) {
    y <- lu0 - rexp(1); L <- u0 - runif(1) * w0; R <- L + w0
    for (i in 1:10) { if (lpt(L) < y) break; L <- L - w0 }
    for (i in 1:10) { if (lpt(R) < y) break; R <- R + w0 }
    repeat { u1 <- runif(1, L, R); l1 <- lpt(u1); if (l1 >= y) return(c(u1, l1)); if (u1 < u0) L <- u1 else R <- u1 }
  }
  props <- t(vapply(seq_len(REPS), function(r) {
    u <- 0; lu <- l0
    for (m in seq_len(M_SLICE)) { o <- slice(u, lu, wsd); u <- o[1]; lu <- o[2] }
    t1 <- transport(u); rr <- resid_at(t1)
    c(u = u, acc = min(1, exp(sum(rr))), rsum = sum(rr), rsd = sd(rr), rmax = max(abs(rr)))
  }, numeric(5)))
  list(d = d, mob = c(mob_mean, mob_sd, fix_sd, gau_sd), jacc = jacc, props = props, el_nodes = el_nodes,
       res_cur_sd = sd(res_cur), res_cur_max = max(abs(res_cur)), inv_err = inv_err)
}

# ---- hybrid: transverse by the Gaussian frame's quotient map, line exact ---
# w = U y moves by the affine map between the frames' marginals of y; t by
# the exact line CDF from the line at w to the line at w'.
Ub <- lapply(seq_len(n), function(s) qr.Q(qr(cbind(Ell[, s], diag(p))))[, -1, drop = FALSE])
frame_of <- function(st) { Vi <- solve(st$V)
  lapply(seq_len(n), function(s) { Q <- H[, , s] + Vi; list(m = drop(solve(Q, LIN[, s] + Vi %*% st$mu)), S = solve(Q)) }) }
hybrid_jumps <- function(d, jumps) {
  alpha <- matrix(ch0$samples$alpha[, , d], p)
  mu <- as.numeric(ch0$samples$theta_mu[, d]); V <- as.matrix(ch0$samples$theta_var[, , d])
  t_cur <- colSums(Ell * alpha); W <- alpha - sweep(Cl, 2, t_cur, `*`)
  cur_ll <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(alpha[, s], 1); colnames(X) <- pn; X })))
  line_at <- function(Wm) {   # pieces per subject for lines through columns of Wm
    fe <- t(vapply(seq_len(n), function(s) feasible(s, Wm[, s]), numeric(2)))
    nodes <- lapply(seq_len(n), function(s) {
      lo2 <- max(tc[s] - SPAN * tsd[s], fe[s, 1]); hi2 <- min(tc[s] + SPAN * tsd[s], fe[s, 2])
      if (!is.finite(lo2) || !is.finite(hi2) || hi2 <= lo2) return(NULL)
      list(tn = seq(lo2, hi2, length.out = K), extra = c(if (is.finite(fe[s, 1])) fe[s, 1] - tsd[s], if (is.finite(fe[s, 2])) fe[s, 2] + tsd[s]), lo = fe[s, 1], hi = fe[s, 2]) })
    if (any(vapply(nodes, is.null, TRUE))) return(NULL)
    y <- ll_rows(lapply(seq_len(n), function(s) { X <- t(Wm[, s] + outer(Cl[, s], c(nodes[[s]]$tn, nodes[[s]]$extra))); colnames(X) <- pn; X }))
    lapply(seq_len(n), function(s) { nd <- nodes[[s]]; yy <- y[[s]]; Kn <- length(nd$tn)
      yfl_lo <- if (is.finite(nd$lo)) yy[Kn + 1] else NA; yfl_hi <- if (is.finite(nd$hi)) yy[Kn + 1 + is.finite(nd$lo)] else NA
      y0 <- yy[seq_len(Kn)]; y0[!is.finite(y0)] <- -1e6
      build_pieces(nd$tn, y0, nd$lo, nd$hi, yfl_lo, yfl_hi) })
  }
  P0 <- line_at(W)
  st0 <- state_at(mu, V, 0); f0 <- frame_of(st0)
  Vi0 <- solve(V); ld0 <- determinant(V)$modulus
  q0 <- lapply(seq_len(n), function(s) line_quad(W[, s], Cl[, s], st0, Vi0, ld0))
  lz0 <- vapply(seq_len(n), function(s) lse(piece_logmass(P0[[s]], q0[[s]])), 1)
  r0 <- cur_ll - vapply(seq_len(n), function(s) log_ltilde(P0[[s]], t_cur[s]), 1)
  lp0 <- lp_group(st0)
  vapply(jumps, function(z) {
    u <- z * wsd; st1 <- state_at(mu, V, u); lp1 <- tryCatch(lp_group(st1), error = function(e) -Inf)
    if (!is.finite(lp1)) return(0)
    f1 <- frame_of(st1); W1 <- W; lA <- numeric(n)
    for (s in seq_len(n)) {
      M <- crossprod(Ub[[s]], diag(p) - Cl[, s] %*% t(Ell[, s]))
      L0 <- t(chol(M %*% f0[[s]]$S %*% t(M))); L1 <- t(chol(M %*% f1[[s]]$S %*% t(M)))
      y0 <- drop(crossprod(Ub[[s]], W[, s]))
      y1 <- drop(M %*% f1[[s]]$m) + drop(L1 %*% forwardsolve(L0, y0 - drop(M %*% f0[[s]]$m)))
      W1[, s] <- drop(Ub[[s]] %*% y1); lA[s] <- sum(log(diag(L1))) - sum(log(diag(L0)))
    }
    P1 <- line_at(W1); if (is.null(P1)) return(0)
    Vi1 <- solve(st1$V); ld1 <- determinant(st1$V)$modulus
    q1 <- lapply(seq_len(n), function(s) line_quad(W1[, s], Cl[, s], st1, Vi1, ld1))
    lz1 <- vapply(seq_len(n), function(s) lse(piece_logmass(P1[[s]], q1[[s]])), 1)
    t1 <- vapply(seq_len(n), function(s) { cs <- cdf_side(P0[[s]], q0[[s]], t_cur[s]); inv_cdf(P1[[s]], q1[[s]], cs$side, cs$lf) }, 1)
    A1 <- W1 + sweep(Cl, 2, t1, `*`)
    ll1 <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(A1[, s], 1); colnames(X) <- pn; X })))
    r1 <- ll1 - vapply(seq_len(n), function(s) log_ltilde(P1[[s]], t1[s]), 1)
    lr <- lp1 + jac_u(u) - lp0 + sum(lz1 - lz0) + sum(lA) + sum(r1 - r0)
    min(1, exp(lr)) }, 1)
}
frame_jumps <- function(d, jumps) {
  alpha <- matrix(ch0$samples$alpha[, , d], p)
  mu <- as.numeric(ch0$samples$theta_mu[, d]); V <- as.matrix(ch0$samples$theta_var[, , d])
  cur <- ns$.emc_shear_state(mu, V, ctx$packed, prior, spec); cur_hier <- ns$.emc_shear_log_hier(alpha, cur)
  cur_ll <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(alpha[, s], 1); colnames(X) <- pn; X })))
  v <- ns$.emc_shear_coord_dir(plan, kk); tt <- jumps * wsd / (if (kind == "lsd") 1 else 1)
  pr <- ns$.emc_shear_candidates(cur, v, tt, alpha, ctx$packed, prior, spec, ctx$walls)
  vapply(seq_along(tt), function(i) { z <- pr[[i]]; if (is.null(z)) return(0)
    ll <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(z$alpha[, s], 1); colnames(X) <- pn; X })))
    min(1, exp(sum(ll) - sum(cur_ll) + z$hier + z$state$lp + z$log_jac - cur_hier - cur$lp)) }, 1)
}

# ---- fixed-frame shear rounds on the same draws -----------------------------
shear_draw <- function(d) {
  alpha <- matrix(ch0$samples$alpha[, , d], p)
  mu <- as.numeric(ch0$samples$theta_mu[, d]); V <- as.matrix(ch0$samples$theta_var[, , d])
  cur <- ns$.emc_shear_state(mu, V, ctx$packed, prior, spec)
  cur_hier <- ns$.emc_shear_log_hier(alpha, cur)
  cur_ll <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(alpha[, s], 1); colnames(X) <- pn; X })))
  v <- ns$.emc_shear_coord_dir(plan, kk)
  step <- exp(plan$log_step[kk])
  E <- numeric(REPS); el <- numeric(REPS)
  for (r in seq_len(REPS)) {
    ta <- step * rnorm(1); t <- ta + step * rnorm(ns$.EMC_SHEAR_CANDIDATES)
    pr <- ns$.emc_shear_candidates(cur, v, t, alpha, ctx$packed, prior, spec, ctx$walls)
    ok <- which(!vapply(pr, is.null, TRUE))
    rows <- lapply(seq_len(n), function(s) { X <- t(vapply(pr[ok], function(z) z$alpha[, s], numeric(p))); colnames(X) <- pn; X })
    t0 <- proc.time()[3]; L <- ll_rows(rows); el[r] <- proc.time()[3] - t0
    Lm <- matrix(unlist(L), length(ok))
    lw <- rep(-Inf, length(t))
    lw[ok] <- colSums(t(Lm)) - sum(cur_ll) + vapply(pr[ok], function(z) z$hier + z$state$lp + z$log_jac - cur_hier - cur$lp, 1)
    w <- c(0, lw); w <- exp(w - max(w)); w <- w / sum(w)
    # E in units of the coordinate: the step u = t v moves the coordinate by t
    # (log sd) or by t (mean direction, dirs column scale)
    E[r] <- sum(w[-1] * t^2)
  }
  list(E = mean(E), el = mean(el))
}

cat(sprintf("%s: chain %d, %s; window sd %.4f; K=%d span=%g dir=%s spacing=%s; n=%d p=%d\n",
            basename(a[1]), j, label, wsd, K, SPAN, DIR, SPACING, n, p))
if (kv("hybrid", "0") == "1") {
  jumps <- c(-1, -0.5, 0.5, 1); HA <- NULL; FA <- NULL
  for (d in draws) { h <- hybrid_jumps(d, jumps); fr <- frame_jumps(d, jumps); HA <- rbind(HA, h); FA <- rbind(FA, fr)
    cat(sprintf("draw %4d jump acc (-1 -.5 +.5 +1): hybrid %s | fixed frames %s\n", d, paste(sprintf("%.2f", h), collapse = " "), paste(sprintf("%.2f", fr), collapse = " "))) }
  cat("mean hybrid:", round(colMeans(HA), 2), " fixed:", round(colMeans(FA), 2), "\n"); quit(save = "no")
}
cat(sprintf("line span (sd_t) median %.3g; walls on %d subjects\n", median(tsd), sum(lengths(support) > 0)))
tot <- NULL
for (d in draws) {
  t0 <- proc.time()[3]
  o <- run_draw(d); sh <- shear_draw(d)
  pr <- o$props
  esjd_line <- mean(pr[, "acc"] * pr[, "u"]^2) / wsd^2
  esjd_shear <- sh$E / wsd^2
  cat(sprintf("draw %4d: mobility mean %+.2f sd %.2f (held subjects %.2f, all collapsed Gaussian %.2f) | jump acc -1 %.2f -.5 %.2f +.5 %.2f +1 %.2f | slice props: acc mean %.2f med %.2f, resid sum sd %.2f, per-subj sd %.3f max %.2f | cur resid sd %.3f max %.2f | ESJD/wsd^2 line %.3f shear %.3f | ll call el nodes %.2fs shear8 %.2fs | inv err %.1e (%.0fs)\n",
    d, o$mob[1], o$mob[2], o$mob[3], o$mob[4], o$jacc[1], o$jacc[2], o$jacc[3], o$jacc[4], mean(pr[, "acc"]), median(pr[, "acc"]),
    sd(pr[, "rsum"]), mean(pr[, "rsd"]), max(pr[, "rmax"]), o$res_cur_sd, o$res_cur_max, esjd_line, esjd_shear, o$el_nodes, sh$el, o$inv_err, proc.time()[3] - t0))
  tot <- rbind(tot, c(o$mob[2], o$mob[3], o$mob[4], o$jacc, mean(pr[, "acc"]), esjd_line, esjd_shear))
}
colnames(tot) <- c("mob_sd", "held_sd", "gauss_sd", "j-1", "j-.5", "j+.5", "j+1", "acc", "esjd_line", "esjd_shear")
cat("\nmeans over draws:\n"); print(round(colMeans(tot), 3))
