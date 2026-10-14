# Emulation of a collapsed group update over conditional importance-sampling
# clouds (one-step particle Gibbs).
#
# Extended target: each subject s holds M particles, one of them (index k_s)
# its value, the others drawn from a fixed proposal g_s:
#   pi(phi) prod_s [ L(a_s,k) N(a_s,k | phi) / M  prod_{i != k} g_s(a_s,i) ].
# Summing over k, phi given the clouds has density
#   pi(phi) prod_s (1/M) sum_i L(a_si) N(a_si | phi) / g_s(a_si),
# so a round is: fresh particles for every subject (one likelihood call),
# any pi~-reversible update of phi (here slice steps along one coordinate),
# then k_s ~ L N(. | phi') / g_s. Every step is an exact Gibbs step; there is
# no correction term. g_s is fixed per block (built from earlier draws only).
#
# Rscript is_collapse.R <emc.rds> <chain> lsd:<par>|mean:<k>|mean:slow <draws|auto>
#   [M=16,32,64] [g=gauss|t5] [infl=1.2] [rounds=20] [NC=8] [m=5]
suppressPackageStartupMessages(library(EMC2, lib.loc = if (nzchar(Sys.getenv("EMC_LIB"))) Sys.getenv("EMC_LIB") else NULL))
ns <- asNamespace("EMC2")
a <- commandArgs(TRUE)
kv <- function(key, def) { m <- grep(paste0("^", key, "="), a, value = TRUE); if (length(m)) sub(".*=", "", m[1]) else def }
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(a[1])), quiet = TRUE)
j <- as.integer(a[2]); coord <- a[3]
Ms <- as.integer(strsplit(kv("M", "16,32,64"), ",")[[1]]); GT <- kv("g", "gauss"); INFL <- as.numeric(kv("infl", 1.2))
ROUNDS <- as.integer(kv("rounds", 20)); NC <- as.integer(kv("NC", 8)); M_SLICE <- as.integer(kv("m", 5))
set.seed(as.integer(kv("seed", 1)))
ch0 <- emc[[j]]; p <- ch0$n_pars; n <- ch0$n_subjects; pn <- ch0$par_names; last <- ch0$samples$idx
draws <- if (a[4] == "auto") round(seq(last - 7 * 60, last, by = 60)) else as.integer(strsplit(a[4], ",")[[1]])
idx <- ns$.emc_chain_proposal_window(emc, seq_len(last), "sample"); idx <- idx[idx <= last]
spec <- ns$.emc_shear_v_spec(ch0); prior <- ch0$prior
plan <- ch0$shear$plan
kind <- sub(":.*", "", coord); what <- sub(".*:", "", coord)
if (kind == "lsd") {
  l <- if (what %in% pn) match(what, pn) else as.integer(what)
  wsd <- sd(0.5 * log(ch0$samples$theta_var[l, l, idx]))
  state_at <- function(mu, V, u) { d <- rep(1, p); d[l] <- exp(u); list(mu = mu, V = V * outer(d, d)) }
  jac_u <- function(u) (spec$size[l] + 1) * u; label <- sprintf("log sd(%s)", pn[l])
} else {
  kk <- if (what == "slow") which.max(plan$lam) else as.integer(what); dvec <- plan$dirs[, kk]
  wsd <- sd(drop(crossprod(dvec, matrix(ch0$samples$theta_mu[, idx], p))) / sum(dvec^2))
  state_at <- function(mu, V, u) list(mu = mu + u * dvec, V = V); jac_u <- function(u) 0
  label <- sprintf("mean direction %d (lambda %.3f)", kk, plan$lam[kk])
}
lp_group <- function(st) ns$.emc_shear_log_prior_mu(st$mu, prior) + ns$.emc_shear_log_prior_v(st$V, spec)
# proposals from draws before the first emulated draw
gidx <- idx[idx < min(draws) - 1L]
G <- lapply(seq_len(n), function(s) { X <- t(matrix(ch0$samples$alpha[, s, gidx], p))
  S <- cov(X) * INFL^2; list(m = colMeans(X), R = chol(S)) })
DF <- if (GT == "t5") 5 else Inf
rg <- function(s, k) { z <- matrix(rnorm(k * p), k); if (is.finite(DF)) z <- z / sqrt(rchisq(k, DF) / DF)
  sweep(z %*% G[[s]]$R, 2, G[[s]]$m, `+`) }
lg <- function(s, X) { z <- backsolve(G[[s]]$R, t(X) - G[[s]]$m, transpose = TRUE); q <- colSums(z^2)
  ld <- -sum(log(diag(G[[s]]$R)))
  if (is.finite(DF)) ld + lgamma((DF + p) / 2) - lgamma(DF / 2) - p / 2 * log(DF * pi) - (DF + p) / 2 * log1p(q / DF)
  else ld - p / 2 * log(2 * pi) - q / 2 }
llb <- function(X, s) ns$.emc_burn_ll_batch(X, ch0$data[[s]], ch0$model, pn)
ll_rows <- function(rows) unlist(parallel::mclapply(seq_len(n), function(s) list(llb(rows[[s]], s)), mc.cores = NC), recursive = FALSE)
lse_cols <- function(X) { m <- apply(X, 2, max); m + log(colSums(exp(sweep(X, 2, m)))) }
# per-particle log N(a | state) for an n-list of M x p particle matrices -> M x n
logN <- function(P, st) { R <- chol(st$V); ld <- -sum(log(diag(R))) - p / 2 * log(2 * pi)
  vapply(P, function(X) { z <- backsolve(R, t(X) - st$mu, transpose = TRUE); ld - colSums(z^2) / 2 }, numeric(nrow(P[[1]]))) }
ug <- seq(-4, 4, by = 0.25) * wsd
sd_of <- function(lg) { w <- exp(lg - max(lg)); w <- w / sum(w); sqrt(sum(w * ug^2) - sum(w * ug)^2) / wsd }

cat(sprintf("%s: chain %d, %s; window sd %.4f; n=%d p=%d; g=%s infl %.2f from %d draws\n",
            basename(a[1]), j, label, wsd, n, p, GT, INFL, length(gidx)))
for (d in draws) {
  alpha <- matrix(ch0$samples$alpha[, , d], p); mu0 <- ch0$samples$theta_mu[, d]; V0 <- ch0$samples$theta_var[, , d]
  cur_ll <- unlist(ll_rows(lapply(seq_len(n), function(s) { X <- matrix(alpha[, s], 1); colnames(X) <- pn; X })))
  held <- sd_of(vapply(ug, function(u) { st <- state_at(mu0, V0, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
    if (!is.finite(g)) -Inf else g + jac_u(u) + sum(logN(lapply(seq_len(n), function(s) matrix(alpha[, s], 1)), st)) }, 1))
  out <- sprintf("draw %4d: held %.2f |", d, held)
  for (M in Ms) {
    # a chain of rounds from this draw; the first round's pi~ width also reported
    mu <- mu0; V <- V0; A <- alpha; LL <- cur_ll; us <- numeric(ROUNDS); w1 <- NA; ess <- NULL; el <- 0
    for (r in seq_len(ROUNDS)) {
      P <- lapply(seq_len(n), function(s) { X <- rbind(A[, s], rg(s, M - 1L)); colnames(X) <- pn; X })
      t0 <- proc.time()[3]; Y <- ll_rows(lapply(P, function(X) X[-1, , drop = FALSE])); el <- el + proc.time()[3] - t0
      Lm <- rbind(LL, vapply(Y, identity, numeric(M - 1)))         # M x n
      Gm <- vapply(seq_len(n), function(s) lg(s, P[[s]]), numeric(M))
      Lm[!is.finite(Lm)] <- -Inf
      base <- Lm - Gm
      lpt <- function(u) { st <- state_at(mu, V, u); g <- tryCatch(lp_group(st), error = function(e) -Inf)
        if (!is.finite(g)) return(-Inf); g + jac_u(u) + sum(lse_cols(base + logN(P, st))) }
      if (r == 1) w1 <- sd_of(vapply(ug, lpt, 1))
      u <- 0; lu <- lpt(0)
      for (m in seq_len(M_SLICE)) {
        y <- lu - rexp(1); L <- u - runif(1) * wsd; R <- L + wsd
        for (i in 1:10) { if (lpt(L) < y) break; L <- L - wsd }; for (i in 1:10) { if (lpt(R) < y) break; R <- R + wsd }
        repeat { u1 <- runif(1, L, R); l1 <- lpt(u1); if (l1 >= y) { u <- u1; lu <- l1; break }; if (u1 < u) L <- u1 else R <- u1 }
      }
      st <- state_at(mu, V, u); W <- base + logN(P, st)
      W <- exp(sweep(W, 2, apply(W, 2, max))); W <- sweep(W, 2, colSums(W), `/`)
      ess <- c(ess, 1 / colSums(W^2))
      k <- vapply(seq_len(n), function(s) sample.int(M, 1, prob = W[, s]), 1L)
      A <- vapply(seq_len(n), function(s) P[[s]][k[s], ], numeric(p)); LL <- Lm[cbind(k, seq_len(n))]
      mu <- st$mu; V <- st$V
      us[r] <- u
    }
    out <- paste0(out, sprintf(" M=%d: pi~ sd %.2f, ESJD/wsd^2 %.3f, net drift sd %.2f, ESS/subj med %.1f min %.1f, call %.2fs |",
      M, w1, mean(us^2) / wsd^2, sd(cumsum(us)) / wsd, median(ess), quantile(ess, 0.05), el / ROUNDS))
  }
  cat(out, "\n")
}
