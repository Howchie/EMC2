# Block k predicts block k+1: ACF at lags 1,2,5,10 and tau, per coordinate and chain.
suppressMessages(library(EMC2, lib.loc = "/data/work/EMC2_dev_oo/.lib"))
SP <- Sys.getenv("SP"); B <- 100L; L <- c(1, 2, 5, 10)
tau1 <- function(r) { r <- pmin(r, 1 - 1e-6); (1 + r) / (1 - r) }
# lag-k autocorrelation per chain against the pooled-over-chains variance
lagr <- function(X, idx, k, pv) sapply(seq_along(X), function(j) { x <- X[[j]][idx, , drop = FALSE]; n <- nrow(x)
  1 - colMeans((x[-(1:k), , drop = FALSE] - x[1:(n - k), , drop = FALSE])^2) / (2 * pv) })
# two-rate (Prony) fit from lags 1-3; falls back to AR(1) on the measured lag-1
fit2 <- function(r1, r2, r3) {
  A <- matrix(c(r1, 1, r2, r1), 2, byrow = TRUE); cc <- tryCatch(solve(A, c(r2, r3)), error = function(e) NULL)
  ar1 <- list(l = c(r1, r1), w = 0)
  if (is.null(cc) || r1 <= 0) return(ar1)
  d <- cc[1]^2 + 4 * cc[2]; if (!is.finite(d) || d < 0) return(ar1)
  l <- sort((cc[1] + c(-1, 1) * sqrt(d)) / 2)
  if (any(l <= -0.99) || any(l >= 1) || abs(diff(l)) < 1e-6) return(ar1)
  w <- (r1 - l[1]) / (l[2] - l[1]); if (!is.finite(w) || w < 0 || w > 1) return(ar1)
  list(l = l, w = w)
}
acf2 <- function(f, k) (1 - f$w) * f$l[1]^k + f$w * f$l[2]^k
tau2 <- function(f) (1 - f$w) * tau1(f$l[1]) + f$w * tau1(f$l[2])
# realised tau of a block: 1 + 2 sum of pooled-variance ACF up to the first non-positive pair (Geyer)
taur <- function(X, idx, pv) { K <- 40; R <- sapply(1:K, function(k) rowMeans(lagr(X, idx, k, pv)))
  apply(R, 1, function(r) { s <- 1; for (m in seq(1, K - 1, 2)) { pr <- r[m] + r[m + 1]; if (pr <= 0) break; s <- s + 2 * pr }; s }) }
run <- function(name, X) {
  T <- nrow(X[[1]]); nb <- T %/% B; out <- NULL
  for (b in seq_len(nb - 1)) {
    i0 <- (b - 1) * B + 1:B; i1 <- b * B + 1:B
    pv0 <- apply(do.call(rbind, lapply(X, function(x) x[i0, ])), 2, var)
    pv1 <- apply(do.call(rbind, lapply(X, function(x) x[i1, ])), 2, var)
    r <- lapply(1:3, function(k) lagr(X, i0, k, pv0))
    real <- lapply(L, function(k) lagr(X, i1, k, pv1))
    tr <- taur(X, i1, pv1)
    for (j in seq_along(X)) for (c in seq_len(ncol(X[[1]]))) {
      r1 <- r[[1]][c, j]; m <- 1 - r1; lo <- max(0, 1 - m * (1 + 2 * sqrt(2 / (B - 1))))
      f <- fit2(r1, r[[2]][c, j], r[[3]][c, j])
      out <- rbind(out, data.frame(set = name, blk = b, ch = j, coord = c,
        A1 = lo, B1 = r1, C1 = acf2(f, 1), A2 = lo^2, B2 = r1^2, C2 = acf2(f, 2), A5 = lo^5, B5 = r1^5, C5 = acf2(f, 5),
        A10 = lo^10, B10 = r1^10, C10 = acf2(f, 10), tA = tau1(lo), tB = tau1(r1), tC = tau2(f),
        R1 = real[[1]][c, j], R2 = real[[2]][c, j], R5 = real[[3]][c, j], R10 = real[[4]][c, j], tR = tr[c]))
    }
  }
  out
}
group_X <- function(chs, idx, D) lapply(chs, function(ch) { s <- ch$samples
  cbind(t(D %*% s$theta_mu[, idx]), t(EMC2:::.emc_alloc_lsd_draws(s$theta_var, idx))) })
res <- list()
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; D <- solve(emc[[1]]$shear$plan$dirs)
st <- emc[[1]]$samples$stage
res$pm_burn <- run("PM burn", group_X(emc, which(st == "burn"), D))
for (f in c("old_1", "new_1", "slow_1")) { r <- readRDS(file.path(SP, "scr", paste0(f, ".rds")))
  res[[f]] <- run(paste("PM", f), lapply(r$draws, function(g) cbind(g[, 1:p] %*% t(D), g[, p + 1:p]))) }
fo <- EMC2:::restore_duplicates(readRDS("/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds"))
res$forst <- run("Forstmann", group_X(fo, 1:1000, diag(fo[[1]]$n_pars)))
d <- do.call(rbind, res); saveRDS(d, file.path(SP, "predtest.rds"))
sm <- function(x) sprintf("%+.3f/%.3f", median(x), median(abs(x)))
for (s in unique(d$set)) { x <- d[d$set == s, ]
  cat(sprintf("%-12s n=%4d | bias/MAE lag1 A %s B %s C %s | lag5 A %s B %s C %s | lag10 A %s B %s C %s\n", s, nrow(x),
    sm(x$A1 - x$R1), sm(x$B1 - x$R1), sm(x$C1 - x$R1), sm(x$A5 - x$R5), sm(x$B5 - x$R5), sm(x$C5 - x$R5),
    sm(x$A10 - x$R10), sm(x$B10 - x$R10), sm(x$C10 - x$R10)))
  cat(sprintf("%-12s        log tau ratio pred/real (median, MAE): A %s  B %s  C %s\n", "",
    sm(log(x$tA / x$tR)), sm(log(x$tB / x$tR)), sm(log(x$tC / x$tR))))
  # the coordinates that matter: slowest decile by realised tau
  q <- x$tR >= quantile(x$tR, 0.9)
  cat(sprintf("%-12s        slowest 10%%: real tau median %.1f, log ratio A %s B %s C %s\n", "", median(x$tR[q]),
    sm(log(x$tA[q] / x$tR[q])), sm(log(x$tB[q] / x$tR[q])), sm(log(x$tC[q] / x$tR[q]))))
}
