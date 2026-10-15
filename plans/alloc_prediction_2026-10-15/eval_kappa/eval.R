# One block predicts the next block under the plan built from it.
SP <- "SCRATCH"
P <- readRDS(file.path(SP, "pairs.rds")); p <- P$p
KMAX <- 50L
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-6); (1 + r) / (1 - r) }
# variogram ACF, lags 1..K, per chain (array chain x coord x lag) against the pooled variance
acfs <- function(X, idx, K = KMAX) { pv <- apply(do.call(rbind, lapply(X, function(x) x[idx, ])), 2, var)
  a <- array(NA, c(length(X), ncol(X[[1]]), K))
  for (j in seq_along(X)) { x <- X[[j]][idx, ]; n <- nrow(x)
    for (k in 1:K) a[j, , k] <- 1 - colMeans((x[-(1:k), , drop = FALSE] - x[1:(n - k), , drop = FALSE])^2) / (2 * pv) }
  a }
geyer <- function(r, M = 30, mono = FALSE) { r <- r[seq_len(min(M, length(r)))]; prev <- Inf; s <- 1; for (m in seq(1, length(r) - 1, 2)) { pr <- r[m] + r[m + 1]; if (pr <= 0) break; if (mono) { pr <- min(pr, prev); prev <- pr }; s <- s + 2 * pr }; s }
unlo <- function(lo, n) { c <- 2 * sqrt(2 / max(1, n)); ifelse(lo > 0, 1 - (1 - lo) / (1 + c), 1 - 1 / (1 + c)) }
# Fits of the measured block ACF r[1..K] (chain-averaged), each returning acf(k) after the plan sets lag-1 to p.
# slow-root models: acf(k) = (1-w) phi^k + w lam^k; the plan moves lam only.
plan_slow <- function(w, phi, p) { w <- min(max(w, 0), 1)
  if (w < 1e-6) return(list(w = 0, phi = p, lam = 0))
  lam <- min(max((p - (1 - w) * phi) / w, 0), 0.999); list(w = w, phi = phi, lam = lam) }
mix_acf <- function(m, k) (1 - m$w) * m$phi^k + m$w * m$lam^k
mix_tau <- function(m) (1 - m$w) * tau1(m$phi) + m$w * tau1(m$lam)
fit_tail <- function(r, lags = 2:10) { y <- r[lags]; ok <- y > 0.02
  if (sum(ok) < 3 || r[1] <= 0.02) return(NULL)
  cf <- coef(lm(log(y[ok]) ~ lags[ok])); lam <- min(exp(cf[2]), 0.999); w <- min(exp(cf[1]), 1)
  if (lam <= 0) return(NULL)
  phi <- if (w < 1) min(max((r[1] - w * lam) / (1 - w), 0), lam) else 0
  if (w * lam > r[1]) { w <- r[1] / lam; phi <- 0 }
  list(w = w, phi = phi, lam = lam) }
fit_two <- function(r, m = 5) { if (r[1] <= 0.02 || r[m] <= 0.02) return(NULL)
  lam <- min((r[m] / r[1])^(1 / (m - 1)), 0.999); w <- min(r[1] / lam, 1); list(w = w, phi = 0, lam = lam) }
fit_prony <- function(r) { A <- matrix(c(r[1], 1, r[2], r[1]), 2, byrow = TRUE)
  cc <- tryCatch(solve(A, c(r[2], r[3])), error = function(e) NULL); if (is.null(cc) || r[1] <= 0) return(NULL)
  d <- cc[1]^2 + 4 * cc[2]; if (!is.finite(d) || d < 0) return(NULL)
  l <- sort((cc[1] + c(-1, 1) * sqrt(d)) / 2); if (any(l <= -0.99) || any(l >= 1) || abs(diff(l)) < 1e-6) return(NULL)
  w <- (r[1] - l[1]) / (l[2] - l[1]); if (!is.finite(w) || w < 0 || w > 1) return(NULL); list(w = w, phi = l[1], lam = l[2]) }
run <- function(set, X, plans, scoreN = 100L, longR = NULL) {
  out <- list()
  its <- sort(unique(unlist(lapply(plans, function(h) vapply(h, `[[`, numeric(1), "iteration")))))
  for (R in its) {
    h <- lapply(plans, function(hh) { f <- Filter(function(x) x$iteration == R, hh); if (length(f)) f[[1]] })
    if (any(vapply(h, is.null, logical(1)))) next
    B <- h[[1]]$pairs[["group"]]; if (B < 100) next
    i0 <- (R - B + 1):R; i1 <- R + 1:(if (R %in% longR) 200L else scoreN)
    if (max(i1) > nrow(X[[1]]) || any(is.na(X[[1]][i0, 1]))) next
    a0 <- acfs(X, i0); a1 <- acfs(X, i1); r0 <- apply(a0, c(2, 3), mean)
    for (j in seq_along(X)) { hj <- h[[j]]; if (is.null(hj$pred_group)) next
      F <- ifelse(hj$group_rho > 1e-9, hj$pred_group / hj$group_rho, NA)
      rpt <- unlo(hj$group_rho, B)
      for (c in seq_len(ncol(X[[1]]))) {
        pl <- hj$pred_group[c]; pp <- if (is.na(F[c])) pl else min(rpt[c] * F[c], 0.999)
        rr <- r0[c, ]; rj <- a0[j, c, ]; real <- a1[j, c, ]
        kap <- geyer(rr) / tau1(rr[1])               # chain-pooled measured excess
        kapj <- geyer(rj) / tau1(rj[1])
        mods <- list(T = fit_tail(rr), S5 = fit_two(rr, 5), S2 = fit_two(rr, 2), PR = fit_prony(rr))
        row <- data.frame(set = set, R = R, ch = j, coord = c, kind = if (c <= ncol(X[[1]]) / 2) "mu" else "lsd",
          p_lo = pl, p = pp, r1_meas = rj[1], r1_pool = rr[1], tau_meas = geyer(rr),
          tR = geyer(real), R1 = real[1], R5 = real[5], R10 = real[10], R20 = real[20],
          sR10 = 1 + 2 * sum(real[1:10]), sR25 = 1 + 2 * sum(real[1:25]),
          tL = tau1(pl), tB = tau1(pp), tK = kap * tau1(pp), tKj = kapj * tau1(pp),
          tK10t = (1 + 2 * sum(rr[1:10])) / tau1(rr[1]) * tau1(pp), tK25t = (1 + 2 * sum(rr[1:25])) / tau1(rr[1]) * tau1(pp), tK50g = geyer(rr, 50) / tau1(rr[1]) * tau1(pp), tK50m = max(1, geyer(rr, 50, TRUE) / tau1(rr[1])) * tau1(pp), tK10g = geyer(rr, 10) / tau1(rr[1]) * tau1(pp), tK50gj = geyer(rj, 50) / tau1(rj[1]) * tau1(pp), tA = tau1(pp) + max(0, geyer(rr) - tau1(rr[1])), s25_meas = 1 + 2 * sum(rr[1:25]), s25_measj = 1 + 2 * sum(rj[1:25]), tau_measj = geyer(rj))
        for (nm in names(mods)) { m <- mods[[nm]]
          mp <- if (is.null(m)) list(w = 1, phi = 0, lam = pp) else plan_slow(m$w, m$phi, pp)
          ak <- mix_acf(mp, 1:25)
          row[[paste0("t", nm)]] <- mix_tau(mp); row[[paste0("s10", nm)]] <- 1 + 2 * sum(ak[1:10])
          row[[paste0("s25", nm)]] <- 1 + 2 * sum(ak); row[[paste0("a5", nm)]] <- ak[5]; row[[paste0("a10", nm)]] <- ak[10]
          row[[paste0("w", nm)]] <- mp$w
          m0 <- if (is.null(m)) list(w = 1, phi = 0, lam = rr[1]) else m; row[[paste0("s25at", nm)]] <- 1 + 2 * sum(mix_acf(m0, 1:25)) }
        ab <- pp^(1:25); row$s10B <- 1 + 2 * sum(ab[1:10]); row$s25B <- 1 + 2 * sum(ab); row$a5B <- ab[5]; row$a10B <- ab[10]
        out[[length(out) + 1]] <- row } } }
  do.call(rbind, out) }
d <- run("PM burn", P$pm_burn$X, P$pm_burn$plans)
for (nm in names(P$arms)) d <- rbind(d, run(paste("arm", nm), P$arms[[nm]]$X, P$arms[[nm]]$plans, longR = 1251))
saveRDS(d, file.path(SP, "eval.rds"))
