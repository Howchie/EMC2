SP <- commandArgs(TRUE)[1]
d <- readRDS(file.path(SP, "eval.rds")); P <- readRDS(file.path(SP, "pairs.rds"))
fac <- function(h) ifelse(h$group_rho > 1e-9, h$pred_group / h$group_rho, NA)
get_F <- function(plans) do.call(rbind, lapply(seq_along(plans), function(j) do.call(rbind, lapply(plans[[j]], function(h)
  if (!is.null(h$pred_group)) data.frame(ch = j, R = h$iteration, coord = seq_along(h$pred_group), F = fac(h), rounds = h$rounds)))))
Fs <- rbind(cbind(set = "PM burn", get_F(P$pm_burn$plans)),
  do.call(rbind, lapply(names(P$arms), function(nm) cbind(set = paste("arm", nm), get_F(P$arms[[nm]]$plans)))))
Fo <- Fs; Fo$R <- Fo$R + ifelse(Fo$set == "PM burn", 100, ifelse(Fo$R == 1151, 100, 100)); names(Fo)[names(Fo) %in% c("F", "rounds")] <- c("F_old", "rounds_old")
x <- merge(merge(d, Fs), Fo)
x$dx <- log(x$F / x$F_old)
x <- x[is.finite(x$dx), ]
cat("rows with a previous plan:", nrow(x), " sd of log plan change", round(sd(x$dx), 2), "\n")
for (sub in list(all = rep(TRUE, nrow(x)), big = abs(x$dx) > 0.1)) {
  y <- x[sub, ]
  cat(sprintf(" n=%d: slope of log(pred/real Geyer) on log(F_new/F_old):", nrow(y)),
    paste(sapply(c("tB", "tK", "tA", "tT"), function(m) { z <- log(y[[m]] / y$tR); ok <- is.finite(z)
      sprintf("%s %+.2f (se %.2f)", m, coef(lm(z[ok] ~ y$dx[ok]))[2], summary(lm(z[ok] ~ y$dx[ok]))$coef[2, 2]) }), collapse = "  "), "\n")
  cat("   same with truncated sums lag<=25:", paste(sapply(c("B", "T", "S5"), function(m) { z <- log(y[[paste0("s25", m)]] / pmax(y$sR25, 1)); ok <- is.finite(z)
      sprintf("%s %+.2f", m, coef(lm(z[ok] ~ y$dx[ok]))[2]) }), collapse = "  "), "\n")
  # direct: realised change in log truncated sum vs plan change (chain-pooled measured sum at b)
  dr <- log(pmax(y$sR25, 1) / y$s25_meas); ok <- is.finite(dr)
  cat(sprintf("   realised log(sR25_{b+1}/s25_b) on plan change: slope %+.2f (se %.2f)\n",
    coef(lm(dr[ok] ~ y$dx[ok]))[2], summary(lm(dr[ok] ~ y$dx[ok]))$coef[2, 2]))
}
# what each model implies for that slope: d log tau / d log F at the measured point
y <- x; pp <- y$r1_pool
el_B <- 2 * pp / ((1 + pp) * (1 - pp))       # d log tau1 / d log rho
cat("implied elasticity d log tau / d log rho1, median: AR(1) ", round(median(el_B), 2),
    " A (additive excess) ", round(median(el_B * (1 + pp) / (1 - pp) / pmax(y$tau_meas, 1)), 2), "\n")
cat("\n--- response to the change in expected shear rounds (the plan's decision, no measurement noise)\n")
x$dn <- x$rounds - x$rounds_old
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-6); (1 + r) / (1 - r) }
exc <- pmax(0, x$tau_meas - tau1(x$r1_pool))
pc <- list(B = log(tau1(x$p) / tau1(x$r1_pool)), A = log((tau1(x$p) + exc) / (tau1(x$r1_pool) + exc)), T = log(x$s25T / x$s25atT))
real <- list(geyer = log(x$tR / x$tau_meas), s25 = log(pmax(x$sR25, 1) / x$s25_meas))
for (sub in c("all", "mu", "lsd", "slow half")) {
  ok <- switch(sub, all = TRUE, mu = x$kind == "mu", lsd = x$kind == "lsd", "slow half" = x$tau_meas >= median(x$tau_meas)) & abs(x$dn) > 0.01
  y <- x[ok, ]
  sl <- function(z) { k <- is.finite(z[ok]); b <- summary(lm(z[ok][k] ~ y$dn[k]))$coef[2, 1:2]; sprintf("%+.2f(%.2f)", b[1], b[2]) }
  cat(sprintf("%-9s n=%4d sd(dn) %.2f | slope on d rounds: real geyer %s real s25 %s | pred B/K %s  A %s  T(s25) %s\n", sub, nrow(y), sd(y$dn),
    sl(real$geyer), sl(real$s25), sl(pc$B), sl(pc$A), sl(pc$T)))
}
