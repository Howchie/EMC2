#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EMC2)
})

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[[1]]) else 123
n <- if (length(args) >= 2) as.integer(args[[2]]) else 30000
out_png <- if (length(args) >= 3) args[[3]] else "WorkingTests/rdmgbm_density_cdf_check.png"

set.seed(seed)

# Sample sensible positive natural-scale parameters.
# Draw s first, then v so v/s > 0.5 with margin (positive GBM log-drift).
s_draw <- runif(1, 0.70, 1.40)
par_row <- c(
  v = runif(1, 1.05 * s_draw, 2.40 * s_draw),
  B = runif(1, 0.15, 1.25),
  A = runif(1, 0.05, 0.65),
  t0 = runif(1, 0.10, 0.35),
  s = s_draw,
  lambda_g = 0,
  lambda_k = 0
)

# Repeat one parameter set over trials for the model rfun.
pars <- matrix(rep(par_row, n), nrow = n, byrow = TRUE,
               dimnames = list(NULL, names(par_row)))
attr(pars, "ok") <- rep(TRUE, nrow(pars))

model <- RDMGBM()
lR <- factor(rep("acc", n), levels = "acc")
sim <- model$rfun(data = list(lR = lR), pars = pars)
rt <- sim$rt
rt <- rt[is.finite(rt)]

if (length(rt) < 100L) {
  stop("Too few finite RTs were generated; adjust parameter ranges.")
}

# Data-driven support for plotting.
q_lo <- as.numeric(quantile(rt, 0.001))
q_hi <- as.numeric(quantile(rt, 0.999))
t_lo <- max(min(rt), q_lo)
t_hi <- q_hi
if (!(t_hi > t_lo)) {
  t_lo <- min(rt)
  t_hi <- max(rt)
}
t_vec <- seq(t_lo, t_hi, length.out = 500)

pars_one <- matrix(par_row, nrow = 1, dimnames = list(NULL, names(par_row)))
pdf_vec <- model$dfun(t_vec, pars_one)
cdf_vec <- model$pfun(t_vec, pars_one)

# Numerical sanity checks.
pdf_area <- sum(diff(t_vec) * (head(pdf_vec, -1) + tail(pdf_vec, -1)) / 2)
ks <- suppressWarnings(ks.test(rt, function(x) {
  model$pfun(x, pars_one)
})$statistic)

# Plot.
png(out_png, width = 1500, height = 650, res = 120)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.2, 1.2))

hist(rt,
     breaks = "FD",
     freq = FALSE,
     border = "white",
     col = "grey85",
     main = "RDMGBM: Histogram vs Theoretical Density",
     xlab = "RT")
lines(t_vec, pdf_vec, col = "#0072B2", lwd = 2.5)
legend("topright", bty = "n", lwd = c(8, 2.5),
       col = c("grey85", "#0072B2"),
       legend = c("Simulated histogram", "Theoretical density"))

plot(ecdf(rt),
     main = "RDMGBM: ECDF vs Theoretical CDF",
     xlab = "RT", ylab = "F(t)",
     verticals = TRUE, do.points = FALSE,
     col = "grey35", lwd = 1.1)
lines(t_vec, cdf_vec, col = "#D55E00", lwd = 2.5)
legend("bottomright", bty = "n", lwd = c(1.1, 2.5),
       col = c("grey35", "#D55E00"),
       legend = c("Empirical CDF", "Theoretical CDF"))

mtext(sprintf("seed=%d | n=%d | pdf_area≈%.4f | KS≈%.4f", seed, length(rt), pdf_area, as.numeric(ks)),
      side = 1, line = -1.2, outer = FALSE, cex = 0.9)

dev.off()

cat("RDMGBM parameter draw (natural scale):\n")
print(round(par_row, 4))
cat(sprintf("Finite RTs: %d\n", length(rt)))
cat(sprintf("PDF integral over plotted range: %.6f\n", pdf_area))
cat(sprintf("KS statistic vs theoretical CDF: %.6f\n", as.numeric(ks)))
cat(sprintf("Plot written to: %s\n", out_png))
