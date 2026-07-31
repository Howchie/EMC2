# How much resolution does an extrapolated pair need at each alpha?
# Scores mean |d log pdf| per trial over the central 96% (the per-trial
# log-likelihood error the sampler pays) against a converged reference, on a
# pinned domain so the ladder really is a refinement ladder.  Output is the nx
# needed to hold a fixed error target as a function of alpha.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths()))  # RLFLIB: library to test
suppressMessages(library(EMC2))
pde <- EMC2:::rlf_fht_pdf_cdf_vec
solve1 <- function(tq, v, a, b0, z0, nx, L)
  pde(tq, v = v, sigma = 1, alpha = a, b0 = b0, z0 = z0,
      nx = nx, nt = 4L * nx, adaptive = FALSE, lower_extent = L)

alphas <- c(1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
cells  <- expand.grid(v = c(0.5, 1, 3), b0 = c(1, 2), KEEP.OUT.ATTRS = FALSE)
LAD <- c(48L, 64L, 72L, 96L, 128L, 144L, 192L, 256L, 288L)
NXREF <- 1024L
res <- NULL
for (a in alphas) for (i in seq_len(nrow(cells))) {
  g <- cells[i, ]; z0 <- 0.5 * g$b0
  L <- -solve1((g$b0 - z0) / g$v, g$v, a, g$b0, z0, 128L, -1)$x_lo
  tp <- exp(seq(log(0.02 * (g$b0 - z0) / g$v), log(60 * (g$b0 - z0) / g$v),
                length.out = 200))
  cp <- solve1(tp, g$v, a, g$b0, z0, 256L, L)$cdf
  tq <- approx(cp / max(cp), tp, xout = seq(0.02, 0.98, length.out = 25))$y
  tq <- tq[is.finite(tq)]
  P <- NULL; TM <- NULL
  for (nx in LAD) {
    tt <- system.time(r <- solve1(tq, g$v, a, g$b0, z0, nx, L))["elapsed"]
    P <- rbind(P, r$pdf); TM <- c(TM, as.numeric(tt))
  }
  rownames(P) <- names(TM) <- as.character(LAD)
  rf <- solve1(tq, g$v, a, g$b0, z0, NXREF, L)
  err <- function(pd) mean(abs(log(pmax(pd, 1e-300) / rf$pdf)))
  # p = 1 extrapolation of an (nx, 1.5 nx) pair, using measured dx ratios
  pairs <- list(c("48","72"), c("64","96"), c("96","144"), c("128","192"),
                c("192","288"))
  for (pr in pairs) {
    r <- 1.5
    ext <- (r * P[pr[2], ] - P[pr[1], ]) / (r - 1)
    res <- rbind(res, data.frame(alpha = a, v = g$v, b0 = g$b0,
      nx = as.integer(pr[1]), kind = "pair", err = err(ext),
      sec = TM[pr[1]] + TM[pr[2]]))
  }
  for (nx in c("64","96","128","192","288"))
    res <- rbind(res, data.frame(alpha = a, v = g$v, b0 = g$b0,
      nx = as.integer(nx), kind = "raw", err = err(P[nx, ]), sec = TM[nx]))
  cat(sprintf("alpha %.2f v %.1f b0 %.0f done\n", a, g$v, g$b0)); flush.console()
}
saveRDS(res, "nxlaw.rds")

agg <- aggregate(cbind(err, sec) ~ alpha + nx + kind, res, mean)
cat("\n=== mean |d log pdf| per trial, extrapolated pairs (nx, 1.5nx) ===\n")
w <- reshape(agg[agg$kind == "pair", c("alpha","nx","err")], idvar = "alpha",
             timevar = "nx", direction = "wide")
w <- w[order(w$alpha), ]
cat(sprintf("%-7s %9s %9s %9s %9s %9s\n", "alpha", "48+72", "64+96", "96+144",
            "128+192", "192+288"))
for (j in seq_len(nrow(w)))
  cat(sprintf("%-7.2f %9.5f %9.5f %9.5f %9.5f %9.5f\n", w$alpha[j],
      w$err.48[j], w$err.64[j], w$err.96[j], w$err.128[j], w$err.192[j]))

cat("\n=== raw solves, same metric ===\n")
wr <- reshape(agg[agg$kind == "raw", c("alpha","nx","err")], idvar = "alpha",
              timevar = "nx", direction = "wide")
wr <- wr[order(wr$alpha), ]
cat(sprintf("%-7s %9s %9s %9s %9s %9s\n", "alpha", "64", "96", "128", "192", "288"))
for (j in seq_len(nrow(wr)))
  cat(sprintf("%-7.2f %9.5f %9.5f %9.5f %9.5f %9.5f\n", wr$alpha[j],
      wr$err.64[j], wr$err.96[j], wr$err.128[j], wr$err.192[j], wr$err.288[j]))

cat("\n=== nx needed by an extrapolated pair to hold a target error ===\n")
cat(sprintf("%-7s %12s %12s %12s\n", "alpha", "e<=0.010", "e<=0.005", "e<=0.002"))
pa <- agg[agg$kind == "pair", ]
for (a in sort(unique(pa$alpha))) {
  s <- pa[pa$alpha == a, ]; s <- s[order(s$nx), ]
  need <- function(tol) {
    # log-log interpolate err(nx) and invert
    ok <- s$err > 0
    if (all(s$err[ok] <= tol)) return(min(s$nx))
    if (all(s$err[ok] > tol)) return(NA_real_)
    exp(approx(log(s$err[ok]), log(s$nx[ok]), xout = log(tol))$y)
  }
  cat(sprintf("%-7.2f %12.0f %12.0f %12.0f\n", a, need(0.01), need(0.005),
              need(0.002)))
}
cat("\nlawdone\n")
