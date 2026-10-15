d <- readRDS(file.path(commandArgs(TRUE)[1], "eval.rds"))
d$grp <- ifelse(d$set == "PM burn", "burn", "arms")
# coordinate-level slowness: mean over every block of the set, measured and realised
cs <- aggregate(cbind(v = (s25_meas + sR25) / 2) ~ set + coord, d, mean)
cs$slow <- ave(cs$v, cs$set, FUN = function(v) v >= sort(v, TRUE)[4])
d <- merge(d, cs[, c("set", "coord", "slow")])
sm <- function(x) { x <- x[is.finite(x)]; sprintf("%+.2f/%.2f", median(x), median(abs(x))) }
M <- c("L", "B", "K", "Kj", "A", "T", "S5", "S2", "PR")
for (g in c("burn", "arms")) for (slow in c(FALSE, TRUE)) {
  x <- d[d$grp == g & (!slow | d$slow), ]
  cat(sprintf("\n== %s %s n=%d  real Geyer tau med %.1f sR25 med %.1f | lag1 p-R1 %s\n", g, if (slow) "4 slowest coords/set" else "all", nrow(x),
    median(x$tR), median(x$sR25), sm(x$p - x$R1)))
  cat("  tau vs Geyer:", paste(sprintf("%s %s", M, sapply(M, function(m) sm(log(x[[paste0("t", m)]] / x$tR)))), collapse = "  "), "\n")
  M2 <- c("B", "T", "S5", "S2", "PR")
  cat("  sum<=25     :", paste(sprintf("%s %s", M2, sapply(M2, function(m) sm(log(x[[paste0("s25", m)]] / pmax(x$sR25, 1))))), collapse = "  "), "\n")
  cat("  acf10       :", paste(sprintf("%s %s", M2, sapply(M2, function(m) sm(x[[paste0("a10", m)]] - x$R10))), collapse = "  "), "\n")
  # response: does the predicted change from block b to the plan's block match the realised change?
  x$dr <- log(pmax(x$sR25, 1) / x$s25_measj)
  pr <- list(B = log((1 + 2 * (x$p / (1 - x$p)) * (1 - x$p^25)) / (1 + 2 * (x$r1_meas / (1 - x$r1_meas)) * (1 - x$r1_meas^25))),
    T = log(x$s25T / x$s25atT), S5 = log(x$s25S5 / x$s25atS5))
  cat("  spread of predicted change sd:", paste(sprintf("%s %.2f", names(pr), sapply(pr, sd, na.rm = TRUE)), collapse = " "),
      " | slope real~pred:", paste(sprintf("%s %.2f", names(pr), sapply(pr, function(z) coef(lm(x$dr ~ z))[2])), collapse = " "), "\n")
}
