suppressMessages(library(posterior))
SP <- commandArgs(TRUE)[1]; d <- readRDS(file.path(SP, "eval.rds")); P <- readRDS(file.path(SP, "pairs.rds"))
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-6); (1 + r) / (1 - r) }
kap <- d$tK50g / d$tB
d$tKc <- pmax(1, kap) * d$tB; d$tA50 <- d$tB + pmax(0, kap - 1) * tau1(d$r1_pool)
d$grp <- ifelse(d$set == "PM burn", "burn", "arms")
cs <- aggregate(cbind(v = (s25_meas + sR25) / 2) ~ set + coord, d, mean)
cs$slow <- ave(cs$v, cs$set, FUN = function(v) v >= sort(v, TRUE)[4]); d <- merge(d, cs[, c("set", "coord", "slow")])
M <- c("B", "Kc", "K50m", "A", "A50")
sm <- function(x) { x <- x[is.finite(x)]; sprintf("%+.2f/%.2f", median(x), median(abs(x))) }
for (g in c("burn", "arms")) for (slow in c(FALSE, TRUE)) { x <- d[d$grp == g & (!slow | d$slow), ]
  cat(sprintf("%-4s %-6s next-block:", g, if (slow) "slow4" else "all"), paste(sprintf("%s %s", M, sapply(M, function(m) sm(log(x[[paste0("t", m)]] / x$tR)))), collapse = "  "), "\n") }
cat("run-level: median log(pred/run) all | slowest 6 | spearman\n")
for (nm in names(P$arms)) { X <- P$arms[[nm]]$X; arr <- simplify2array(lapply(X, function(x) x[1052:1451, ]))
  tr <- apply(arr, 2, function(a) 1200 / ess_basic(a)); x <- d[d$set == paste("arm", nm), ]
  ag <- aggregate(x[, paste0("t", M)], list(coord = x$coord), function(v) exp(mean(log(v)))); ag$run <- tr[ag$coord]; q <- ag$run >= sort(ag$run, TRUE)[6]
  cat(sprintf("%-8s", nm), paste(sprintf("%s %+.2f|%+.2f|%.2f", M, sapply(M, function(m) median(log(ag[[paste0("t", m)]] / ag$run))),
    sapply(M, function(m) median(log(ag[[paste0("t", m)]][q] / ag$run[q]))), sapply(M, function(m) cor(ag[[paste0("t", m)]], ag$run, method = "spearman"))), collapse = "  "), "\n") }
cat("share of rows with kappa < 1:", round(mean(kap < 1), 2), "\n")
