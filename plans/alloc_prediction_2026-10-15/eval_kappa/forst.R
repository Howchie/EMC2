suppressMessages(library(EMC2, lib.loc = "/data/work/EMC2_dev_oo/.lib"))
fo <- EMC2:::restore_duplicates(readRDS("/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds"))
st <- fo[[1]]$samples$stage; cat("stages:", paste(names(table(st)), table(st)), " alloc history:", length(fo[[1]]$alloc$history), "\n")
idx <- which(st == "sample")
X <- lapply(fo, function(ch) { s <- ch$samples; cbind(t(s$theta_mu[, idx]), t(EMC2:::.emc_alloc_lsd_draws(s$theta_var, idx))) })
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-6); (1 + r) / (1 - r) }
acfp <- function(i) { pv <- apply(do.call(rbind, lapply(X, function(x) x[i, ])), 2, var)
  Reduce(`+`, lapply(X, function(x) { x <- x[i, ]; n <- nrow(x); sapply(1:50, function(k) 1 - colMeans((x[-(1:k), , drop = FALSE] - x[1:(n - k), , drop = FALSE])^2) / (2 * pv)) })) / length(X) }
geyer <- function(r, M = 50) { r <- r[1:M]; s <- 1; for (m in seq(1, M - 1, 2)) { pr <- r[m] + r[m + 1]; if (pr <= 0) break; s <- s + 2 * pr }; s }
out <- NULL; nb <- length(idx) %/% 100
for (b in 1:(nb - 1)) { a0 <- acfp((b - 1) * 100 + 1:100); a1 <- acfp(b * 100 + 1:100)
  for (c in 1:nrow(a0)) out <- rbind(out, data.frame(b = b, c = c, kap = geyer(a0[c, ]) / tau1(a0[c, 1]), tB = tau1(a0[c, 1]), tK = geyer(a0[c, ]), tR = geyer(a1[c, ], 30))) }
cat("Forstmann kappa50 quantiles:", round(quantile(out$kap, c(.1, .5, .9)), 2), " | log(pred/next) B", sprintf("%+.2f/%.2f", median(log(out$tB / out$tR)), median(abs(log(out$tB / out$tR)))),
    " K50g", sprintf("%+.2f/%.2f", median(log(out$tK / out$tR)), median(abs(log(out$tK / out$tR)))), "\n")
suppressMessages(library(posterior)); arr <- simplify2array(X); tr <- apply(arr, 2, function(a) prod(dim(a)) / ess_basic(a))
ag <- aggregate(cbind(tB, tK) ~ c, out, function(v) exp(mean(log(v)))); ag$run <- tr[ag$c]
cat("vs run-level tau: median log(pred/run) B", round(median(log(ag$tB / ag$run)), 2), " K50g", round(median(log(ag$tK / ag$run)), 2),
    " spearman B", round(cor(ag$tB, ag$run, method = "s"), 2), " K", round(cor(ag$tK, ag$run, method = "s"), 2), " run tau range", round(range(tr), 1), "\n")
