suppressMessages(library(EMC2))
f <- commandArgs(TRUE)[1]; stg <- commandArgs(TRUE)[2]
load(f); st <- emc[[1]]$samples$stage; idx <- tail(which(st == stg), 500)
a <- lapply(emc, function(ch) ch$samples$alpha[, , idx, drop = FALSE])   # par x subj x iter
m <- simplify2array(lapply(a, function(x) apply(x, c(1, 2), mean)))     # par x subj x chain
s <- simplify2array(lapply(a, function(x) apply(x, c(1, 2), sd)))
gap <- apply(m, c(1, 2), function(z) diff(range(z))) / apply(s, c(1, 2), mean)
cat(f, stg, ": subject x par cells with chain gap > 3 within-SD:", sum(gap > 3), "of", length(gap), "\n")
w <- sort(apply(gap, 2, max), decreasing = TRUE)[1:8]
print(round(w, 1))
for (sj in names(w)[1:3]) { p <- names(which.max(gap[, sj])); cat(sj, p, "chain means", round(m[p, sj, ], 3), "within SD", round(mean(s[p, sj, ]), 3), "\n") }
