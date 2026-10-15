# Agreement attribution: each move family's movement towards the other chains, per block.
# Rscript agree.R pm <arm_path_2.rds> | full <emc.rds> <lib>
args <- commandArgs(TRUE); FAM <- c("particle", "group", "shear", "collapsed", "hmc"); nf <- 5
famcost <- function(h) c(particle = sum(h$particle_rate * h$c_part, na.rm = TRUE), group = h$C_rest,
  shear = sum(h$rounds) * h$c_round + h$shear_fixed, collapsed = h$collapsed_rounds * h$collapsed_cost,
  hmc = sum(h$rate * h$c_hmc, na.rm = TRUE) + h$hmc_fixed)
# X0: per chain n x d states at the start of each iteration; P: per chain n x (nf*d3) paths; cols: which of d3 path columns map to X0 columns.
score <- function(X0, P, d3, cols, half = NULL) {
  nch <- length(X0); n <- nrow(X0[[1]]); t <- if (is.null(half)) seq_len(n) else half
  V <- apply(do.call(rbind, X0), 2, var)
  A <- E <- matrix(0, nf, length(cols), dimnames = list(FAM, NULL)); gap2 <- 0
  net <- identical(Sys.getenv("NET"), "1")
  for (j in 1:nch) {
    gap <- Reduce(`+`, X0[-j]) / (nch - 1) - X0[[j]]
    if (net) gap <- matrix(colMeans(gap[t, , drop = FALSE]), nrow(gap), ncol(gap), byrow = TRUE)
    gap2 <- gap2 + colMeans(gap[t, , drop = FALSE]^2) / V / nch
    for (f in 1:nf) { D <- P[[j]][t, (f - 1) * d3 + cols, drop = FALSE]
      A[f, ] <- A[f, ] + colMeans(D * gap[t, , drop = FALSE]) / V / nch
      E[f, ] <- E[f, ] + colMeans(D^2) / (2 * V) / nch } }
  list(A = A, E = E, gap2 = gap2)
}
blocks <- list()
if (args[1] == "pm") {
  r <- readRDS(args[2]); ridge <- readRDS(file.path(dirname(args[2]), "ridge_u.rds"))
  r$plans <- lapply(r$plans, function(h) { it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] })
  p <- ncol(r$alpha_mean[[1]]); nch <- length(r$draws); d3 <- 3 * p
  Zst <- lapply(1:nch, function(j) cbind(r$draws[[j]], r$alpha_mean[[j]]))
  its <- vapply(r$plans[[1]], `[[`, numeric(1), "iteration"); n0 <- min(its) - 100
  paths <- lapply(1:nch, function(j) c(lapply(r$plans[[j]], `[[`, "path"), list(r$last_path[[j]])))
  ends <- c(its, n0 + nrow(Zst[[1]])); starts <- c(n0 + 1, head(ends, -1) + 1)
  for (b in seq_along(ends)) { rows <- (starts[b]:ends[b]) - n0; P <- lapply(paths, `[[`, b)
    if (any(vapply(P, is.null, logical(1))) || any(vapply(P, nrow, integer(1)) != length(rows))) next
    keep <- rows > 1; P <- lapply(P, function(x) x[keep, , drop = FALSE])
    X0 <- lapply(Zst, function(z) z[rows[keep] - 1, ])
    hp <- if (b > 1) r$plans[[1]][[b - 1]] else NULL
    blocks[[length(blocks) + 1]] <- list(X0 = X0, P = P, cost = if (is.null(hp)) rep(NA, nf) else famcost(hp), stage = "sample", R = ends[b]) }
  dirs <- list(ridge = c(ridge, numeric(p)))
} else {
  suppressMessages(library(EMC2, lib.loc = args[3])); ns <- asNamespace("EMC2")
  emc <- ns$restore_duplicates(readRDS(args[2])); nch <- length(emc); p <- emc[[1]]$n_pars; d3 <- 3 * p
  lsd <- function(V) 0.5 * log(pmax(diag(V), .Machine$double.xmin))
  st <- function(ch, k) t(sapply(k, function(t) c(ch$samples$theta_mu[, t], lsd(ch$samples$theta_var[, , t]),
    rowMeans(matrix(ch$samples$alpha[, , t], p)))))
  hist <- lapply(emc, function(ch) { h <- ch$alloc$history; it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] })
  its <- Reduce(intersect, lapply(hist, function(h) vapply(h, `[[`, numeric(1), "iteration")))
  for (R in its) { hs <- lapply(hist, function(h) Filter(function(x) x$iteration == R, h)[[1]])
    P <- lapply(hs, `[[`, "path"); if (any(vapply(P, is.null, logical(1)))) next
    n <- min(vapply(P, nrow, integer(1))); if (n < 50) next
    P <- lapply(P, function(x) tail(x, n)); rows <- (R - n + 1):R
    if (R > emc[[1]]$samples$idx || min(rows) < 2) next
    X0 <- lapply(emc, st, k = rows - 1L)
    prev <- Filter(function(x) x$iteration < R, hist[[1]])
    blocks[[length(blocks) + 1]] <- list(X0 = X0, P = P, cost = if (length(prev)) famcost(prev[[length(prev)]]) else rep(NA, nf),
      stage = emc[[1]]$samples$stage[R], R = R) }
  dirs <- list()
}
# Collective coordinates: half the group mean and half the subjects' mean, then the log sds.
M <- rbind(cbind(diag(p) / 2, matrix(0, p, p), diag(p) / 2), cbind(matrix(0, p, p), diag(p), matrix(0, p, p)))
blocks <- lapply(blocks, function(bk) { bk$X0 <- lapply(bk$X0, function(x) x %*% t(M))
  bk$P <- lapply(bk$P, function(x) do.call(cbind, lapply(1:nf, function(f) x[, (f - 1) * d3 + 1:d3] %*% t(M)))); bk })
d3 <- 2 * p; cols <- 1:(2 * p)
if (length(dirs)) dirs <- lapply(dirs, function(v) v[1:(2 * p)])
cat(sprintf("%s: %d blocks, %d coordinates\n", basename(args[2]), length(blocks), 2 * p))
res <- lapply(blocks, function(bk) { n <- nrow(bk$X0[[1]]); h1 <- 1:floor(n / 2); h2 <- (floor(n / 2) + 1):n
  c(list(all = score(bk$X0, bk$P, d3, cols), h1 = score(bk$X0, bk$P, d3, cols, h1), h2 = score(bk$X0, bk$P, d3, cols, h2)), bk[c("cost", "stage", "R")]) })
for (b in seq_along(res)) { x <- res[[b]]; A <- x$all$A; E <- x$all$E
  slow <- order(-x$all$gap2)[1:4]
  # Per family, over coordinates: half-block agreement of A, and A against the lag-one credit E.
  hc <- sapply(FAM, function(f) suppressWarnings(cor(x$h1$A[f, ], x$h2$A[f, ])))
  ae <- sapply(FAM, function(f) suppressWarnings(cor(A[f, ], E[f, ])))
  sh_gap <- rowSums(A[, slow, drop = FALSE]); sh_all <- rowSums(A)
  cat(sprintf("blk %2d %-6s R=%d gap2 med %.2f max %.2f | share on 4 most-apart coords: %s | all coords: %s\n   half-block cor: %s | cor(A,E): %s | per cost (most apart): %s\n",
    b, x$stage, x$R, median(x$all$gap2), max(x$all$gap2),
    paste(sprintf("%s %+.2f", substr(FAM, 1, 4), sh_gap / sum(abs(sh_gap))), collapse = " "),
    paste(sprintf("%+.2f", sh_all / sum(abs(sh_all))), collapse = " "),
    paste(sprintf("%.2f", hc), collapse = " "), paste(sprintf("%.2f", ae), collapse = " "),
    paste(sprintf("%.4f", sh_gap / x$cost), collapse = " ")))
}
if (length(res) > 1) { cat("consecutive-block correlation of A over coordinates, per family:\n")
  for (b in 2:length(res)) cat(sprintf("  %d->%d: %s\n", b - 1, b, paste(sprintf("%s %.2f", substr(FAM, 1, 4),
    sapply(FAM, function(f) suppressWarnings(cor(res[[b - 1]]$all$A[f, ], res[[b]]$all$A[f, ])))), collapse = " "))) }
if (length(dirs)) for (nm in names(dirs)) { v <- dirs[[nm]]
  for (b in seq_along(blocks)) { bk <- blocks[[b]]; Xv <- lapply(bk$X0, function(x) matrix(x %*% v))
    Pv <- lapply(bk$P, function(x) do.call(cbind, lapply(1:nf, function(f) x[, (f - 1) * d3 + cols] %*% v)))
    s <- score(Xv, Pv, 1, 1); a <- s$A[, 1]
    cat(sprintf(" %s blk %d gap2 %.2f: agreement share %s | per cost %s\n", nm, b, s$gap2, paste(sprintf("%s %+.2f", substr(FAM, 1, 4), a / sum(abs(a))), collapse = " "),
      paste(sprintf("%.4f", a / bk$cost), collapse = " "))) } }
