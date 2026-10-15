# Check whether a joint Dirichlet surrogate assigns slow motion to the right kernels.
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
scratch <- file.path(root, ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/19c4c6df-25fb-4858-8ef0-fdc4d9b6acea/scratchpad")
outdir <- file.path(root, "plans/alloc_prediction_2026-10-15/eval_operator")
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", file.path(root, ".Rlib_alloc_review"))))
ns <- asNamespace("EMC2")
families <- c("particle", "group", "shear", "collapsed", "hmc")
blocks <- list()
path <- readRDS(file.path(root, "plans/alloc_prediction_2026-10-15/handoff_2026-10-15/arms/arm_path_2.rds"))
hs <- lapply(path$plans, function(h) h[!duplicated(vapply(h, `[[`, numeric(1), "iteration"), fromLast = TRUE)])
p <- ncol(path$alpha_mean[[1]])
pm_paths <- lapply(seq_along(hs), function(j)
  do.call(rbind, c(lapply(hs[[j]], `[[`, "path"), list(path$last_path[[j]]))))
for (b in seq_len(nrow(path$draws[[1]]) %/% 100L)) {
  rows <- (b - 1L) * 100L + seq_len(100L)
  P <- lapply(pm_paths, function(z) z[rows, , drop = FALSE])
  if (any(vapply(P, is.null, logical(1)))) next
  X <- Map(function(g, a) cbind(g[rows, ], a[rows, ]), path$draws, path$alpha_mean)
  if (all(vapply(P, nrow, integer(1)) == 100L)) blocks[[paste0("PM_", b)]] <- list(X = X, P = P, p = p)
}
for (name in c("full_T30", "full_forst")) {
  emc <- ns$restore_duplicates(readRDS(file.path(scratch, paste0(name, "_emc.rds"))))
  p <- emc[[1]]$n_pars
  hs <- lapply(emc, function(ch) ch$alloc$history[!duplicated(vapply(ch$alloc$history, `[[`, numeric(1), "iteration"), fromLast = TRUE)])
  its <- Reduce(intersect, lapply(hs, function(h) vapply(h, `[[`, numeric(1), "iteration")))
  for (it in its) {
    P <- lapply(hs, function(h) h[[which(vapply(h, `[[`, numeric(1), "iteration") == it)]]$path)
    if (any(vapply(P, is.null, logical(1)))) next
    n <- min(vapply(P, nrow, integer(1)))
    if (n < 50L || it > emc[[1]]$samples$idx) next
    rows <- (it - n + 1L):it
    state <- function(ch, k) t(vapply(k, function(j) c(ch$samples$theta_mu[, j],
      0.5 * log(diag(ch$samples$theta_var[, , j])), rowMeans(ch$samples$alpha[, , j])), numeric(3L * p)))
    X <- lapply(emc, state, k = rows)
    X0 <- lapply(emc, state, k = rows - 1L)
    P <- lapply(P, tail, n)
    if (any(vapply(seq_along(X), function(j) max(abs(Reduce(`+`, lapply(seq_along(families), function(f) P[[j]][, (f - 1L) * 3L * p + seq_len(3L * p)])) - (X[[j]] - X0[[j]]))) > 1e-8, logical(1)))) next
    blocks[[paste0(name, "_", emc[[1]]$samples$stage[it], "_", it)]] <- list(X = X, P = P, p = p)
  }
}
records <- list()
for (name in names(blocks)) {
  b <- blocks[[name]]
  X <- b$X; P <- b$P; p <- b$p; d <- 3L * p
  n <- nrow(X[[1]])
  sc <- sqrt(diag(var(do.call(rbind, X)))); sc[sc < 1e-12 | !is.finite(sc)] <- 1
  X <- lapply(X, function(x) sweep(x, 2L, sc, "/"))
  P <- lapply(P, function(x) sweep(x, 2L, rep(sc, length(families)), "/"))
  for (part in c("whole", "first_half", "last_half")) {
    rows <- switch(part, whole = seq_len(n), first_half = seq_len(n %/% 2L), last_half = (n %/% 2L + 1L):n)
    C <- var(do.call(rbind, lapply(X, function(x) x[rows, , drop = FALSE])))
    S <- lapply(seq_along(families), function(f) {
      a <- do.call(rbind, lapply(P, function(z) z[rows, (f - 1L) * d + seq_len(d), drop = FALSE]))
      crossprod(a) / nrow(a)
    })
    D <- Reduce(`+`, S)
    ed <- eigen(D, symmetric = TRUE)
    Di <- ed$vectors %*% (t(ed$vectors) / pmax(ed$values, 1e-8 * max(ed$values)))
    A <- Di %*% C
    taus <- 4 * colSums(C * A) / pmax(diag(C), 1e-12) - 1
    k <- which.max(taus[seq_len(2L * p)])
    ec <- eigen(C, symmetric = TRUE)
    keep <- ec$values > 1e-8 * max(ec$values)
    W <- sweep(ec$vectors[, keep, drop = FALSE], 2L, sqrt(ec$values[keep]), "/")
    ee <- eigen(crossprod(W, D %*% W), symmetric = TRUE)
    a <- drop(W %*% ee$vectors[, length(ee$values)])
    shares <- vapply(S, function(s) drop(crossprod(a, s %*% a)), numeric(1))
    shares <- shares / sum(shares)
    # Independent reference: transport in the selected coordinate over the
    # opposite half; no Markov surrogate is assumed for this measurement.
    rr <- switch(part, whole = seq_len(n), first_half = (n %/% 2L + 1L):n, last_half = seq_len(n %/% 2L))
    tf <- numeric(length(families))
    for (j in seq_along(P)) {
      dz <- sapply(seq_along(families), function(f) drop(P[[j]][rr, (f - 1L) * d + seq_len(d), drop = FALSE] %*% a))
      total <- rowSums(dz)
      L <- 10L
      smooth <- vapply(seq_along(total), function(t) {
        l <- -(L - 1L):(L - 1L); ok <- t + l >= 1L & t + l <= length(total)
        sum((1 - abs(l[ok]) / L) * total[t + l[ok]])
      }, numeric(1))
      tf <- tf + colMeans(dz * smooth) / length(P)
    }
    tf <- tf / sum(tf)
    records[[length(records) + 1L]] <- data.frame(block = name, part = part, coord = unname(k), tau = unname(taus[k]),
      family = families, predicted_share = unname(shares), measured_share = unname(tf))
  }
}
res <- do.call(rbind, records)
write.csv(res, file.path(outdir, "family_operator.csv"), row.names = FALSE)
for (dataset in c("PM", "full_T30_burn", "full_T30_sample", "full_forst")) {
  for (part in c("whole", "first_half", "last_half")) {
    x <- res[startsWith(res$block, dataset) & res$part == part, ]
    if (!nrow(x)) next
    cat(dataset, part, "\n")
    print(aggregate(cbind(predicted_share, measured_share) ~ family, x, median), row.names = FALSE, digits = 3)
  }
}
