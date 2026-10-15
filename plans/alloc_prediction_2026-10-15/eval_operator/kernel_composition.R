# Estimate each reversible kernel's linear action, then compose in sampler order.
source("plans/alloc_prediction_2026-10-15/eval_operator/family_operator.R")
suppressPackageStartupMessages(library(posterior))
compose_fit <- function(X, P, p) {
  raw <- do.call(rbind, X)
  d <- ncol(raw)
  scale <- sqrt(diag(var(raw))); scale[!is.finite(scale) | scale <= 1e-12] <- 1
  Z <- lapply(X, function(x) sweep(x, 2L, scale, "/"))
  C <- var(do.call(rbind, Z))
  ec <- eigen(C, symmetric = TRUE)
  keep <- ec$values > 1e-8 * max(ec$values)
  U <- ec$vectors[, keep, drop = FALSE]
  L <- sweep(U, 2L, sqrt(ec$values[keep]), "*")
  W <- sweep(U, 2L, sqrt(ec$values[keep]), "/")
  S <- lapply(seq_along(families), function(f) {
    a <- do.call(rbind, lapply(P, function(z) sweep(z[, (f - 1L) * d + seq_len(d), drop = FALSE], 2L, scale, "/")))
    crossprod(a) / nrow(a)
  })
  kernels <- lapply(S, function(s) {
    e <- eigen(crossprod(W, s %*% W) / 2, symmetric = TRUE)
    list(E = e$vectors, retention = pmin(1, pmax(-0.99, 1 - e$values)),
         clipped = sum(e$values > 1.99))
  })
  propagate <- function(rates = rep(1, length(families))) {
    A <- diag(ncol(L))
    for (f in seq_along(families)) {
      k <- kernels[[f]]
      x <- rates[f]; whole <- floor(x); frac <- x - whole
      retention <- k$retention^whole * (1 - frac + frac * k$retention)
      A <- (k$E %*% (t(k$E) * retention)) %*% A
    }
    J <- tryCatch(solve(diag(nrow(A)) - A), error = function(e) NULL)
    tv <- if (is.null(J)) rep(Inf, d) else pmax(1, rowSums((L %*% (2 * J - diag(nrow(J)))) * L) / pmax(diag(C), 1e-12))
    list(tau = tv, acf = function(lag) {
      Ak <- diag(nrow(A))
      for (k in seq_len(lag)) Ak <- Ak %*% A
      rowSums((L %*% Ak) * L) / pmax(diag(C), 1e-12)
    })
  }
  D <- Reduce(`+`, S)
  ed <- eigen(D, symmetric = TRUE)
  Di <- ed$vectors %*% (t(ed$vectors) / pmax(ed$values, 1e-8 * max(ed$values)))
  additive_tau <- pmax(1, 4 * rowSums((C %*% Di) * C) / pmax(diag(C), 1e-12) - 1)
  list(propagate = propagate, additive_tau = additive_tau,
       clipped = sum(vapply(kernels, `[[`, numeric(1), "clipped")))
}
records <- list(); counterfactuals <- list()
all_names <- names(blocks)
for (b in seq_along(blocks)) {
  block <- blocks[[b]]; name <- all_names[b]; ng <- 2L * block$p
  dataset <- sub("_[^_]+$", "", name)
  future_names <- all_names[seq_along(all_names) > b & startsWith(all_names, paste0(dataset, "_"))]
  fit <- compose_fit(block$X, block$P, block$p)
  pred <- fit$propagate()
  for (f in c(1L, 3L, 4L, 5L)) {
    for (rate in c(0, 0.5, 1, 2)) {
      rates <- rep(1, 5L); rates[f] <- rate
      pr <- fit$propagate(rates)
      counterfactuals[[length(counterfactuals) + 1L]] <- data.frame(block = name, family = families[f], multiplier = rate,
        max_group_tau = max(pr$tau[seq_len(ng)]), median_group_tau = median(pr$tau[seq_len(ng)]))
    }
  }
  if (!length(future_names)) next
  test <- lapply(seq_along(block$X), function(j) head(do.call(rbind, lapply(blocks[future_names], function(bb) bb$X[[j]])), 300L))
  first <- lapply(test, head, 100L)
  variance <- diag(var(do.call(rbind, first)))
  real10 <- 1 - Reduce(`+`, lapply(first, function(x) {
    n <- nrow(x); colMeans((x[11:n, ] - x[seq_len(n - 10L), ])^2)
  })) / length(first) / (2 * variance)
  truth <- vapply(seq_len(ng), function(k) {
    a <- do.call(cbind, lapply(test, function(x) x[, k]))
    length(a) / posterior::ess_basic(a)
  }, numeric(1))
  slow <- order(truth, decreasing = TRUE)[seq_len(min(8L, ng))]
  for (method in c("composition", "additive")) records[[length(records) + 1L]] <- data.frame(block = name, dataset = dataset,
    coord = seq_len(ng), method = method, future_draws = nrow(test[[1]]), slow = seq_len(ng) %in% slow,
    pred_tau = if (method == "composition") pred$tau[seq_len(ng)] else fit$additive_tau[seq_len(ng)],
    true_tau = truth, pred10 = if (method == "composition") pred$acf(10L)[seq_len(ng)] else NA_real_,
    true10 = unname(real10[seq_len(ng)]), clipped = fit$clipped, row.names = NULL)
}
res <- do.call(rbind, records)
write.csv(res, file.path(outdir, "composition_prediction.csv"), row.names = FALSE)
write.csv(do.call(rbind, counterfactuals), file.path(outdir, "composition_counterfactuals.csv"), row.names = FALSE)
for (ds in unique(res$dataset)) for (method in unique(res$method)) {
  x <- res[res$dataset == ds & res$method == method, ]
  cat(ds, method, "ratio all", round(median(x$pred_tau / x$true_tau), 2), "slow8", round(median(x$pred_tau[x$slow] / x$true_tau[x$slow]), 2),
      "acf10 MAE", round(median(abs(x$pred10 - x$true10)), 3), "clipped", unique(x$clipped), "\n")
}
