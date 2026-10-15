# Compare scalar and joint linear predictors on held-out blocks of saved draws.
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
outdir <- file.path(root, "plans/alloc_prediction_2026-10-15/eval_operator")
scratch <- file.path(root, ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/19c4c6df-25fb-4858-8ef0-fdc4d9b6acea/scratchpad")
suppressPackageStartupMessages({library(EMC2, lib.loc = Sys.getenv("EMC_LIB", file.path(root, ".Rlib_alloc_review"))); library(posterior)})
ns <- asNamespace("EMC2")
ev <- new.env(parent = ns)
sys.source(file.path(outdir, "evaluated_allocation.R"), ev)
datasets <- list()
add_fit <- function(name, emc, stages) {
  p <- emc[[1]]$n_pars
  for (stage in stages) {
    idx <- which(emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)] == stage)
    if (length(idx) < 200L) next
    datasets[[paste(name, stage)]] <<- lapply(emc, function(ch) {
      s <- ch$samples
      cbind(t(s$theta_mu[, idx, drop = FALSE]), t(ev$.emc_alloc_lsd_draws(s$theta_var, idx)),
            t(apply(s$alpha[, , idx, drop = FALSE], c(1L, 3L), mean)))
    })
  }
}
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
add_fit("PM", emc, "burn")
pm_dirs <- solve(emc[[1]]$shear$plan$dirs)
for (name in c("full_T30", "full_forst"))
  add_fit(name, ns$restore_duplicates(readRDS(file.path(scratch, paste0(name, "_emc.rds")))), c("burn", "sample"))
for (name in c("pair_new_2", "pair_new_3", "pair_point_2", "pair_point_3", "arm_path_2")) {
  r <- readRDS(file.path(root, "plans/alloc_prediction_2026-10-15/handoff_2026-10-15/arms", paste0(name, ".rds")))
  datasets[[name]] <- if (!is.null(r$alpha_mean)) Map(cbind, r$draws, r$alpha_mean) else r$draws
}
for (name in names(datasets)[grepl("^PM |^pair_|^arm_path", names(datasets))]) {
  datasets[[paste0(name, " dirs")]] <- lapply(datasets[[name]], function(x) {
    p <- nrow(pm_dirs)
    cbind(x[, seq_len(p)] %*% t(pm_dirs), x[, -(seq_len(p)), drop = FALSE])
  })
}
tau <- function(r) (1 + r) / pmax(1 - r, 1e-10)
lag_r <- function(X, lag) {
  v <- diag(var(do.call(rbind, X)))
  d <- Reduce(`+`, lapply(X, function(x) {
    n <- nrow(x)
    colMeans((x[(lag + 1L):n, , drop = FALSE] - x[seq_len(n - lag), , drop = FALSE])^2)
  })) / length(X)
  1 - d / (2 * v)
}
joint <- function(X, shrink = 0, lag = 1L) {
  xx <- do.call(rbind, X)
  sd <- sqrt(diag(var(xx)))
  sd[!is.finite(sd) | sd <= 0] <- 1
  Z <- lapply(X, function(x) sweep(x, 2L, sd, "/"))
  C <- var(do.call(rbind, Z))
  de <- do.call(rbind, lapply(Z, function(x) {
    n <- nrow(x)
    x[(lag + 1L):n, , drop = FALSE] - x[seq_len(n - lag), , drop = FALSE]
  }))
  S <- crossprod(de) / nrow(de)
  S <- (1 - shrink) * S + shrink * diag(diag(S))
  ec <- eigen(C, symmetric = TRUE)
  keep <- ec$values > 1e-8 * max(ec$values)
  U <- ec$vectors[, keep, drop = FALSE]
  W <- sweep(U, 2L, sqrt(ec$values[keep]), "/")
  e <- eigen(crossprod(W, S %*% W) / 2, symmetric = TRUE)
  rates <- pmin(1 - 1e-9, pmax(0, 1 - e$values))^(1 / lag)
  modes <- sweep(U, 2L, sqrt(ec$values[keep]), "*") %*% e$vectors
  weights <- modes^2 / pmax(diag(C), 1e-12)
  list(tau = drop(weights %*% tau(rates)), acf = function(k) drop(weights %*% rates^k))
}
quadratic <- function(X, pcs = 5L) {
  xx <- do.call(rbind, X)
  m <- colMeans(xx); sd <- sqrt(diag(var(xx))); sd[sd <= 1e-12] <- 1
  Z <- lapply(X, function(x) sweep(sweep(x, 2L, m), 2L, sd, "/"))
  e <- eigen(var(do.call(rbind, Z)), symmetric = TRUE)
  E <- e$vectors[, seq_len(min(pcs, ncol(xx))), drop = FALSE]
  lapply(seq_along(X), function(i) {
    a <- Z[[i]] %*% E
    products <- lapply(seq_len(ncol(a)), function(j) a[, j] * a[, j:ncol(a), drop = FALSE])
    cbind(X[[i]], do.call(cbind, products))
  })
}
records <- list()
for (name in names(datasets)) {
  X <- datasets[[name]]
  d <- ncol(X[[1]])
  # Augmented statistics append p subject means to the 2p group coordinates.
  augmented <- grepl("^PM |^arm_path_2|^full_", name)
  ng <- if (augmented) 2L * d / 3L else d
  T <- min(vapply(X, nrow, integer(1)))
  if (T < 200L) next
  for (b in seq_len(T %/% 100L - 1L)) {
    tr <- (b - 1L) * 100L + seq_len(100L)
    te <- b * 100L + seq_len(100L)
    train <- lapply(X, function(x) x[tr, , drop = FALSE])
    test <- lapply(X, function(x) x[te, seq_len(ng), drop = FALSE])
    real10 <- lag_r(test, 10L)
    realtau <- ev$.emc_alloc_geyer_tau(test)
    future_idx <- (b * 100L + 1L):min(T, b * 100L + 300L)
    reference <- lapply(X, function(x) x[future_idx, seq_len(ng), drop = FALSE])
    trun <- vapply(seq_len(ng), function(k) {
      a <- do.call(cbind, lapply(reference, function(x) x[, k]))
      length(a) / posterior::ess_basic(a)
    }, numeric(1))
    slow <- order(trun, decreasing = TRUE)[seq_len(min(8L, ng))]
    scalar <- pmax(0, pmin(1 - 1e-9, lag_r(lapply(train, function(x) x[, seq_len(ng), drop = FALSE]), 1L)))
    preds <- list(scalar = list(tau = tau(scalar), acf = function(k) scalar^k))
    for (shrink in c(0, 0.01, 0.1)) {
      preds[[paste0("joint_group_", shrink)]] <- joint(lapply(train, function(x) x[, seq_len(ng), drop = FALSE]), shrink)
      if (augmented) preds[[paste0("joint_aug_", shrink)]] <- joint(train, shrink)
    }
    for (lag in c(5L, 10L, 25L)) {
      preds[[paste0("joint_group_lag", lag)]] <- joint(lapply(train, function(x) x[, seq_len(ng), drop = FALSE]), lag = lag)
      if (augmented) preds[[paste0("joint_aug_lag", lag)]] <- joint(train, lag = lag)
      rs <- pmax(0, pmin(1 - 1e-9, lag_r(lapply(train, function(x) x[, seq_len(ng), drop = FALSE]), lag)))^(1 / lag)
      # Capture each loop iteration's rates in its own closure.
      preds[[paste0("scalar_lag", lag)]] <- local({rr <- rs; list(tau = tau(rr), acf = function(k) rr^k)})
    }
    for (pcs in c(3L, 5L)) {
      z <- quadratic(train, pcs)
      preds[[paste0("joint_quad", pcs)]] <- joint(z)
    }
    for (method in names(preds)) {
      pr <- preds[[method]]
      records[[length(records) + 1L]] <- data.frame(dataset = name, block = b, method = method, coord = seq_len(ng),
        slow = seq_len(ng) %in% slow, future_draws = length(future_idx), pred_tau = pr$tau[seq_len(ng)], real_tau = realtau, run_tau = trun,
        pred10 = pr$acf(10L)[seq_len(ng)], real10 = real10)
    }
  }
  cat(name, "blocks", T %/% 100L, "dimensions", d, "\n")
}
res <- do.call(rbind, records)
write.csv(res, file.path(outdir, "prediction.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(res, list(res$dataset, res$method), drop = TRUE), function(x)
  data.frame(dataset = x$dataset[1], method = x$method[1],
    acf10_mae = median(abs(x$pred10 - x$real10)), acf10_bias = median(x$pred10 - x$real10),
    tau_ratio = median(x$pred_tau / x$real_tau), slow_tau_ratio = median(x$pred_tau[x$slow] / x$real_tau[x$slow]),
    run_ratio = median(x$pred_tau / x$run_tau), slow_run_ratio = median(x$pred_tau[x$slow] / x$run_tau[x$slow]))))
write.csv(summary, file.path(outdir, "prediction_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
