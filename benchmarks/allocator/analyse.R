# Calibration tables for one calibrate.R run: every rebuild's predicted lag-one
# rho per coordinate against what the following block realised.
#   EMC_LIB=<lib> Rscript benchmarks/allocator/analyse.R <run_dir> [label]
# Writes <run_dir>/calib.rds (one row per chain x coordinate x rebuild) and
# prints the summary by coordinate kind.
#
# Columns:
#   pred      the rebuild's prediction under its plan (lower bound on the rest
#             of the iteration, as solved);
#   pred_pt   the same plan with the rest of the iteration at its point
#             estimate (the confidence bound removed);
#   pred0     the rest of the iteration with nothing bought (subjects under
#             their particle step alone);
#   rest_real the rest of the iteration's lag-one rho in the block (trace);
#   real      the whole iteration's in the block, by the allocator's
#             estimator (ESJD against the variance the rebuild used);
#   real_acf  the standard lag-one autocorrelation of the block's draws;
#   iact      the block's integrated autocorrelation time (single chain).
args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(EMC2, lib.loc = Sys.getenv("EMC_LIB"))
  library(posterior)
})
run_dir <- args[1]
label <- if (length(args) >= 2) args[2] else basename(run_dir)
run <- readRDS(file.path(run_dir, "run.rds"))
lsd <- function(V) 0.5 * log(pmax(diag(as.matrix(V)), .Machine$double.xmin))
tau_of <- function(r) { r <- pmin(r, 1 - 1e-9); (1 + r) / (1 - r) }
# Undo the confidence bound: rho_lo = 1 - m (1 + 2 sqrt(2 / n)).
point <- function(rho_lo, n) ifelse(rho_lo > 0 & rho_lo < 1,
                                    1 - (1 - rho_lo) / (1 + 2 * sqrt(2 / max(1, n))), rho_lo)
rho_esjd <- function(x, var) 1 - mean(diff(x)^2) / (2 * var)
rho_acf <- function(x) { x <- x - mean(x); sum(x[-1] * x[-length(x)]) / sum(x^2) }
iact <- function(x) { e <- ess_basic(matrix(x, ncol = 1)); if (is.finite(e) && e > 0) (length(x)) / e else NA_real_ }

starts <- c(0L, run$offsets)
ends <- c(run$offsets, run$iters)
rows <- list()
for (k in seq_along(starts)) {
  f <- file.path(run_dir, sprintf("alloc_%02d.rds", k))
  if (!file.exists(f)) next
  cap <- readRDS(f)
  nxt <- readRDS(file.path(run_dir, sprintf("pre_%02d.rds", k + 1L)))
  blk <- (starts[k] + 1L):(ends[k] + 1L)  # draw indices incl. the state at the rebuild
  for (j in seq_along(cap$records)) {
    rec <- cap$records[[j]]
    m <- cap$measure[[which(vapply(cap$measure, `[[`, integer(1), "j") == j)[1L]]]
    ms <- m$ms
    p <- nrow(run$draws[[j]]$mu)
    n <- dim(run$draws[[j]]$alpha)[2]
    dirs <- if (is.matrix(m$dirs) && isTRUE(ms$shear)) m$dirs else diag(p)
    Dinv <- solve(dirs)
    D <- run$draws[[j]]
    Z <- rbind(Dinv %*% D$mu[, blk, drop = FALSE], D$lsd[, blk, drop = FALSE])
    tr <- nxt$emc[[j]]$ws1_trace
    # Rest-of-iteration jumps in the block, in this rebuild's coordinates.
    Gd <- if (length(tr$group)) do.call(cbind, tr$group) else NULL
    Gd <- if (!is.null(Gd)) rbind(Dinv %*% Gd[seq_len(p), , drop = FALSE], Gd[p + seq_len(p), , drop = FALSE])
    n_g <- 2L * p
    xg <- rec$rounds
    for (c in seq_len(n_g)) {
      v <- ms$g_var[c]
      rows[[length(rows) + 1L]] <- data.frame(
        fixture = label, k = k, chain = j, kind = if (c <= p) "mu_dir" else "lsd",
        coord = c, subject = NA_integer_, bought = xg[c], collapsed = rec$collapsed_rounds,
        pred = rec$pred_group[c],
        pred_pt = if (rec$pred0_group[c] > 0) rec$pred_group[c] / rec$pred0_group[c] *
          point(rec$pred0_group[c], ms$pairs[["group"]]) else rec$pred_group[c],
        pred0 = rec$pred0_group[c], rest_in = point(ms$g_rho[c], ms$pairs[["group"]]),
        g = ms$g_fac[c], a = ms$g_sat[c],
        rest_real = if (!is.null(Gd)) 1 - mean(Gd[c, ]^2) / (2 * v) else NA_real_,
        real = rho_esjd(Z[c, ], v), real_acf = rho_acf(Z[c, ]), iact = iact(Z[c, -1L]),
        var_ratio = var(Z[c, ]) / v)
    }
    # Subjects: whitened alpha.
    Ad <- if (length(tr$subjects)) lapply(tr$subjects, `[[`, "d") else NULL
    for (s in seq_len(n)) {
      W <- ms$W[[s]]
      X <- W %*% matrix(D$alpha[, s, blk], p)
      Ra <- if (!is.null(Ad)) W %*% vapply(Ad, function(d) d[, s], numeric(p)) else NULL
      for (c in seq_len(p)) {
        v <- ms$s_var[s, c]
        pr <- rec$pred_subject[s, c]
        p0 <- rec$pred0_subject[s, c]
        rows[[length(rows) + 1L]] <- data.frame(
          fixture = label, k = k, chain = j, kind = "subject", coord = c, subject = s,
          bought = rec$rate[s], collapsed = rec$particle_rate[s],
          pred = pr,
          pred_pt = if (ms$s_rho[s, c] > 0) pr / ms$s_rho[s, c] *
            point(ms$s_rho[s, c], ms$pairs[["subjects"]]) else pr,
          pred0 = p0, rest_in = point(ms$s_rho[s, c], ms$pairs[["subjects"]]),
          g = if (all(is.finite(ms$p_fac[s, ]))) ms$p_fac[s, c] else NA_real_,
          a = NA_real_,
          rest_real = if (!is.null(Ra)) 1 - mean(Ra[c, ]^2) / (2 * v) else NA_real_,
          real = rho_esjd(X[c, ], v), real_acf = rho_acf(X[c, ]), iact = iact(X[c, -1L]),
          var_ratio = var(X[c, ]) / v)
      }
    }
  }
  rm(nxt); gc(verbose = FALSE)
}
tab <- do.call(rbind, rows)
saveRDS(tab, file.path(run_dir, "calib.rds"))

summ <- function(t) {
  e <- t$real - t$pred
  ept <- t$real - t$pred_pt
  lt <- log(tau_of(pmax(t$real, 0)) / tau_of(t$pred_pt))
  data.frame(n = nrow(t), pred = mean(t$pred), pred_pt = mean(t$pred_pt),
             real = mean(t$real), real_acf = mean(t$real_acf),
             bias = mean(e), bias_pt = mean(ept), sd_pt = sd(ept),
             logtau_ratio = median(lt), iact_over_tau = median(t$iact / tau_of(t$pred_pt), na.rm = TRUE),
             rest_in = mean(t$rest_in), rest_real = mean(t$rest_real, na.rm = TRUE))
}
cat("\n== ", label, "\n")
for (kd in c("mu_dir", "lsd", "subject")) {
  t <- tab[tab$kind == kd, ]
  if (!nrow(t)) next
  cat(kd, "\n")
  print(round(summ(t), 3))
}
# The coordinates the minimax was about: per chain and rebuild, the slowest
# predicted coordinate against the slowest realised.
worst <- do.call(rbind, lapply(split(tab, list(tab$k, tab$chain), drop = TRUE), function(t) {
  wp <- which.max(t$pred_pt); wr <- which.max(t$real)
  data.frame(k = t$k[1], chain = t$chain[1],
             pred_worst = paste(t$kind[wp], t$coord[wp], t$subject[wp]), pred_rho = t$pred_pt[wp],
             real_at_pred = t$real[wp],
             real_worst = paste(t$kind[wr], t$coord[wr], t$subject[wr]), real_rho = t$real[wr],
             pred_at_real = t$pred_pt[wr])
}))
print(worst, digits = 3)
