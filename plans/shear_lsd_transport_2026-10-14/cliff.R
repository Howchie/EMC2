suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
grid <- c(-2, -1, 0, 0.5, 1, 1.5, 2)
cases <- list(c(1, 39), c(1, 21), c(1, 58), c(2, 60), c(3, 53))
for (cs in cases) {
  j <- cs[1]; ss <- cs[2]; ch <- emc[[j]]; s <- ch$samples; i <- s$idx; p <- ch$n_pars
  k <- p + ch$shear$curve$l
  idx <- ns$.emc_chain_proposal_window(emc, seq_len(i), "sample")
  ctx <- ns$.emc_shear_context(ch); fits_of <- function(st) ctx$packed
  mu <- as.numeric(s$theta_mu[, i]); V <- s$theta_var[, , i]; alpha <- matrix(s$alpha[, , i], p)
  spec <- ns$.emc_shear_v_spec(ch); walls <- ns$.emc_shear_walls_of(ch); pn <- rownames(s$alpha)
  path <- ns$.emc_shear_coord_path(ch, idx, k, mu, V, grid)
  for (g in seq_along(grid)) path[[g]]$chol <- chol(path[[g]]$V)
  states <- lapply(path, function(st) ns$.emc_shear_state(st$mu, st$V, fits_of(st), ch$prior, spec))
  g0 <- which(grid == 0)
  tr <- lapply(states, function(st) ns$.emc_shear_transport(alpha, states[[g0]], st, walls))
  X <- t(vapply(tr, function(x) x$alpha[, ss], numeric(p))); colnames(X) <- pn
  ll <- ns$.emc_burn_ll_batch(X, ch$data[[ss]], ch$model, pn)
  ldn <- vapply(seq_along(grid), function(g) { z <- backsolve(path[[g]]$chol, X[g, ] - path[[g]]$mu, transpose = TRUE); -0.5 * sum(z^2) - sum(log(diag(path[[g]]$chol))) }, 1)
  jac <- vapply(tr, function(x) x$jac[ss], 1)
  # natural-scale parameters per design cell
  mp <- tryCatch(ns$map_p(ns$add_constants(X, ch$model()$constants), ch$data[[ss]], ch$model()), error = function(e) NULL)
  wall_margin <- if (!is.null(walls) && walls$has[ss]) drop(X %*% walls$U[, ss]) - walls$b[ss] else NA
  wn <- if (!is.null(walls) && walls$has[ss]) paste(pn[abs(walls$U[, ss]) > 1e-8], collapse = "+") else "-"
  cat(sprintf("\nchain %d subj %d (wall normal on: %s)\n", j, ss, wn))
  print(round(cbind(grid = grid, ll = ll - ll[g0], hier = ldn - ldn[g0], jac = jac, wall_margin = wall_margin,
    X[, c("B", "B_Espeed", "B_lRright", "A", "t0")], `exp(t0)` = exp(X[, "t0"]), `exp(B)` = exp(X[, "B"]), `exp(B+BEs)` = exp(X[, "B"] + X[, "B_Espeed"]), `exp(A)` = exp(X[, "A"])), 3))
}
