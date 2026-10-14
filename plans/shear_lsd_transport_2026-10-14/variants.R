# Cheaper Laplace-following frames versus the per-state refit ceiling.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
draws <- as.integer(strsplit(commandArgs(TRUE)[2], ",")[[1]]); j <- as.integer(commandArgs(TRUE)[3])
NC <- 5L
grid <- ns$.EMC_SHEAR_PILOT_GRID
ch0 <- emc[[j]]; p <- ch0$n_pars; cv <- ch0$shear$curve; l <- cv$l; k <- p + l; n <- ch0$n_subjects
pn <- rownames(ch0$samples$alpha)
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch0$samples$idx), "sample")
ctx0 <- ns$.emc_shear_context(ch0)
# conditional mode under N(mu, V) implied by a (lin, prec) likelihood frame
cmode <- function(lin, P, Vi, mu) solve(P + Vi, lin + Vi %*% mu)
# deterministic refit: start at the fixed frame's predicted mode
refit_all <- function(ch, mu, V) {
  Vi <- solve(V)
  fits <- parallel::mclapply(seq_len(n), function(s) {
    x0 <- drop(cmode(ctx0$packed$lin[, s], ctx0$packed$cube[, , s], Vi, mu))
    ns$.emc_shear_node_fit(ch, s, mu, V, x0) }, mc.cores = NC)
  ns$.emc_shear_pack(lapply(seq_len(n), function(s) if (is.null(fits[[s]])) ctx0$fits[[s]] else fits[[s]]), p)
}
# one Newton correction of the fixed frame's predicted mode, Hessian kept
newton1 <- function(ch, mu, V) {
  Vi <- solve(V); out <- ctx0$packed
  res <- parallel::mclapply(seq_len(n), function(s) {
    P <- ctx0$packed$cube[, , s]; m <- drop(cmode(ctx0$packed$lin[, s], P, Vi, mu))
    h <- 1e-4; X <- rbind(t(m + diag(h, p)), t(m - diag(h, p)))
    ll <- ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, pn)
    g <- (ll[1:p] - ll[p + 1:p]) / (2 * h)
    if (!all(is.finite(g))) return(ctx0$packed$lin[, s])
    m2 <- m + drop(solve(P + Vi, g - Vi %*% (m - mu)))
    drop((P + Vi) %*% m2 - Vi %*% mu) }, mc.cores = NC)
  for (s in seq_len(n)) out$lin[, s] <- res[[s]]
  out
}
# damped Newton with the block's Hessian: up to NIT iterations with step
# halving on the conditional objective; deterministic in (mu, V)
NIT <- as.integer(Sys.getenv("NIT", "3")); nev <- new.env(); nev$v <- integer(0)
newtonk <- function(ch, mu, V) {
  Vi <- solve(V); out <- ctx0$packed
  res <- parallel::mclapply(seq_len(n), function(s) {
    P <- ctx0$packed$cube[, , s]; H <- P + Vi; m <- drop(cmode(ctx0$packed$lin[, s], P, Vi, mu)); ev <- 0L
    llf <- function(X) { ev <<- ev + nrow(X); ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, pn) }
    obj <- function(x, ll) ll - 0.5 * drop(crossprod(x - mu, Vi %*% (x - mu)))
    f <- obj(m, llf(matrix(m, 1))); h <- 1e-4
    for (it in seq_len(NIT)) {
      ll <- llf(rbind(t(m + diag(h, p)), t(m - diag(h, p))))
      g <- (ll[1:p] - ll[p + 1:p]) / (2 * h) - drop(Vi %*% (m - mu))
      if (!all(is.finite(g))) break
      step <- drop(solve(H, g)); a <- 1; acc <- FALSE
      for (bt in 1:6) { x <- m + a * step; fx <- obj(x, llf(matrix(x, 1))); if (is.finite(fx) && fx > f) { m <- x; f <- fx; acc <- TRUE; break }; a <- a / 2 }
      if (!acc) break
    }
    list(lin = drop(H %*% m - Vi %*% mu), ev = ev) }, mc.cores = NC)
  for (s in seq_len(n)) out$lin[, s] <- res[[s]]$lin
  nev$v <- c(nev$v, sum(vapply(res, `[[`, 1L, "ev")))
  out
}
# screened refit: refit only subjects whose Gaussian frame mispredicts the
# likelihood change from the block reference to the state's predicted mode
ref_mu <- as.numeric(ch0$samples$theta_mu[, ch0$samples$idx - 1L]); ref_Vi <- solve(ch0$samples$theta_var[, , ch0$samples$idx - 1L])
qf <- function(s, x) sum(ctx0$packed$lin[, s] * x) - 0.5 * drop(crossprod(x, ctx0$packed$cube[, , s] %*% x))
m_ref <- lapply(seq_len(n), function(s) drop(cmode(ctx0$packed$lin[, s], ctx0$packed$cube[, , s], ref_Vi, ref_mu)))
ll_ref <- unlist(parallel::mclapply(seq_len(n), function(s) ns$.emc_burn_ll_batch(matrix(m_ref[[s]], 1), ch0$data[[s]], ch0$model, pn), mc.cores = NC))
TAU <- as.numeric(Sys.getenv("TAU", "0.5")); nref <- new.env(); nref$v <- integer(0)
screened <- function(ch, mu, V) {
  Vi <- solve(V); out <- ctx0$packed
  bad <- unlist(parallel::mclapply(seq_len(n), function(s) {
    m <- drop(cmode(ctx0$packed$lin[, s], ctx0$packed$cube[, , s], Vi, mu))
    ll <- ns$.emc_burn_ll_batch(matrix(m, 1), ch$data[[s]], ch$model, pn)
    !is.finite(ll) || abs((ll - ll_ref[s]) - (qf(s, m) - qf(s, m_ref[[s]]))) > TAU }, mc.cores = NC))
  ss <- which(bad); nref$v <- c(nref$v, length(ss))
  if (length(ss)) {
    fits <- parallel::mclapply(ss, function(s) ns$.emc_shear_node_fit(ch, s, mu, V,
      drop(cmode(ctx0$packed$lin[, s], ctx0$packed$cube[, , s], Vi, mu))), mc.cores = NC)
    for (a in seq_along(ss)) if (!is.null(fits[[a]])) { out$lin[, ss[a]] <- fits[[a]]$lin; out$cube[, , ss[a]] <- fits[[a]]$prec }
  }
  out
}
# all-subject atlas: refits on the rebuild state's lsd path, linear in the lsd offset
reb <- ch0$samples$idx - 1L
chr <- ch0; chr$samples$idx <- reb
apath <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
pr <- ns$.emc_shear_coord_path(chr, idx, k, as.numeric(chr$samples$theta_mu[, reb]), chr$samples$theta_var[, , reb], apath)
lsd_ref <- 0.5 * log(chr$samples$theta_var[l, l, reb]); lsd_sc <- cv$lsd_scale
t_node <- vapply(pr, function(st) (0.5 * log(st$V[l, l]) - lsd_ref) / lsd_sc, 1)
anodes <- NULL
atlas_at <- function(mu, V) {
  t <- min(max((0.5 * log(V[l, l]) - lsd_ref) / lsd_sc, min(t_node)), max(t_node))
  i2 <- max(2L, min(length(t_node), which(t_node >= t)[1])); i1 <- i2 - 1L; w <- (t - t_node[i1]) / (t_node[i2] - t_node[i1])
  list(lin = (1 - w) * anodes[[i1]]$lin + w * anodes[[i2]]$lin, cube = (1 - w) * anodes[[i1]]$cube + w * anodes[[i2]]$cube)
}
cat(sprintf("chain %d, coord %s, rebuild-state draw %d, atlas node t = %s\n", j, pn[l], reb, paste(round(t_node, 2), collapse = " ")))
for (d in draws) {
  ch <- ch0; ch$samples$idx <- d
  ctx <- ns$.emc_shear_context(ch)
  memo <- new.env()
  cache <- function(tag, f) function(st) { key <- paste(tag, signif(c(st$mu, diag(as.matrix(st$V))), 10), collapse = ","); if (is.null(memo[[key]])) memo[[key]] <- f(ch, st$mu, st$V); memo[[key]] }
  t0 <- proc.time()[3]
  sp <- ns$.emc_shear_orbit_spread(ch, idx, k, list(
    function(st) ctx$packed, function(st) ns$.emc_shear_fits_at(ctx, st$mu, st$V),
    cache("refit", refit_all), cache("screen", screened), cache("newtonk", newtonk)), NC)
  cat(sprintf("draw %4d (%3.0f s): total orbit spread  fixed %9.1f | current curve %8.1f | refit %6.1f | screened %7.1f [refits/state %s] | newton%d %7.1f [ll evals/state %s]\n",
    d, proc.time()[3] - t0, sum(sp[, 1]), sum(sp[, 2]), sum(sp[, 3]), sum(sp[, 4]), paste(range(nref$v), collapse="-"), NIT, sum(sp[, 5]), paste(range(nev$v), collapse="-"))); nref$v <- integer(0); nev$v <- integer(0)
}
