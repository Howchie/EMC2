# Side-by-side identity check of the C++ ensemble kernels against the R
# reference on a real run: every pool draw and every move is computed both
# ways from the same inputs and RNG state; the chain continues on the C++
# result. Reports selection mismatches, RNG-state mismatches and the largest
# relative numeric differences.
#   EMC_LIB=<new lib> EMC_IN=<fit.rds> ITER=200 Rscript sbs_ens.R [out.rds]
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
emc <- readRDS(Sys.getenv("EMC_IN"))
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
emc <- emc[as.integer(Sys.getenv("CHAIN", "1"))]; class(emc) <- "emc"
iter <- as.integer(Sys.getenv("ITER", "200"))
s0 <- emc[[1]]
cat(sprintf("type %s p %d subjects %d\n", s0$type, s0$n_pars, s0$n_subjects))
rel <- function(a, b) { a <- as.numeric(a); b <- as.numeric(b); if (!length(a)) return(0); max(abs(a - b) / pmax(1, abs(b))) }
acc <- new.env()
acc$draw <- c(calls = 0, rows = 0, lq = 0, L = 0, rng_mismatch = 0)
acc$move <- c(calls = 0, sel_mismatch = 0, ll = 0, tmu = 0, tvar = 0, a_half = 0, wess = 0, rng_mismatch = 0, cpp_path = 0)
seed_get <- function() get(".Random.seed", envir = globalenv())
seed_set <- function(x) assign(".Random.seed", x, envir = globalenv())

draw_pools_ref <- function(sampler, chains_mu, chains_var) {
  cfg <- ns$.emc_ensemble_options()
  n <- sampler$n_subjects; pn <- sampler$par_names
  q <- lapply(seq_len(n), function(s) ns$.emc_ensemble_q(as.numeric(chains_mu[[s]]), as.matrix(chains_var[[s]]), cfg))
  lq <- vector("list", n)
  rows <- lapply(seq_len(n), function(s) {
    fresh <- ns$.emc_ensemble_draw_q(cfg$pool - 1L, q[[s]], cfg)
    lq[[s]] <<- ns$.emc_ensemble_log_q_m2(attr(fresh, "m2"), q[[s]], cfg)
    attr(fresh, "m2") <- NULL; colnames(fresh) <- pn; fresh
  })
  list(q = q, rows = rows, lq = lq)
}
g_draw <- ns$.emc_ensemble_draw_pools
w_draw <- function(sampler, stage, chains_mu, chains_var, iteration = 0L) {
  seed0 <- seed_get()
  out <- g_draw(sampler, stage, chains_mu, chains_var, iteration)
  if (is.null(out)) return(out)
  seed_new <- seed_get()
  seed_set(seed0)
  ref <- draw_pools_ref(sampler, chains_mu, chains_var)
  seed_ref <- seed_get()
  d <- acc$draw; d[["calls"]] <- d[["calls"]] + 1
  n <- sampler$n_subjects
  d[["rows"]] <- max(d[["rows"]], max(vapply(seq_len(n), function(s) rel(out$rows[[s]], ref$rows[[s]]), 0)))
  d[["lq"]] <- max(d[["lq"]], max(vapply(seq_len(n), function(s) rel(out$lq[[s]], ref$lq[[s]]), 0)))
  d[["L"]] <- max(d[["L"]], max(vapply(seq_len(n), function(s) rel(out$L[, , s], ref$q[[s]]$L), 0)))
  d[["rng_mismatch"]] <- d[["rng_mismatch"]] + !identical(seed_ref, seed_new)
  acc$draw <- d
  seed_set(seed_new)
  out
}
environment(w_draw) <- environment(g_draw)
assignInNamespace(".emc_ensemble_draw_pools", w_draw, ns)

# Pure-R group updates for the reference path.
gibbs_diag_gamma_ref <- function(sampler, alpha) {
  last <- ns$last_sample_diag_gamma(sampler$samples)
  prior <- sampler$prior; tvinv <- diag(last$tvinv); n_pars <- sum(!sampler$nuisance); n <- sampler$n_subjects
  P0 <- prior$theta_mu_invar
  if (all(abs(P0[row(P0) != col(P0)]) <= 1e-12 * max(abs(diag(P0))))) {
    var_mu <- 1 / (n * tvinv + diag(P0)); mean_mu <- var_mu * (rowSums(alpha) * tvinv + prior$theta_mu_mean * diag(P0))
    tmu <- rnorm(n_pars, mean_mu, sd = sqrt(var_mu))
  } else {
    prec <- P0 + diag(n * tvinv, n_pars); b <- rowSums(alpha) * tvinv + drop(P0 %*% prior$theta_mu_mean)
    R_post <- chol(prec); tmu <- backsolve(R_post, forwardsolve(t(R_post), b) + rnorm(n_pars))
  }
  names(tmu) <- sampler$par_names[!sampler$nuisance]
  shape <- prior$shape + n / 2; rate <- prior$rate + rowSums((alpha - tmu)^2) / 2
  tvinv <- rgamma(n_pars, shape = shape, rate = rate)
  list(tmu = tmu, tvar = diag(1 / tvinv, n_pars), tvinv = diag(tvinv, n_pars), alpha = alpha)
}
g_gibbs_std <- ns$gibbs_step_standard; g_gibbs_dg <- ns$gibbs_step_diag_gamma
gibbs_std_r <- ns$.gibbs_step_standard_r
set_gibbs <- function(cpp) {
  fs <- if (cpp) g_gibbs_std else gibbs_std_r
  environment(fs) <- environment(g_gibbs_std); assignInNamespace("gibbs_step_standard", fs, ns)
  fd <- if (cpp) g_gibbs_dg else gibbs_diag_gamma_ref
  environment(fd) <- environment(g_gibbs_dg); assignInNamespace("gibbs_step_diag_gamma", fd, ns)
}
acc$gibbs <- c(calls = 0, tmu = 0, tvar = 0, a_half = 0, rng_mismatch = 0)
g_gibbs <- ns$gibbs_step
w_gibbs <- function(sampler, alpha, type, ...) {
  if (!(type %in% c("standard", "diagonal-gamma")) || sys.nframe() > 60) return(g_gibbs(sampler, alpha, type, ...))
  seed0 <- seed_get()
  out <- g_gibbs(sampler, alpha, type, ...)
  seed_new <- seed_get()
  seed_set(seed0); set_gibbs(FALSE)
  ref <- g_gibbs(sampler, alpha, type, ...)
  set_gibbs(TRUE)
  seed_ref <- seed_get()
  m <- acc$gibbs; m[["calls"]] <- m[["calls"]] + 1
  m[["tmu"]] <- max(m[["tmu"]], rel(out$tmu, ref$tmu)); m[["tvar"]] <- max(m[["tvar"]], rel(out$tvar, ref$tvar))
  if (!is.null(ref$a_half)) m[["a_half"]] <- max(m[["a_half"]], rel(out$a_half, ref$a_half))
  m[["rng_mismatch"]] <- m[["rng_mismatch"]] + !identical(seed_ref, seed_new)
  acc$gibbs <- m
  seed_set(seed_new)
  out
}
environment(w_gibbs) <- environment(g_gibbs)
assignInNamespace("gibbs_step", w_gibbs, ns)

g_move <- ns$.emc_ensemble_move
set_cpp <- function(on) {
  f <- if (on) g_in_cpp else function(sampler) FALSE
  environment(f) <- environment(g_in_cpp)
  assignInNamespace(".emc_ensemble_rounds_in_cpp", f, ns)
}
g_in_cpp <- ns$.emc_ensemble_rounds_in_cpp
w_move <- function(sampler, pars, proposals, stage, wpool, wpool_ctx, wpool_part, pools) {
  if (is.null(pools)) return(g_move(sampler, pars, proposals, stage, wpool, wpool_ctx, wpool_part, pools))
  out <- g_move(sampler, pars, proposals, stage, wpool, wpool_ctx, wpool_part, pools)
  set_cpp(FALSE); set_gibbs(FALSE)
  ref <- g_move(sampler, pars, proposals, stage, wpool, wpool_ctx, wpool_part, pools)
  set_cpp(TRUE); set_gibbs(TRUE)
  p <- sampler$n_pars
  m <- acc$move; m[["calls"]] <- m[["calls"]] + 1
  m[["cpp_path"]] <- m[["cpp_path"]] + g_in_cpp(sampler)
  m[["sel_mismatch"]] <- m[["sel_mismatch"]] + !identical(out$proposals[seq_len(p), ], ref$proposals[seq_len(p), ])
  m[["ll"]] <- max(m[["ll"]], rel(out$proposals[p + 1L, ], ref$proposals[p + 1L, ]))
  m[["tmu"]] <- max(m[["tmu"]], rel(out$pars$tmu, ref$pars$tmu))
  m[["tvar"]] <- max(m[["tvar"]], rel(out$pars$tvar, ref$pars$tvar))
  if (!is.null(ref$pars$a_half)) m[["a_half"]] <- max(m[["a_half"]], rel(out$pars$a_half, ref$pars$a_half))
  m[["wess"]] <- max(m[["wess"]], rel(out$sampler$ensemble_stats$subj_wess, ref$sampler$ensemble_stats$subj_wess))
  m[["rng_mismatch"]] <- m[["rng_mismatch"]] + !identical(out$sampler$rng$gibbs, ref$sampler$rng$gibbs)
  acc$move <- m
  out
}
environment(w_move) <- environment(g_move)
assignInNamespace(".emc_ensemble_move", w_move, ns)

set.seed(1)
emc <- suppressMessages(run_emc(emc, stage = "sample", cores_per_chain = as.integer(Sys.getenv("CORES", "9")),
                                cores_for_chains = 1, verbose = FALSE, step_size = iter, max_tries = 1,
                                stop_criteria = list(iter = iter, max_gd = Inf, min_unique = 0, min_es = 0)))
cat("draw_pools:\n"); print(acc$draw)
cat("move:\n"); print(acc$move)
cat("gibbs:\n"); print(acc$gibbs)
ok <- acc$draw[["rng_mismatch"]] == 0 && acc$move[["rng_mismatch"]] == 0 && acc$move[["sel_mismatch"]] == 0 && acc$gibbs[["rng_mismatch"]] == 0
cat(if (ok) "IDENTITY OK\n" else "IDENTITY FAILED\n")
out <- commandArgs(TRUE)
if (length(out)) saveRDS(list(draw = acc$draw, move = acc$move, gibbs = acc$gibbs, type = s0$type, p = s0$n_pars, n = s0$n_subjects, iter = iter), out[1])
