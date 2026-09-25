## Minimal witness for the between-run variance discrepancy in the saved fits.
## Reads the two existing fit files; it does not rerun the expensive sampler.
## REPRO_FITS_DIR may point to a copy of the project Fits directory.
suppressMessages(library(EMC2))

fits_dir <- Sys.getenv("REPRO_FITS_DIR", "/data/work/PM/NirvanaHons_Nback/Fits")
files <- c(
  "Match_Control_Exp1_rdm7s_mass_noA_t0R_vs2_mt_d42.RData",
  "Match_Control_Exp1_rdm7s_mass_noA_t0R_UT0_vs2_mt_d42.RData"
)
fits <- lapply(file.path(fits_dir, files), function(path) {
  if (!file.exists(path)) stop("Missing fit: ", path)
  box <- new.env(parent = emptyenv())
  load(path, envir = box)
  box$emc
})
stopifnot(all(vapply(fits, length, integer(1)) == 3L))

same_prior <- all(vapply(names(fits[[1]][[1]]$prior), function(name)
  identical(fits[[1]][[1]]$prior[[name]], fits[[2]][[1]]$prior[[name]]), logical(1)))
same_trials <- all(vapply(seq_along(fits[[1]][[1]]$data), function(s) {
  a <- fits[[1]][[1]]$data[[s]]
  b <- fits[[2]][[1]]$data[[s]]
  identical(names(a), names(b)) && all(vapply(names(a), function(name)
    identical(a[[name]], b[[name]]), logical(1)))
}, logical(1)))
cat("Same prior values:", same_prior, "\nSame trial columns:", same_trials, "\n")
stopifnot(same_prior, same_trials)

## Check the same parameter vectors in both current likelihood implementations.
max_ll_gap <- 0
for (s in c(1L, 7L, 22L, 42L)) {
  p <- rbind(t(fits[[1]][[1]]$samples$alpha[, s, c(1L, 500L, 1000L)]),
             t(fits[[2]][[1]]$samples$alpha[, s, c(1L, 500L, 1000L)]))
  colnames(p) <- fits[[1]][[1]]$par_names
  ll <- lapply(fits, function(x) EMC2:::calc_ll_manager(
    p, x[[1]]$data[[s]], x[[1]]$model))
  max_ll_gap <- max(max_ll_gap, max(abs(ll[[1]] - ll[[2]])))
}
cat("Largest tested log-likelihood gap:", max_ll_gap, "\n")

sd_chains <- lapply(fits, function(x) lapply(x, function(ch) {
  v <- ch$samples$theta_var
  ans <- vapply(seq_len(dim(v)[3L]), function(t) sqrt(diag(v[, , t])),
                numeric(dim(v)[1L]))
  rownames(ans) <- x[[1L]]$par_names
  ans
}))
medians <- lapply(sd_chains, function(chains) apply(do.call(cbind, chains), 1L, median))
within <- vapply(fits, function(x)
  max(unlist(gd_summary(x, selection = "sigma2")), na.rm = TRUE), numeric(1))
all_chains <- unlist(sd_chains, recursive = FALSE)
pooled <- coda::mcmc.list(lapply(all_chains, function(x) coda::mcmc(t(x))))
between <- coda::gelman.diag(pooled, multivariate = FALSE,
                             autoburnin = FALSE)$psrf[, 1L]

## Ask the actual stopping function whether including sigma2 would have
## prolonged either saved fit at its 1000th sample draw.
would_stop <- vapply(fits, function(x) {
  progress <- EMC2:::check_progress(
    x, stage = "sample", iter = 1000L,
    stop_criteria = list(iter = 1000L, max_gd = 1.1,
                         selection = c("alpha", "mu", "sigma2"),
                         omit_mpsrf = TRUE),
    max_tries = 20L, step_size = 100L, n_cores = 27L,
    verbose = FALSE, n_blocks = 1L, rhat_version = "old")
  isTRUE(progress$done)
}, logical(1))

cat("Within-fit max sigma2 Rhat:", paste(within, collapse = ", "), "\n")
cat("Would stop at 1000 with sigma2 included:",
    paste(would_stop, collapse = ", "), "\n")
cat("Six-chain max SD Rhat:", max(between), "\n")
key <- "v.d_match_LoadHigh"
cat(key, "SD medians:", paste(round(vapply(medians, `[[`, numeric(1), key), 3),
                                collapse = ", "), "\n")
stopifnot(all(would_stop), max(between) > 1.1)
