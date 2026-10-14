# Check the implemented map rather than just its intended surrogate density.
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
fits <- readRDS(args[1])$fits
ch <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[2])), quiet = TRUE)[[1]]
p <- ch$n_pars; n <- ch$n_subjects
root <- file.path("plans", "shear_lsd_transport_2026-10-14")
load_functions <- function(file, names) {
  for (ex in parse(file)) if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
                            is.symbol(ex[[2]]) && as.character(ex[[2]]) %in% names) eval(ex, globalenv())
}
load_functions(file.path(root, "line_collapse.R"), c("build_pieces", "ldphi", "piece_logmass", "lse", "cdf_side", "inv_cdf", "log_ltilde"))
load_functions(file.path(root, "ridge_first.R"), c("make_state", "map_one"))
d <- ch$samples$idx; mu <- ch$samples$theta_mu[, d]; V <- ch$samples$theta_var[, , d]
scale <- rep(1, p); scale[match("B_Espeed", ch$par_names)] <- exp(0.2)
V1 <- V * outer(scale, scale)
for (s in c(1, 2, 21, 39, 58, 75)) {
  a <- ch$samples$alpha[, s, d]
  st0 <- make_state(mu, V, s); st1 <- make_state(mu, V1, s)
  fw <- map_one(a, st0, st1, s); bw <- map_one(fw$a, st1, st0, s)
  h <- 1e-6
  D <- vapply(seq_len(p), function(j) {
    delta <- replace(numeric(p), j, h)
    (map_one(a + delta, st0, st1, s)$a - map_one(a - delta, st0, st1, s)$a) / (2 * h)
  }, numeric(p))
  jac_error <- abs(as.numeric(determinant(D, logarithm = TRUE)$modulus) - fw$jac)
  inverse_error <- max(abs(bw$a - a))
  identity_error <- max(abs(map_one(a, st0, st0, s)$a - a))
  stopifnot(jac_error < 1e-5, inverse_error < 1e-10, identity_error < 1e-10)
  cat(sprintf("subject %d: Jacobian error %.3g, inverse error %.3g, identity error %.3g\n", s, jac_error, inverse_error, identity_error))
}
