.summarise_ess <- function(fit, label, n_draws, rebuild_seconds = NA_real_,
                           sample_seconds = NA_real_, cores_per_chain = NA_integer_) {
  subject <- "52091"
  params <- c("A", "B", "B_lR_target", "t0_LoadHigh", "t0_LoadLow")
  par_names <- dimnames(fit[[1L]]$samples$theta_mu)[[1L]]
  p <- dim(fit[[1L]]$samples$alpha)[1L]
  stopifnot(length(par_names) == p, all(params %in% par_names))
  key <- array(NA_real_, c(n_draws, length(fit), 2L * length(params)),
    dimnames = list(NULL, NULL,
      c(paste0("mu:", params), paste0(subject, ":", params))))
  group <- array(NA_real_, c(n_draws, length(fit), p),
    dimnames = list(NULL, NULL, par_names))
  person <- array(NA_real_, c(n_draws, length(fit), p),
    dimnames = list(NULL, NULL, par_names))
  for (j in seq_along(fit)) {
    ch <- fit[[j]]
    idx <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n_draws)
    stopifnot(length(idx) == n_draws)
    s <- match(subject, dimnames(ch$samples$alpha)[[2L]])
    stopifnot(!is.na(s))
    group[, j, ] <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    person[, j, ] <- t(matrix(ch$samples$alpha[, s, idx], nrow = p, ncol = n_draws))
    for (k in seq_along(params)) {
      par <- params[k]
      key[, j, k] <- ch$samples$theta_mu[par, idx]
      key[, j, length(params) + k] <- ch$samples$alpha[par, s, idx]
    }
  }
  kd <- posterior::as_draws_array(key)
  ed <- posterior::summarise_draws(kd, posterior::ess_bulk,
                                    posterior::ess_tail, posterior::rhat)
  names(ed)[match("variable", names(ed))] <- "variable"
  names(ed)[grep("ess_bulk$", names(ed))] <- "bulk"
  names(ed)[grep("ess_tail$", names(ed))] <- "tail"
  names(ed)[grep("rhat$", names(ed))] <- "rhat"
  gd <- posterior::as_draws_array(group)
  pd <- posterior::as_draws_array(person)
  list(label = label, draws_per_chain = n_draws,
    cores_for_chains = length(fit), cores_per_chain = cores_per_chain,
    rebuild_seconds = rebuild_seconds, sample_seconds = sample_seconds,
    total_seconds = if (all(is.finite(c(rebuild_seconds, sample_seconds))))
      rebuild_seconds + sample_seconds else NA_real_,
    collapsed_rounds = vapply(fit, function(ch)
      if (is.null(ch$collapsed$plan$rounds)) 0 else ch$collapsed$plan$rounds,
      numeric(1)),
    collapsed_scale = vapply(fit, function(ch)
      if (is.null(ch$collapsed$plan$scale)) 1 else ch$collapsed$plan$scale,
      numeric(1)),
    collapsed_block = lapply(fit, function(ch) ch$collapsed$plan$block),
    collapsed_block_acceptance = vapply(fit, function(ch) {
      b <- ch$collapsed$plan$block
      if (is.null(b) || is.null(b$n) || b$n == 0L) NA_real_ else b$acc / b$n
    }, numeric(1)),
    calibration = lapply(fit, function(ch) ch$collapsed$plan$calibration),
    key_ess = ed,
    group_mu_bulk = posterior::summarise_draws(gd, posterior::ess_bulk)[[2L]],
    group_mu_tail = posterior::summarise_draws(gd, posterior::ess_tail)[[2L]],
    subject_52091_bulk = posterior::summarise_draws(pd, posterior::ess_bulk)[[2L]],
    subject_52091_tail = posterior::summarise_draws(pd, posterior::ess_tail)[[2L]],
    retained_draws = key)
}
