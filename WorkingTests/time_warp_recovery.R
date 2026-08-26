# Bounded recovery study for the ballistic operational-time warp.
#
# Run from an installed EMC2 package with:
#   Rscript WorkingTests/time_warp_recovery.R
# Set EMC2_SKIP_WARP_FIT=1 to run simulation, the eta profile, and the
# finite-difference identifiability diagnostic without the short MCMC fits.
# The fits deliberately use rt_resolution = NULL: make_emc's default 1/60 s
# floor can bias a parameter whose signal is concentrated at u = rt - t0.

if (!requireNamespace("EMC2", quietly = TRUE))
  stop("Install the EMC2 package before running this recovery study.")
library(EMC2)

set.seed(260826)
run_mcmc <- !identical(Sys.getenv("EMC2_SKIP_WARP_FIT"), "1")
n_trials_per_subject <- 180L
true_eta <- -0.6

recovery_models <- list(
  BAwR = list(
    model = BAwR,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1,
                   kappa ~ 1, p ~ 1, eta ~ 1),
    values = c(mu = .25, sigma = log(.35), B = log(1), A = log(.2),
               t0 = log(.1), kappa = log(.25), p = log(1), eta = true_eta),
    spread = "sigma"
  ),
  LBA = list(
    model = LBA,
    formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, eta ~ 1),
    values = c(v = 2, sv = log(.35), B = log(1), A = log(.2),
               t0 = log(.1), eta = true_eta),
    spread = "sv"
  )
)

# One row per subject is enough for design construction; make_data expands
# those rows into repeated trials and the latent race accumulators.
template <- data.frame(
  subjects = factor(c("s1", "s2")),
  S = factor(c("stim", "stim"), levels = "stim"),
  R = factor(c(NA_character_, NA_character_), levels = c("left", "right"))
)

make_recovery_design <- function(spec) {
  design(data = template, Rlevels = c("left", "right"), model = spec$model,
         matchfun = function(d) as.character(d$S) == as.character(d$lR),
         formula = spec$formula, report_p_vector = FALSE)
}

# This is intentionally an independent, small race likelihood used only for
# the eta profile. It calls each model's public dfun/pfun and therefore does not
# reuse the optimizer or the compiled likelihood under study.
natural_pars <- function(kind, eta, n = 2L) {
  vals <- if (kind == "BAwR") {
    c(mu = .25, sigma = .35, b = 1.2, A = .2, t0 = .1,
      kappa = .25, p = 1)
  } else {
    c(v = 2, sv = .35, b = 1.2, A = .2, t0 = .1)
  }
  vals <- c(vals, eta = eta)
  matrix(rep(unname(vals), length.out = n * length(vals)), nrow = n,
         byrow = TRUE, dimnames = list(NULL, names(vals)))
}

profile_log_likelihood <- function(dat, kind, eta) {
  model <- if (kind == "BAwR") BAwR() else LBA()
  pars <- natural_pars(kind, eta, 2L)
  total <- 0
  for (i in seq_len(nrow(dat))) {
    if (is.finite(dat$rt[i]) && !is.na(dat$R[i])) {
      winner <- if (as.character(dat$R[i]) == "left") 1L else 2L
      loser <- 3L - winner
      density <- model$dfun(dat$rt[i], pars[winner, , drop = FALSE])
      survivor <- 1 - model$pfun(dat$rt[i], pars[loser, , drop = FALSE])
      total <- total + log(max(density * survivor, .Machine$double.xmin))
    } else {
      survivor <- vapply(seq_len(2L), function(j)
        1 - model$pfun(Inf, pars[j, , drop = FALSE]), numeric(1))
      total <- total + sum(log(pmax(survivor, .Machine$double.xmin)))
    }
  }
  total
}

finite_difference_block <- function(ll) {
  h <- .02
  l0 <- ll(0)
  lp <- ll(h)
  lm <- ll(-h)
  c(score_at_eta0 = (lp - lm) / (2 * h),
    hessian_at_eta0 = (lp - 2 * l0 + lm) / h^2,
    eta0_locally_identifiable = is.finite(l0) && is.finite(lp) &&
      is.finite(lm) && (lp - 2 * l0 + lm) < 0)
}
posterior_correlations <- function(fit, spread) {
  empty <- data.frame(subject = character(), cor_eta_t0 = numeric(),
                      cor_eta_B = numeric(), cor_eta_spread = numeric())
  if (is.null(fit))
    return(empty)
  ch <- fit[[1L]]
  keep <- which(ch$samples$stage == "sample")
  if (!length(keep))
    return(empty)
  alpha <- ch$samples$alpha[, , keep, drop = FALSE]
  pn <- dimnames(alpha)[[1L]]
  needed <- c("eta", "t0", "B", spread)
  if (!all(needed %in% pn))
    return(empty)
  out <- lapply(seq_len(dim(alpha)[2L]), function(s) {
    vals <- t(matrix(alpha[needed, s, , drop = FALSE],
                     nrow = length(needed)))
    colnames(vals) <- needed
    cc <- cor(vals, use = "pairwise.complete.obs")["eta", ]
    data.frame(subject = as.character(ch$subjects[s]),
               cor_eta_t0 = unname(cc[["t0"]]),
               cor_eta_B = unname(cc[["B"]]),
               cor_eta_spread = unname(cc[[spread]]))
  })
  do.call(rbind, out)
}

summary_rows <- list()
for (kind in names(recovery_models)) {
  spec <- recovery_models[[kind]]
  des <- make_recovery_design(spec)
  p_true <- sampled_pars(des, doMap = FALSE)
  p_true[names(spec$values)] <- spec$values
  # p_true includes eta ~ 1 explicitly; this is also the starting value used
  # by make_emc before the profile optimizer below is run.
  dat <- make_data(p_true, des, data = template,
                   expand = n_trials_per_subject)
  dat$subjects <- droplevels(dat$subjects)

  ll <- function(e) profile_log_likelihood(dat, kind, e)
  optimizer_start_eta <- 0
  opt <- optim(optimizer_start_eta, function(e) -ll(e), method = "Brent",
               lower = -2, upper = 2)
  fd <- finite_difference_block(ll)

  posterior <- posterior_correlations(
    if (run_mcmc) {
      emc <- make_emc(dat, des, type = "standard", n_chains = 1,
                      compress = FALSE, rt_resolution = NULL)
      fit(emc, iter = 120,
          stop_criteria = list(sample = list(iter = 120,
                                              max_sample_iter = 240)),
          cores_for_chains = 1, cores_per_chain = 1, max_tries = 5,
          particle_factor = 10, verbose = FALSE)
    } else NULL,
    spec$spread
  )

  cat("\n", kind, " (subjects = ", nlevels(dat$subjects), ")\n", sep = "")
  cat("true eta: ", true_eta, "; optimizer start eta: ", optimizer_start_eta,
      "; profile estimate: ", signif(opt$par, 5), "\n", sep = "")
  cat("finite-difference score/Hessian at eta=0: ",
      paste(signif(fd[c("score_at_eta0", "hessian_at_eta0")], 5), collapse = ", "),
      "; eta=0 locally identifiable: ", fd[["eta0_locally_identifiable"]], "\n", sep = "")
  if (nrow(posterior)) print(posterior, row.names = FALSE)
  else cat("posterior correlations: fit skipped (set EMC2_SKIP_WARP_FIT=0 to run)\n")

  summary_rows[[kind]] <- data.frame(
    model = kind, true_eta = true_eta, optimizer_start_eta = optimizer_start_eta,
    profile_eta = unname(opt$par), score_eta0 = unname(fd[["score_at_eta0"]]),
    hessian_eta0 = unname(fd[["hessian_at_eta0"]]),
    eta0_locally_identifiable = unname(fd[["eta0_locally_identifiable"]]),
    posterior_subjects = nrow(posterior), row.names = NULL
  )
}

cat("\nRECOVERY SUMMARY\n")
print(do.call(rbind, summary_rows), row.names = FALSE)
