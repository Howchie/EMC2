# Focused observations supporting sampling_refactor_review_2026-10-03.md.
# Run from the repository root. Uses current R source with the existing binary.
# Scoped mocks are restored by testthat; no sampler source is modified.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
cat("REVIEW PROBES: current R worktree, pre-existing src/EMC2.so\n")

# F12: shape-based slicing mistakes a static matrix for sample history.
store <- list(stage = rep("sample", 3), idx = 2L,
              theta_mu = matrix(0, 3, 3), last_theta_var_inv = diag(3))
scratch <- EMC2:::.emc_ensemble_store(store)
cat("Static precision shape in scratch (expected 3 3):",
    dim(scratch$last_theta_var_inv), "\n")
# The actual fill before Gibbs repairs it, so this is not a current crash.
store$theta_var <- array(diag(3), c(3, 3, 3))
store$a_half <- matrix(1, 3, 3)
scratch <- EMC2:::.emc_ensemble_store(store)
scratch <- EMC2:::fill_samples_standard(scratch,
  list(tmu = rep(0, 3), tvar = diag(3), tvinv = diag(3), a_half = rep(1, 3)),
  proposals = NULL, j = 1L, n_pars = 3L)
cat("Static precision shape after standard fill:",
    dim(scratch$last_theta_var_inv), "\n")

# F8: poor-Laplace gate ignores infinite k and NA k with support failure.
k <- c(0.2, Inf, NA_real_, 0.2)
outside <- c(0, 0, 0.9, 0.9)
cat("Current Laplace poor gate:",
    is.finite(k) & (k > 0.7 | outside > 0.5), "\n")

# F7: serial SEM conditions on the intercept; workers receive subject means.
em <- c(a = 1, b = 1)
attr(em, "B") <- diag(0.5, 2)
attr(em, "given") <- c(0, 0)
pars <- list(tmu = c(a = 1, b = 2), subj_mu = matrix(c(4, 6), 2, 1),
             subj_var = diag(2), tvar = diag(2),
             alpha = matrix(c(a = 0, b = 0), 2, 1,
                            dimnames = list(c("a", "b"), NULL)))
settings <- list(list(epsilon = rep(1, 3), mix = c(0, 0, 0, 1),
                      n_particles = 1000L, iter = 0L))
tune <- list(components = c(1L, 1L), shared_ll_idx = c(1L, 1L))
captured <- NULL
mockll <- function(proposals, ...) {
  captured <<- proposals
  rep(0, nrow(proposals))
}
testthat::with_mocked_bindings({
  set.seed(1093)
  invisible(EMC2:::new_particle(1L, NULL, settings, eff_mu = em,
    eff_var = diag(2), chains_mu = c(0, 0), chains_var = diag(2), prev_ll = 0,
    parameters = pars, stage = "sample", type = "SEM", tune = tune))
  serial <- captured
  set.seed(1093)
  invisible(EMC2:::new_particle(1L, NULL, settings, eff_mu = em,
    eff_var = diag(2), chains_mu = c(0, 0), chains_var = diag(2), prev_ll = 0,
    parameters = NULL, stage = "sample", type = "SEM", tune = tune,
    current_alpha = c(a = 0, b = 0), population_mu = c(a = 4, b = 6),
    population_var = diag(2)))
  cat("SEM worker minus serial fresh proposal means:",
      colMeans(captured[-1, ] - serial[-1, ]), "\n")
}, calc_ll_pooled = mockll, .package = "EMC2")

# F9: raw eigenvalue regularization is not equivariant to coordinate units.
S <- matrix(1, 2, 2)
D <- diag(c(1, 1e6))
reg1 <- EMC2:::.emc_regularize_covariance(S)
reg2 <- solve(D) %*% EMC2:::.emc_regularize_covariance(D %*% S %*% D) %*% solve(D)
cat("Regularized correlation, original/rescaled units:",
    cov2cor(reg1)[1, 2], cov2cor(reg2)[1, 2], "\n")

# F5: population moment formulas, then actual implied-mean helper.
cat("SEM covariate variance contribution, current/correct:",
    0.3^2 + 0.8^2, (0.3 + 0.8)^2, "\n")
x <- list(samples = list(theta_mu = matrix(0, 2, 2),
                        B = array(0, c(1, 1, 2)), G = array(0.8, c(1, 1, 2)),
                        K = array(0.3, c(2, 1, 2)), lambda = array(1, c(2, 1, 2))),
          sem_settings = list(covariates = matrix(c(0, 2), 2, 1)))
cat("SEM implied mean helper (each element should be 1.1):\n")
print(EMC2:::get_mu_implied(x, 1:2))

# F4: refitting a capped random subset retains noise even for a fixed history.
set.seed(1369)
p <- 2L
n <- 10000L
mu <- matrix(rnorm(p * n), p)
a <- 0.7 * mu + matrix(rnorm(p * n), p)
sampler <- list(type = "standard", nuisance = c(FALSE, FALSE),
  samples = list(stage = rep("sample", n), idx = n, theta_mu = mu,
                 theta_var = array(rep(diag(p), n), c(p, p, n)),
                 alpha = array(a, c(p, 1, n))))
co <- replicate(100L, {
  z <- EMC2:::extract_samples(sampler, max_n_sample = 750L, n_chains = 1L)
  y <- EMC2:::.emc_mean_conditional(1L, z, p, idx = rep(TRUE, p))
  c(y$eff_var[1, 1], y$B[1, 1])
})
cat("Repeated efficient refits SD (variance/slope):", apply(co, 1, sd), "\n")

# F1: retain a known score draw and observe the actual loading conditional.
# Row 2 has lambda21 free and lambda22 fixed at 1. With y2 = eta2,
# the residual response is zero; the correct conditional mean for lambda21 is 0.
eta <- cbind(c(-2, -1, 0, 1, 2), c(-2, -1, 0, 1, 2))
alpha <- rbind(rep(0, 5), eta[, 2], rep(0, 5))
s <- list(n_subjects = 5L, nuisance = rep(FALSE, 3), n_factors = 2L,
  samples = list(idx = 1L, theta_mu = matrix(0, 3, 1),
    eta = array(eta, c(5, 2, 1)), lambda_untransf = array(0, c(3, 2, 1)),
    psi_inv = matrix(1, 2, 1), epsilon_inv = matrix(1, 3, 1)),
  prior = list(theta_mu_invar = diag(3), theta_mu_mean = rep(0, 3),
    as = 1, bs = 1, theta_lambda_invar = rep(1, 3), ap = 1, bp = 1))
attr(s, "Lambda_mat") <- matrix(c(1, Inf, 0, 0, 1, 0), 3, 2)
s$samples$lambda_untransf[, , 1] <- matrix(c(1, 0, 0, 0, 1, 0), 3, 2)
mock_mvn <- function(n, mean = rep(0, 3), sigma) {
  # Current anchor matrix gives conditional score mean (0, .5 * y2).
  # Return residuals making the final score draw equal the chosen eta exactly.
  if (n == 5L) return(eta - cbind(0, 0.5 * alpha[2, ]))
  if (length(mean) == 3L) return(matrix(0, 1, 3))
  cat("Factor row-2 free loading conditional mean (correct = 0):",
      as.numeric(mean), "\n")
  matrix(as.numeric(mean), 1)
}
testthat::with_mocked_bindings({
  set.seed(2049)
  invisible(EMC2:::gibbs_step_factor(s, alpha))
}, rmvnorm = mock_mvn, .package = "EMC2")

# F2/F3/F13: display algebraic counterexamples; these are not full-chain tests.
old_delta <- c(1, 2, 3)
new_delta <- c(1, 4, 3)
cat("Third cumulative shrinkage after delta2 changes, stale/correct:",
    cumprod(old_delta)[3], cumprod(new_delta)[3], "\n")
cat("Two-factor cycle determinant with B12 = .2, B21 = .3:",
    1 - 0.2 * 0.3, "(depends on both free coefficients)\n")
cat("Student t4: E[M2]/p = 2*c; Var(M2) is infinite.\n")
