# Arithmetic-only review probes. Run from the repository root with Rscript.
# No package install, likelihood evaluations, MCMC transitions, fits, workers,
# checkpoint mutations, or random draws. Assertions pin the observed defects;
# they are review reproductions, not acceptance tests for a future fix.
source("R/move_allocation.R")
source("R/subject_hmc.R")
source("R/shear_move.R")
options(emc2.alloc_hmc_regimes = TRUE)

emit <- function(label, ...) {
  cat("\n", label, "\n", sep = "")
  print(list(...))
}

# 1. Coordinates reachable only through the joint move are excluded.
units <- list(list(kind = "fixed", rho = 0.95),
              list(kind = "fixed", rho = 0.9))
joint <- list(idx = 1:2, g = c(0.3, 0.4), cost = 0.05, cap = 4)
excluded <- .emc_alloc_solve(units, 1, joint = joint)
included <- .emc_alloc_solve(units, 1, joint = joint,
                             movable = c(TRUE, TRUE))
stopifnot(excluded$joint == 0, included$joint > 0)
emit("Collapsed-only reachability", rounds_excluded = excluded$joint,
     rounds_included = included$joint, tau0_excluded = excluded$tau0,
     tau0_all = excluded$tau0_all)

# 2. The measurement's upper clamp invents movement from zero ESJD.
g <- min(1 - 1e-6, max(0, 1 - 0 / (2 * 1)))
zero <- .emc_alloc_solve(list(list(kind = "shear", rho = 1, g = g,
                                  cost = 0.05)), 1)
stopifnot(zero$x == 4)
emit("Measured zero movement", factor = g, rounds = zero$x,
     predicted_tau = zero$tau, baseline_tau = zero$tau0,
     predicted_gain = 1 - zero$objective / zero$objective0)

# 3. A regime plus supplemental HMC models two kernels, but executes one
# Bernoulli transition at the capped sum of their rates.
u <- list(units = list(
  list(kind = "shear", rho = 0.8436096, g = 0.9562901, cost = 0.1575721),
  list(kind = "subject", rho = 0.8590438, p = 0.6462185, c_p = 0.3781685,
       h = 0.1771749, h_dev = 0.1771749, d_rho = 0.5206982,
       c_h = 0.1249796)), C_rest = 1, hmc_fixed = 0.01, joint = NULL,
  r_bar = 0.5206982, g_floor = 0.1)
s <- .emc_alloc_solve_units(u)
r <- s$regime; x <- s$x[2]; h <- u$units[[2]]$h
executed <- min(1, r + x)
modeled_factor <- (1 - r * (1 - h)) * (1 - x * (1 - h))
executed_factor <- 1 - executed * (1 - h)
stopifnot(r > 0, x > 0, abs(modeled_factor - executed_factor) > 0.01)
emit("Regime versus execution", regime = r, supplemental_rate = x,
     executed_rate = executed, modeled_factor = modeled_factor,
     executed_factor = executed_factor,
     modeled_HMC_work = (r + x) * u$units[[2]]$c_h,
     executed_HMC_work = executed * u$units[[2]]$c_h)

# 4. A pure regime selects a subject that the pilot gate does not select.
u <- list(units = list(list(kind = "shear", rho = 0.98, g = 0.99, cost = 10),
  list(kind = "subject", rho = 0.01, p = NULL, c_p = 0,
       h = 0, h_dev = 0, d_rho = 0.8, c_h = 0.1)),
  C_rest = 1, hmc_fixed = 0, joint = NULL, r_bar = 0.8, g_floor = 0.1)
s <- .emc_alloc_solve_units(u)
stopifnot(s$regime == 1, s$x[2] == 0)
emit("Regime-only pilot gate", regime = s$regime, supplemental_rate = s$x[2],
     calibration_wanted = s$x[2] > 0,
     intended_rate = min(1, s$regime + s$x[2]))

# 5. Mixed HMC eligibility changes the population averaged in the numerator.
u <- list(units = list(list(kind = "shear", rho = 0.9, g = 0.7, cost = 0.01),
  list(kind = "subject", rho = 0.8, h = 0.2, h_dev = 0.2,
       d_rho = 0.2, c_h = 0.01),
  list(kind = "subject", rho = 0.8, h = NULL, c_h = 0)),
  r_bar = 0.6, g_floor = 0.5)
projected <- .emc_alloc_regime(u, 0)$units[[1]]$rho
stopifnot(abs(projected - 0.9) > 0.1)
emit("Zero-rate regime must be identity", before = 0.9, after = projected)

# 6. Disabling curves only gates initial atlas creation.
options(emc2.shear_curve = FALSE)
z <- c(-1, -0.5, 0, 0.5, 1, 1.5)
ch <- list(n_pars = 1, samples = list(idx = 1L, theta_mu = matrix(0)),
  shear = list(plan = list(dirs = matrix(1), lam = 0,
                          last = list(n = c(0L, 0L))),
    atlas = list(nodes = list(list(mu = matrix(z, 1), lin = matrix(z, 1),
                                  prec = matrix(1, 1, length(z)))))))
curved <- .emc_shear_build_curves(list(ch), 1L)[[1]]
stopifnot(!is.null(curved$shear$curve))
emit("Curve switch with an existing atlas", option = getOption("emc2.shear_curve"),
     curve_built = !is.null(curved$shear$curve))

# 7. The trigger and allocator use different variance references.
vals <- c(-0.1, 0.1, -0.1, 0.1)
ch <- list(n_pars = 1, samples = list(theta_mu = matrix(vals, 1)),
  shear = list(plan = list(dirs = matrix(1), lam = 0.99,
                          last = list(n = c(10L, 0L), E = c(0.1, 0)))))
selected <- .emc_shear_curve_dirs(ch, 1:4)
stopifnot(!length(selected))
emit("Curve trigger variance", selected_directions = selected,
     g_short_own_window = 1 - (0.1 / 10) / (2 * var(vals)),
     g_pooled_variance_one = 1 - (0.1 / 10) / 2)

# 8. A q+r floor is an expected-count floor under independent selection.
emit("Refresh selection arithmetic", particle_rate = 0.5, HMC_rate = 0.5,
     probability_neither = (1 - 0.5) * (1 - 0.5))

# 9. Lag one alone can reward fast-mode changes that hurt ESS per cost.
weights <- c(0.2, 0.8)
tau_real <- function(lam) sum(weights * (1 + lam) / (1 - lam))
base_lam <- c(0.99, 0.1); candidate_lam <- c(0.99, 0)
base_pred <- .emc_alloc_tau(sum(weights * base_lam))
candidate_pred <- 1.05 * .emc_alloc_tau(sum(weights * candidate_lam))
base_real <- tau_real(base_lam); candidate_real <- 1.05 * tau_real(candidate_lam)
stopifnot(candidate_pred < 0.9 * base_pred, candidate_real > base_real)
emit("Two independent stationary Gaussian AR modes", rho1 = sum(weights * base_lam),
     predicted_tau_cost_base = base_pred, predicted_tau_cost_candidate = candidate_pred,
     real_tau_cost_base = base_real, real_tau_cost_candidate = candidate_real)
