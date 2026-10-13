# A fixed move allocation in place of the measured one: `shear` rounds on
# every shear coordinate, `collapsed` rounds of the collapsed move, and every
# subject with HMC an HMC rate of `hmc` and a particle rate of one minus it,
# whatever the measurements say. Plans then do not depend on timings.
force_allocation <- function(shear = 1, collapsed = 1, hmc = 0.5, env = parent.frame()) {
  testthat::local_mocked_bindings(.emc_alloc_solve_units = function(u) {
    units <- u$units
    has_h <- vapply(units, function(x) identical(x$kind, "subject") && !is.null(x$h),
                    logical(1))
    x <- vapply(units, function(x) if (identical(x$kind, "shear")) shear else 0, numeric(1))
    x[has_h] <- hmc
    list(x = x, q = ifelse(has_h, 1 - hmc, 1),
         joint = if (!is.null(u$joint)) collapsed,
         rho = NA_real_, tau = NA_real_, cost = u$C_rest,
         objective = NA_real_, objective0 = NA_real_, tau0 = NA_real_)
  }, .package = "EMC2", .env = env)
}

# A smaller ensemble update (pool size and rounds) to keep a test fast.
small_ensemble <- function(pool, rounds = 10L, env = parent.frame()) {
  cfg <- EMC2:::.emc_ensemble_options()
  cfg$pool <- as.integer(pool)
  cfg$rounds <- as.integer(rounds)
  testthat::local_mocked_bindings(.emc_ensemble_options = function() cfg,
                                  .package = "EMC2", .env = env)
}
