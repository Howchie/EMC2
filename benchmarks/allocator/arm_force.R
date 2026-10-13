# Counterfactual arm for calibrate.R (ARM_R): every production plan is
# forced away from the solver's choice, and the recorded predictions follow
# the forced plan, so analyse.R scores the model's factors off its optimum.
#   FORCE=up    shear coordinates the plan left at zero rounds get UP_ROUNDS
#               (default 1); subjects with measured HMC and no HMC get rate
#               UP_RATE (default 0.5) on top of any regime;
#   FORCE=down  shear rounds and subject HMC rates above zero are divided by
#               DOWN_DIV (default 3).
# The particle rates and the collapsed rate stay as solved.

arm_hook <- function(emc) {
  ns <- asNamespace("EMC2")
  mode <- Sys.getenv("FORCE", "up")
  up_rounds <- as.numeric(Sys.getenv("UP_ROUNDS", "1"))
  up_rate <- as.numeric(Sys.getenv("UP_RATE", "0.5"))
  div <- as.numeric(Sys.getenv("DOWN_DIV", "3"))
  orig_units <- get(".emc_alloc_units", envir = ns)
  units_f <- function(ms, hmc_state, shear_on, optimistic, B = NULL) {
    u <- orig_units(ms, hmc_state, shear_on, optimistic, B)
    attr(u, "ws1_production") <- !isTRUE(optimistic)
    u
  }
  environment(units_f) <- list2env(list(orig_units = orig_units), parent = environment(orig_units))
  utils::assignInNamespace(".emc_alloc_units", units_f, ns)

  orig_solve <- get(".emc_alloc_solve_units", envir = ns)
  solve_f <- function(u) {
    depth <<- depth + 1L
    on.exit(depth <<- depth - 1L)
    sol <- orig_solve(u)
    # Only the plan a chain runs: the outermost solve of a production unit set.
    if (depth != 1L || !isTRUE(attr(u, "ws1_production"))) return(sol)
    units <- u$units
    reg <- if (is.null(sol$regime)) 0 else sol$regime
    x <- sol$x
    for (i in seq_along(units)) {
      z <- units[[i]]
      if (identical(z$kind, "shear") && is.finite(z$g) && z$g < 1) {
        x_new <- if (mode == "up") { if (x[i] == 0) up_rounds else x[i] } else x[i] / div
        if (x_new != x[i]) {
          a <- if (is.null(z$a)) z$g else z$a
          sol$pred[[i]] <- sol$pred[[i]] * ns$.emc_alloc_factor(z$g, x_new, a) /
            ns$.emc_alloc_factor(z$g, x[i], a)
          x[i] <- x_new
        }
      } else if (identical(z$kind, "subject") && !is.null(z$h)) {
        r_old <- min(1, reg + x[i])
        r_new <- if (mode == "up") { if (r_old == 0) up_rate else r_old } else r_old / div
        if (r_new != r_old) {
          q <- if (is.null(z$p) || length(sol$q) < i) 1 else sol$q[i]
          P <- if (is.null(z$p)) rep(1, length(z$rho)) else z$p
          f <- function(r) vapply(seq_along(z$rho), function(k)
            ns$.emc_alloc_dev_factor(q, r, P[k], z$h[k]), numeric(1))
          sol$pred[[i]] <- sol$pred[[i]] * f(r_new) / pmax(f(r_old), 1e-12)
          # The chain runs regime + x, capped at one.
          x[i] <- r_new - reg
        }
      }
    }
    sol$x <- x
    sol
  }
  environment(solve_f) <- list2env(list(orig_solve = orig_solve, mode = mode, depth = 0L,
                                        up_rounds = up_rounds, up_rate = up_rate, div = div),
                                   parent = environment(orig_solve))
  utils::assignInNamespace(".emc_alloc_solve_units", solve_f, ns)
  cat("arm: FORCE =", mode, "\n")
  emc
}
