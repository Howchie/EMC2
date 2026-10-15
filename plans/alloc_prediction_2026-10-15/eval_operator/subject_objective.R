# Experimental decomposition: subject budgets are priced independently of
# the group's target; their exact combined cost enters the group solve.
subject_objective <- function(u, ev, mode = c("own", "baseline")) {
  mode <- match.arg(mode)
  units <- u$units
  is_subject <- vapply(units, function(x) identical(x$kind, "subject"), logical(1))
  si <- which(is_subject)
  gi <- which(!is_subject)
  n <- length(si)
  hmc_fixed <- if (is.null(u$hmc_fixed)) 0 else u$hmc_fixed
  alloc <- lapply(si, function(i) {
    x <- units[[i]]
    if (mode == "baseline") {
      need <- ev$.emc_alloc_unit_need(x, max(ev$.emc_alloc_subject_base(x)), ev$.EMC_ALLOC_MAX_ROUNDS)
      return(list(x = need$x, q = need$q, cost = need$cost))
    }
    sol <- ev$.emc_alloc_solve(list(x), u$C_rest / max(1L, n), hmc_fixed = hmc_fixed / max(1L, n))
    list(x = sol$x, q = sol$q,
         cost = sol$q * x$c_p + sol$x * x$c_h)
  })
  sx <- vapply(alloc, function(x) x$x, numeric(1))
  sq <- vapply(alloc, function(x) x$q, numeric(1))
  sc <- sum(vapply(alloc, function(x) x$cost, numeric(1)))
  activation <- hmc_fixed * (1 - prod(1 - sx))
  joint <- u$joint
  if (!is.null(joint)) joint$idx <- match(joint$idx, gi)
  sol <- ev$.emc_alloc_solve(units[gi], u$C_rest + sc + activation, joint = joint,
                            shear_fixed = if (is.null(u$shear_fixed)) 0 else u$shear_fixed)
  gx <- sol$x; gq <- sol$q
  sol$x <- numeric(length(units)); sol$q <- rep(1, length(units))
  sol$x[gi] <- gx; sol$x[si] <- sx; sol$q[gi] <- gq; sol$q[si] <- sq
  sol$regime <- 0
  sol$choice <- "0"
  sol$options <- c(none = Inf, "0" = sol$objective)
  sol$pred <- ev$.emc_alloc_predict(units, sol, u$joint)
  sol
}
