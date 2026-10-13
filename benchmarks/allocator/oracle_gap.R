# How far the allocator's plan is from the best plan under its own model
# (oracle.R), on random small unit sets: one shear coordinate and two
# subjects with particle steps, HMC and deviation statistics.
#   EMC_LIB=<lib> Rscript benchmarks/allocator/oracle_gap.R [n_trials] [seed]
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "oracle.R"))
ns <- asNamespace("EMC2")
n_trials <- if (length(args) >= 1) as.integer(args[1]) else 40L
set.seed(if (length(args) >= 2) as.integer(args[2]) else 1L)

rand_problem <- function() {
  subj <- function() {
    d_rho <- runif(1, 0.1, 0.9)
    list(kind = "subject", rho = runif(1, 0.3, 0.98), p = runif(1, 0.1, 0.9),
         c_p = runif(1, 0.05, 0.4), h = runif(1, 0, 0.9), c_h = runif(1, 0.02, 0.4),
         d_rho = d_rho, p_dev = runif(1, 0.1, 0.9), h_dev = runif(1, 0.2, 1),
         q_old = 1, r_old = 0)
  }
  rho_g <- runif(1, 0.5, 0.99)
  list(units = list(list(kind = "shear", rho = rho_g, g = runif(1, 0.2, 0.95), cost = runif(1, 0.01, 0.2)),
                    subj(), subj()),
       C_rest = 1, joint = NULL, hmc_fixed = runif(1, 0, 0.2), shear_fixed = 0,
       g_floor = runif(1, 0, rho_g))
}

rows <- list()
for (t in seq_len(n_trials)) {
  u <- rand_problem()
  sol <- ns$.emc_alloc_solve_units(u)
  rates <- ns$.emc_alloc_plan_rates(u$units, sol)
  # A regime's rate is the floor of every subject's total HMC rate.
  reg <- if (is.null(sol$regime)) 0 else sol$regime
  for (i in 2:3) if (!is.null(u$units[[i]]$h)) rates$r[i] <- min(1, rates$r[i] + reg)
  # Under the model, the plan as run: group credited at its own rates.
  obj_sol <- ws1_model_objective(u, sol$x[1], rates$q[2:3], rates$r[2:3])
  obj0 <- ws1_model_objective(u, 0, c(1, 1), c(0, 0))
  orc <- ws1_oracle(u)
  rows[[t]] <- data.frame(trial = t, solver = obj_sol, claimed = sol$objective, zero = obj0,
                          oracle = orc$obj, gap = obj_sol / orc$obj,
                          oracle_gain = orc$obj / obj0, regime = sol$regime,
                          x = sol$x[1], q = paste(rates$q[2:3], collapse = "/"),
                          r = paste(round(rates$r[2:3], 2), collapse = "/"),
                          ox = orc$x[1], oq = paste(orc$q, collapse = "/"), or = paste(orc$r, collapse = "/"))
}
tab <- do.call(rbind, rows)
print(tab, digits = 3)
cat("\nclear-gain trials (oracle < 0.85 zero):", sum(tab$oracle_gain < 0.85), "\n")
cl <- tab[tab$oracle_gain < 0.85, ]
cat("gap quantiles:", round(quantile(cl$gap, c(0, .5, .9, 1)), 3), "\n")
cat("claimed vs model objective (solver honesty) quantiles:",
    round(quantile(tab$claimed / tab$solver, c(0, .5, 1)), 3), "\n")
