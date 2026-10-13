# Functional check of the between-block plan on a short Forstmann fit.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
des <- design(data = forstmann, model = LNR, matchfun = function(d) d$S == d$lR,
              formula = list(m ~ lM, s ~ 1, t0 ~ 1), verbose = FALSE)
set.seed(1)
emc <- suppressMessages(make_emc(forstmann, des, n_chains = 1))
emc <- suppressWarnings(suppressMessages(fit(emc, cores_per_chain = as.integer(Sys.getenv("CORES", "4")),
                                             cores_for_chains = 1, verbose = FALSE)))
h <- emc[[1]]$ensemble_plan_history
cat(sprintf("%d plan updates\n", length(h)))
for (r in h) cat(sprintf("n %3d M %3d R %2d -> M %3d R %2d | G %.3f J1 %.3f JR %.3f J_tot %.3f C %.4f c_round %.2e c_member %.2e eff %.3f moved %.3f\n",
  r$n, r$M, r$R, r$pool, r$rounds, r$G, r$J[1], r$J[length(r$J)], r$J_tot, r$C, r$c_round, r$c_member, r$eff, r$moved))
print(ns$ensemble_diagnostics(emc))
saveRDS(emc, Sys.getenv("OUT", "/data/work/EMC2_ws5_libs/prof/plan_smoke_emc.rds"))
