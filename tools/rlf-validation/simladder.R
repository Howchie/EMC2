# cdfcmp showed simulator-vs-solver cdf agreement of ~0.002 at alpha 1.5-1.9
# once emc2.rlf_sim_dt is refined, but disagreement up to 0.083 at alpha = 1.1
# that did not improve from dt 1e-3 to 2e-4 -- and the two dt values disagreed
# with each other by 0.03, i.e. the simulator itself is not converged there.
# Ladder dt down and see whether it converges onto the solver (simulator at
# fault) or onto something else (formulation difference).
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths())); suppressMessages(library(EMC2))
source("refcdf.R")
hit <- EMC2:::rlf_hit_times_vec
B <- 1.5; A <- 0.5; N <- 50000
qs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
tg <- c(seq(1e-3, 10, length.out=5000), seq(10.02, 22, length.out=2000))
for (a in c(1.1, 1.3)) for (v in c(2)) {
  C <- ref_cdf(a, v, tg, B, A)
  tq <- approx(C, tg, xout=qs, ties="ordered")$y
  cat(sprintf("\nalpha %.1f v %.0f  (mc se at q=.25 is %.4f)\n", a, v,
              sqrt(0.25*0.75/N)))
  cat(sprintf("  %-12s %s\n", "quantile", paste(sprintf("%9.2f", qs), collapse="")))
  for (dt in c(1e-3, 2e-4, 4e-5)) {
    set.seed(7)
    ft <- hit(rep(v,N), rep(B,N), rep(A,N), rep(1,N), rep(a,N), dt, 22)
    cat(sprintf("  sim dt %-6.0e %s\n", dt,
        paste(sprintf("%9.4f", sapply(tq, function(x) mean(ft <= x)) - qs),
              collapse=""))); flush.console()
  }
}
cat("\nladderdone\n")
