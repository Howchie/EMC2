# Direct comparison of the two things that disagree, with no estimator in
# between: the CMS simulator's empirical single-accumulator first-passage cdf
# against the PDE's.  The simulator draws exact stable increments, so its only
# error is discrete monitoring of the boundary (which biases passage late), and
# that is checked by refining dt.  The PDE side is already known to be converged
# and invariant here.  Scored at fixed probability levels so the comparison is
# in the units the likelihood cares about.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths()))  # RLFLIB: library to test
suppressMessages(library(EMC2))
source("refcdf.R")
hit <- EMC2:::rlf_hit_times_vec
B <- 1.5; A <- 0.5; b0 <- B + A
N <- 200000
qs <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

for (a in c(1.1, 1.5, 1.9)) for (v in c(1, 2)) {
  tg <- c(seq(1e-3, 10, length.out=5000), seq(10.02, 22, length.out=2000))
  C <- ref_cdf(a, v, tg, B, A)
  tq <- approx(C, tg, xout=qs, ties="ordered")$y
  cat(sprintf("\nalpha %.1f  v %.0f   (pde P(finish) = %.4f)\n", a, v,
              C[length(C)]))
  cat(sprintf("  %-10s %s\n", "quantile",
              paste(sprintf("%8.3f", qs), collapse="")))
  cat(sprintf("  %-10s %s\n", "pde t",
              paste(sprintf("%8.3f", tq), collapse="")))
  for (dt in c(1e-3, 2e-4)) {
    set.seed(7)
    # one accumulator, start uniform on [0, A]: rlf_hit_times_vec takes B and A
    ft <- hit(rep(v, N), rep(B, N), rep(A, N), rep(1, N), rep(a, N), dt, 22)
    emp <- sapply(tq, function(x) mean(ft <= x))
    cat(sprintf("  sim dt %-5.0e %s   (P(finish) %.4f)\n", dt,
                paste(sprintf("%8.4f", emp - qs), collapse=""),
                mean(is.finite(ft))))
    flush.console()
  }
}
cat("\ncdfcmpdone\n")
