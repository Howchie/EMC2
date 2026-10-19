# Kernel-only RSUB cost (no likelihood plumbing) across the noise scale s,
# which sets the Skellam sum length for packets (terms ~ sqrt(a + b),
# a = 2 v u / s^2, b = 2 d / s^2).
# Run from the repo root:  Rscript plans/rsub_speed_2026-10-17/bench_kernel.R
lib <- Sys.getenv("EMC_LIB", ".lib/rigs")
suppressMessages(library(EMC2, lib.loc = lib))
set.seed(1)
n <- 1e5; u <- runif(n, .1, 3)
tm <- function(e) { t0 <- proc.time()[3]; for (i in 1:3) e(); unname(proc.time()[3] - t0) / 3 / n * 1e9 }
for (packets in c(FALSE, TRUE)) for (s in c(1, .3, .1)) {
  r <- c(pdf_A = tm(function() EMC2:::drsub(u, 1.3, .9, .5, s, packets, TRUE)),
         surv_A = tm(function() EMC2:::prsub(u, 1.3, .9, .5, s, packets, FALSE, TRUE)),
         pdf_A0 = tm(function() EMC2:::drsub(u, 1.3, .9, 0, s, packets, TRUE)),
         surv_A0 = tm(function() EMC2:::prsub(u, 1.3, .9, 0, s, packets, FALSE, TRUE)))
  cat(if (packets) "packets" else "dense  ", "s", format(s, width = 3), " ns/row:",
      paste(names(r), round(r), collapse = "  "), "\n")
}
