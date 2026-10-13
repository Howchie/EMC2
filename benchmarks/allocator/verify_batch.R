# Checks the batched subject pricing (.emc_alloc_units_need) against the
# per-unit path on every call made while solving captured rebuild inputs, and
# that the plans are those the per-unit path gives (saved by REF=1 runs).
#   EMC_LIB=<lib> Rscript benchmarks/allocator/verify_batch.R <out.rds> <alloc_*.rds files...>
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
out <- args[1]; files <- args[-1]
batched <- exists(".emc_alloc_units_need", envir = ns, inherits = FALSE)
stats <- new.env(); stats$calls <- 0L; stats$max_rel <- 0; stats$q_diff <- 0L; stats$fin_diff <- 0L
if (batched) {
  orig <- ns$.emc_alloc_units_need
  chk <- function(us, movable, rs, max_rounds, refresh = FALSE, layout = NULL) {
    new <- orig(us, movable, rs, max_rounds, refresh, layout)
    old <- lapply(seq_along(us), function(i)
      if (movable[i]) ns$.emc_alloc_unit_need(us[[i]], rs = rs, max_rounds = max_rounds, refresh = refresh)
      else ns$.emc_alloc_unit_still(us[[i]], length(rs)))
    stats$calls <- stats$calls + 1L
    for (i in seq_along(us)) {
      a <- old[[i]]; b <- new[[i]]
      fa <- is.finite(a$cost); fb <- is.finite(b$cost)
      stats$fin_diff <- stats$fin_diff + sum(fa != fb)
      both <- fa & fb
      stats$q_diff <- stats$q_diff + sum(a$q[both] != b$q[both])
      same_q <- both & a$q == b$q
      if (any(same_q)) {
        rel <- max(abs(a$cost[same_q] - b$cost[same_q]) / pmax(1e-300, abs(a$cost[same_q])),
                   abs(a$x[same_q] - b$x[same_q]) / pmax(1e-12, abs(a$x[same_q])))
        stats$max_rel <- max(stats$max_rel, rel)
      }
    }
    new
  }
  environment(chk) <- environment()
  assignInNamespace(".emc_alloc_units_need", chk, ns)
}
sols <- list()
t_solve <- 0
for (f in files) {
  cap <- readRDS(f)
  for (i in seq_along(cap$solve)) {
    t0 <- proc.time()[["elapsed"]]
    s <- ns$.emc_alloc_solve_units(cap$solve[[i]]$u)
    t_solve <- t_solve + proc.time()[["elapsed"]] - t0
    sols[[paste(f, i)]] <- s[c("x", "q", "joint", "regime", "tau", "objective")]
  }
}
cat("calls", stats$calls, "q differences", stats$q_diff, "reach differences", stats$fin_diff,
    "max relative difference", stats$max_rel, "\n")
saveRDS(sols, out)
