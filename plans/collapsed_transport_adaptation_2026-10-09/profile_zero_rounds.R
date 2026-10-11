pkgload::load_all(".", quiet = TRUE, compile = FALSE)
ns <- asNamespace("EMC2")
before <- new.env(parent = ns)
sys.source("/tmp/collapsed_followup_before/collapsed_move.R", envir = before)
x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/calibrated_snapshot.rds")), quiet = TRUE)
ch <- x[[1L]]
i <- ch$samples$idx
ch$collapsed$plan$rounds <- 0.125
set.seed(sum(utf8ToInt("collapsed-zero-round-profile")))
repeat {
  seed <- .Random.seed
  if (runif(1L) >= ch$collapsed$plan$rounds) break
}
ch$rng$gibbs <- seed
pars <- list(tmu = ch$samples$theta_mu[, i], tvar = ch$samples$theta_var[, , i],
             tvinv = solve(ch$samples$theta_var[, , i]))
proposals <- rbind(ch$samples$alpha[, , i], ch$samples$subj_ll[, i])
a <- before$.emc_collapsed_move(ch, pars, proposals, "sample")
b <- EMC2:::.emc_collapsed_move(ch, pars, proposals, "sample")
stopifnot(identical(a, b))
bench <- function(move, n = 100L) {
  started <- proc.time()[["elapsed"]]
  for (k in seq_len(n)) move(ch, pars, proposals, "sample")
  (proc.time()[["elapsed"]] - started) / n
}
cat("identical zero-round output, RNG and counters\n")
cat("previous seconds/call", bench(before$.emc_collapsed_move), "\n")
cat("new seconds/call", bench(EMC2:::.emc_collapsed_move), "\n")
