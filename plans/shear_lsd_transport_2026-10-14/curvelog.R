# Logged sample-stage continuation: why curves are (not) built at each rebuild.
args <- commandArgs(TRUE); ckpt <- args[1]; iters <- as.integer(args[2]); logf <- args[3]; ncpc <- as.integer(args[4])
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
L <- function(fmt, ...) { a <- lapply(list(...), function(x) if (is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)) as.integer(x) else x)
  cat(do.call(sprintf, c(list(fmt), a)), "\n", file = logf, append = TRUE, sep = "") }
put <- function(n, f) { environment(f) <- ns; unlockBinding(n, ns); assign(n, f, ns); lockBinding(n, ns) }
orig <- list()
for (n in c(".emc_shear_curve_lsd", ".emc_shear_curve_dirs", ".emc_shear_curve_pilot", ".emc_shear_build_curves",
            ".emc_shear_curve_validate", ".emc_shear_atlas_subjects")) orig[[n]] <- get(n, ns)
gstat <- function(chain, idx, k_all, v) {
  rec <- chain$shear$plan$last; n <- rec$n[k_all]
  g <- 1 - (rec$E[k_all] / pmax(1L, n)) / (2 * pmax(v, .Machine$double.xmin)); list(n = n, g = g)
}
assign(".orig", orig, globalenv()); assign("L", L, globalenv()); assign("gstat", gstat, globalenv())
put(".emc_shear_curve_lsd", function(chain, idx, stage) {
  r <- get(".orig", globalenv())[[".emc_shear_curve_lsd"]](chain, idx, stage)
  rec <- chain$shear$plan$last
  if (identical(stage, "sample") && is.list(rec) && !is.null(rec$n)) {
    p <- chain$n_pars; v <- apply(.emc_shear_window_lsd(chain, idx), 1L, stats::var)
    s <- get("gstat", globalenv())(chain, idx, p + seq_len(p), v)
    nm <- rownames(chain$samples$theta_mu)
    get("L", globalenv())("LSD it=%d chosen=[%s] | g(n): %s", as.integer(chain$samples$idx), paste(nm[r], collapse = ","),
      paste(sprintf("%s=%.2f(%.0f)", nm, s$g, s$n), collapse = " "))
  }
  r
})
put(".emc_shear_curve_dirs", function(chain, idx) {
  r <- get(".orig", globalenv())[[".emc_shear_curve_dirs"]](chain, idx)
  plan <- chain$shear$plan
  if (is.list(plan$last) && !is.null(plan$lam)) {
    slow <- which(plan$lam >= .EMC_SHEAR_CURVE_LAM)
    Dinv <- tryCatch(solve(plan$dirs), error = function(e) NULL)
    if (!is.null(Dinv) && length(slow)) {
      t <- Dinv %*% chain$samples$theta_mu[, idx, drop = FALSE]; v <- apply(t, 1L, stats::var)
      s <- get("gstat", globalenv())(chain, idx, slow, v[slow])
      get("L", globalenv())("DIRS it=%d chosen=[%s] slow=%d | g(n): %s", chain$samples$idx, paste(r, collapse = ","), length(slow),
        paste(sprintf("d%d=%.2f(%.0f)", slow, s$g, s$n), collapse = " "))
    } else get("L", globalenv())("DIRS it=%d no slow dirs (lam max %.2f)", chain$samples$idx, max(plan$lam))
  }
  r
})
put(".emc_shear_curve_pilot", function(...) {
  r <- get(".orig", globalenv())[[".emc_shear_curve_pilot"]](...)
  get("L", globalenv())("PILOT gain=%.3f (need %.1f)", r, .EMC_SHEAR_CURVE_GAIN); r
})
put(".emc_shear_atlas_subjects", function(ch, idx, ks, n_cores) {
  r <- get(".orig", globalenv())[[".emc_shear_atlas_subjects"]](ch, idx, ks, n_cores)
  get("L", globalenv())("ATLAS it=%d ks=[%s] subjects fitted %d/%d", ch$samples$idx, paste(ks, collapse = ","), length(r), ch$n_subjects); r
})
put(".emc_shear_build_curves", function(emc, idx, stage) {
  out <- get(".orig", globalenv())[[".emc_shear_build_curves"]](emc, idx, stage)
  for (j in seq_along(out)) { ch <- out[[j]]; at <- ch$shear$atlas; cv <- ch$shear$curve
    nn <- if (is.null(at$nodes)) 0 else sum(vapply(at$nodes, function(a) !is.null(a), TRUE))
    get("L", globalenv())("BUILD chain %d it=%d stage=%s atlas_mean=%s subj_with_nodes=%d curve=%s", j, ch$samples$idx, stage, isTRUE(at$mean), nn,
      if (is.null(cv)) "none" else sprintf("q=[%s] l=[%s] kept_after_LOO+PD=%d", paste(cv$q, collapse = ","), paste(cv$l, collapse = ","),
        sum(!vapply(cv$B, is.null, TRUE)))) }
  out
})
put(".emc_shear_curve_validate", function(ch, idx, n_cores = 1L) {
  out <- get(".orig", globalenv())[[".emc_shear_curve_validate"]](ch, idx, n_cores)
  if (!is.null(ch$shear$curve)) {
    k0 <- sum(!vapply(ch$shear$curve$B, is.null, TRUE)); cv <- out$shear$curve
    get("L", globalenv())("VALIDATE it=%d kept %d -> %d", ch$samples$idx, k0, if (is.null(cv)) 0L else sum(!vapply(cv$B, is.null, TRUE)))
  }
  out
})
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(ckpt)), quiet = TRUE)
class(emc) <- "emc"
base <- local({ set.seed(1, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) { base <- parallel::nextRNGStream(base); emc[[j]]$rng <- ns$.emc_chain_rng(base, emc[[j]]$n_subjects); emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions }
L("rebuild schedule: next %s", paste(sapply(emc, function(c) c$sample_proposals_next), collapse = ","))
L("START %s iters %d", basename(ckpt), iters)
t0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample", stop_criteria = list(sample = list(iter = iters)),
  cores_for_chains = length(emc), cores_per_chain = ncpc, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)))
L("END elapsed %.0f", (proc.time() - t0)[3])
