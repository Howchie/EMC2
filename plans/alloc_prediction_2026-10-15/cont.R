# EMC_LIB=<lib> Rscript cont.R <label> <seed> <iters> <out.rds>
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(EMC2, lib.loc = Sys.getenv("EMC_LIB")); library(posterior) })
cat("EMC2 from", find.package("EMC2"), " shrink_g:", exists(".emc_alloc_shrink_g", asNamespace("EMC2")), "\n")
label <- args[1]; seed <- as.integer(args[2]); iters <- as.integer(args[3]); out <- args[4]
ns <- asNamespace("EMC2"); `%||%` <- function(a, b) if (is.null(a)) b else a
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE); class(emc) <- "emc"
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names; ns_sub <- emc[[1]]$n_subjects
group_draws <- function(e, idx) lapply(e, function(ch) cbind(t(ch$samples$theta_mu[, idx, drop = FALSE]),
  t(0.5 * log(apply(ch$samples$theta_var[, , idx, drop = FALSE], 3, diag)))))
n0 <- emc[[1]]$samples$idx
ref <- (n0 - 399):n0                          # last 400 checkpoint draws (burn + adapt)
G <- group_draws(emc, ref)
centre <- Reduce(`+`, lapply(G, colMeans)) / 3; scale <- sqrt(Reduce(`+`, lapply(G, function(g) apply(g, 2, var))) / 3)
base <- local({ set.seed(seed, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) { base <- parallel::nextRNGStream(base); emc[[j]]$rng <- ns$.emc_chain_rng(base, ns_sub)
  emc[[j]]$sample_proposals_next <- emc[[j]]$sample_transitions }
totals <- function(e) lapply(e, function(ch) c(hmc = ch$subject_hmc$totals$attempts %||% 0,
  shear = ch$shear$rounds %||% 0, collapsed = ch$collapsed$rounds %||% 0))
tot0 <- totals(emc)
cg_file <- file.path("/sys/fs/cgroup", sub("^0::", "", readLines("/proc/self/cgroup")[1]), "cpu.stat")
cg_cpu <- function() { x <- tryCatch(readLines(cg_file), error = function(e) character(0)); u <- grep("^usage_usec", x, value = TRUE)
  if (length(u)) as.numeric(sub("usage_usec ", "", u)) / 1e6 else NA_real_ }
c0 <- cg_cpu(); t0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample", stop_criteria = list(sample = list(iter = iters)),
  step_size = 100L, max_tries = 1L, cores_for_chains = 3, cores_per_chain = 9, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)))
el <- proc.time() - t0; cpu <- cg_cpu() - c0
if (!is.finite(cpu)) cpu <- sum(el[c("user.self", "user.child", "sys.self", "sys.child")])
emc <- ns$restore_duplicates(emc)
new <- (n0 + 1):emc[[1]]$samples$idx; score <- new[-(1:100)]
Gn <- group_draws(emc, new); Gs <- lapply(Gn, function(g) g[-(1:100), , drop = FALSE])
arr <- aperm(simplify2array(Gs), c(1, 3, 2)); colnames_q <- c(paste0("mu.", pn), paste0("lsd.", pn))
ess <- setNames(apply(arr, 3, ess_bulk), colnames_q); rh <- setNames(apply(arr, 3, rhat), colnames_q)
lag1 <- setNames(apply(arr, 3, function(a) mean(apply(a, 2, function(x) acf(x, 1, plot = FALSE)$acf[2]))), colnames_q)
# Mean squared lag-10 jump, standardised by the checkpoint variance (all new draws).
jump10 <- setNames(Reduce(`+`, lapply(Gn, function(g) { n <- nrow(g); colMeans((g[-(1:10), ] - g[-((n - 9):n), ])^2) })) / 3 / scale^2, colnames_q)
ridge <- paste0("mu.", c("v_LoadHigh", "v_LoadLow", "K", "t0_LoadHigh", "t0_LoadLow", "sv"))
res <- list(label = label, seed = seed, lib = find.package("EMC2"), cpu = cpu, elapsed = el[["elapsed"]], ess = ess, rhat = rh,
  lag1 = lag1, jump10 = jump10, moves = mapply(function(a, b) b - a, tot0, totals(emc)),
  plans = lapply(emc, function(ch) tail(ch$alloc$history, 4)), draws = Gn)
saveRDS(res, out)
cat(sprintf("%s seed %d: cpu %.0f s elapsed %.0f s\n", label, seed, cpu, el[["elapsed"]]))
print(round(rbind(ess = ess[ridge], per_kcpu = 1000 * ess[ridge] / cpu, lag1 = lag1[ridge], jump10 = jump10[ridge], rhat = rh[ridge]), 3))
cat(sprintf("all group: ESS min %.0f median %.0f  min/kCPU %.2f  max rhat %.3f\n", min(ess), median(ess), 1000 * min(ess) / cpu, max(rh)))
print(res$moves)
