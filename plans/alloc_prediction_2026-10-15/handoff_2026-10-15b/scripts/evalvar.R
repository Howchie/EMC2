# Which single-block statistic, with which variance, predicts a group coordinate's realized tau?
# For each arm: realized tau over draws 101-400 (Geyer, 3 chains, pooled variance).
# Predictors from block b (draws (b-1)*100+1 .. b*100) with variance V_blk (block, pooled over chains),
# V_stage (draws 1..b*100), V_ck (last 400 checkpoint draws + stage so far).
suppressPackageStartupMessages(library(posterior))
A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15"
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
gd <- function(e, idx) lapply(e, function(ch) cbind(t(ch$samples$theta_mu[, idx, drop = FALSE]),
  t(0.5 * log(apply(ch$samples$theta_var[, , idx, drop = FALSE], 3, diag)))))
n0 <- emc[[1]]$samples$idx
CK <- gd(emc, (n0 - 399):n0)
cn <- c(paste0("mu.", pn), paste0("lsd.", pn))
pooled_var <- function(X) apply(do.call(rbind, X), 2, var)
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-9); (1 + r) / (1 - r) }
# mean squared lag-L jump per iteration, chain-averaged
msj <- function(X, L) Reduce(`+`, lapply(X, function(g) colMeans((g[-(1:L), , drop = FALSE] - g[1:(nrow(g) - L), , drop = FALSE])^2))) / length(X)
# transport rho: AR(1) whose lag-L jump matches
trho <- function(T, V, L) pmax(0, 1 - T / (2 * V) )^(1 / L)
geyer_tau <- function(X) { a <- simplify2array(X); apply(a, 2, function(x) nrow(x) * ncol(x) / ess_basic(x)) }
files <- c(old_1 = file.path(A, "arms/old_1.rds"), new_1 = file.path(A, "arms/new_1.rds"),
  new_2 = file.path(A, "handoff_2026-10-15/arms/pair_new_2.rds"), new_3 = file.path(A, "handoff_2026-10-15/arms/pair_new_3.rds"),
  point_1 = file.path(A, "handoff_2026-10-15/arms/arm_point_1.rds"), point_2 = file.path(A, "handoff_2026-10-15/arms/pair_point_2.rds"),
  point_3 = file.path(A, "handoff_2026-10-15/arms/pair_point_3.rds"), kappa_1 = file.path(A, "handoff_2026-10-15/arms/arm_kappa_1.rds"),
  transport_1 = file.path(A, "eval_transport/arm_transport_1.rds"), slow2_1 = file.path(A, "arms/slow2_1.rds"))
res <- list(); allrows <- NULL
for (nm in names(files)) {
  r <- readRDS(files[nm]); G <- r$draws
  run <- geyer_tau(lapply(G, function(g) g[101:400, ])); names(run) <- cn
  slow <- order(-run)[1:8]
  rows <- NULL
  for (b in 1:3) {
    idx <- (b - 1) * 100 + 1:100
    blk <- lapply(G, function(g) g[idx, , drop = FALSE])
    Vb <- pooled_var(blk)
    Vs <- pooled_var(lapply(G, function(g) g[1:(b * 100), , drop = FALSE]))
    Vck <- pooled_var(mapply(function(c, g) rbind(c, g[1:(b * 100), , drop = FALSE]), CK, G, SIMPLIFY = FALSE))
    Vck0 <- pooled_var(CK)
    T1 <- msj(blk, 1); T10 <- msj(blk, 10); T25 <- msj(blk, 25)
    means <- sapply(blk, colMeans); W <- Reduce(`+`, lapply(blk, function(g) apply(g, 2, var))) / 3
    Bm <- apply(means, 1, var)
    est <- list(
      lag1_blk = tau1(1 - T1 / (2 * Vb)), lag1_stage = tau1(1 - T1 / (2 * Vs)), lag1_ck = tau1(1 - T1 / (2 * Vck)),
      L10_blk = tau1(trho(T10, Vb, 10)), L10_stage = tau1(trho(T10, Vs, 10)), L10_ck = tau1(trho(T10, Vck, 10)),
      L25_blk = tau1(trho(T25, Vb, 25)), L25_stage = tau1(trho(T25, Vs, 25)), L25_ck = tau1(trho(T25, Vck, 25)),
      bm = pmax(1, 100 * Bm / W),
      # 'drift' tau: the block's chain-mean spread against the lag-1 jump: tau = 2 (W + B) / T1 with B counted n times? no: the AR(1) at block var plus n B / V_b
      stageG = if (b >= 1) geyer_tau(lapply(G, function(g) g[1:(b * 100), , drop = FALSE])) else NA)
    rows <- rbind(rows, data.frame(arm = nm, blk = b, coord = cn, run = run, slow8 = seq_along(cn) %in% slow,
      Vratio_stage = Vs / Vb, Vratio_ck = Vck / Vb, Vratio_ck0 = Vck0 / Vb,
      as.data.frame(lapply(est, function(e) as.numeric(e)))))
  }
  allrows <- rbind(allrows, rows)
}
saveRDS(allrows, file.path(Sys.getenv("SP"), "evalvar.rds"))
ests <- c("lag1_blk", "lag1_stage", "lag1_ck", "L10_blk", "L10_stage", "L10_ck", "L25_blk", "L25_stage", "L25_ck", "bm", "stageG")
summ <- function(d) sapply(ests, function(e) { x <- d[[e]] / d$run; s <- d$slow8
  c(med_all = median(x, na.rm = TRUE), med_slow8 = median(x[s], na.rm = TRUE), q10_slow8 = quantile(x[s], .1, na.rm = TRUE), q90_slow8 = quantile(x[s], .9, na.rm = TRUE),
    spearman = cor(d[[e]], d$run, method = "spearman", use = "complete")) })
cat("== predicted / realized tau, by block (all arms pooled)\n")
for (b in 1:3) { cat("block", b, "\n"); print(round(t(summ(allrows[allrows$blk == b, ])), 2)) }
cat("\n== variance ratios on the slowest 8 (median over arms), by block\n")
print(aggregate(cbind(Vratio_stage, Vratio_ck, Vratio_ck0) ~ blk, allrows[allrows$slow8, ], median))
cat("\n== realized tau slowest 8 per arm\n")
for (nm in names(files)) { d <- allrows[allrows$arm == nm & allrows$blk == 1 & allrows$slow8, ]; cat(sprintf("%-12s %s\n", nm, paste(sprintf("%s %.0f", d$coord, d$run), collapse = " "))) }
