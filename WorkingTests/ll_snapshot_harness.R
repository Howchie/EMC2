# ll snapshot harness for Phase-2 refactor verification (plan.md commits 9-12).
# Usage:
#   Rscript ll_harness.R capture <out.rds>     # build datasets, compute lls, save
#   Rscript ll_harness.R compare <ref.rds> <new.rds> [tol]
# 'capture' regenerates datasets from fixed seeds, so a later capture on the
# same seeds + compare against the stored reference detects any ll change.
#
# CAUTION: library(EMC2) loads whatever install is first on .libPaths(), which
# is easy to leave stale — a stale install makes ref-vs-new comparisons pass
# vacuously (both captures run the same old code). Install the code under test
# into a dedicated library and run with R_LIBS=<that lib>; the banner below
# prints the resolved package path and .so mtime so captures are auditable.

suppressMessages({library(EMC2); library(dplyr)})
cat(sprintf("EMC2 from %s (.so mtime %s)\n", find.package("EMC2"),
            format(file.info(file.path(find.package("EMC2"), "libs",
                                       paste0("EMC2", .Platform$dynlib.ext)))$mtime)))

N_TRIALS   <- 300
N_PART     <- 100
PART_SD    <- 0.25

.get_designs <- function(dadm, p_types) {
  setNames(lapply(p_types, function(p) {
    d <- attr(dadm, "designs")[[p]]
    d[attr(d, "expand"), , drop = FALSE]
  }), p_types)
}

cpp_ll_mat <- function(p_mat, dadm, model) {
  p_types <- names(model$p_types)
  consts  <- attr(dadm, "constants"); if (is.null(consts)) consts <- NA
  EMC2:::calc_ll_oo(p_mat, dadm, constants = consts,
                    designs   = .get_designs(dadm, p_types),
                    type      = model$c_name,
                    bound     = model$bound,
                    transform = model$transform,
                    pretransforms = model$pre_transform,
                    p_types   = p_types,
                    min_ll    = log(1e-10),
                    trend     = model$trend)
}

make_p_mat <- function(p_vector, seed) {
  set.seed(seed)
  n <- length(p_vector)
  p_mat <- matrix(rep(as.numeric(p_vector), each = N_PART), nrow = N_PART) +
    matrix(rnorm(N_PART * n, sd = PART_SD), nrow = N_PART)
  colnames(p_mat) <- names(p_vector)
  p_mat
}

run_case <- function(name, design_obj, p_vector, seed, n_trials = N_TRIALS, TC = NULL,
                     make_data_args = list()) {
  set.seed(seed)
  args <- c(list(p_vector, design_obj, n_trials = n_trials), make_data_args)
  if (!is.null(TC)) args$TC <- TC
  dat <- do.call(make_data, args)
  emc <- make_emc(dat, design_obj, type = "single", verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_mat <- make_p_mat(p_vector, seed + 1)
  t0 <- proc.time()[["elapsed"]]
  lls <- as.numeric(cpp_ll_mat(p_mat, dadm, model))
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("%-22s n_ll=%d finite=%d time=%.2fs\n", name, length(lls), sum(is.finite(lls)), el))
  list(lls = lls, elapsed = el)
}

capture_all <- function() {
  RNGkind("L'Ecuyer-CMRG")
  out <- list()

  ## --- DDM, censored + truncated -------------------------------------------
  designDDM <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"), model = DDM,
    formula = list(v ~ 0 + S, a ~ 1, t0 ~ 1, Z ~ 1),
    constants = c(s = log(1), sv = log(0), SZ = qnorm(0), st0 = log(0)))
  pv <- sampled_pars(designDDM, doMap = FALSE)
  pv[c("v_Sleft","v_Sright","a","t0","Z")] <- c(-1.2, 1.2, log(2), log(0.2), qnorm(0.6))
  out$DDM <- run_case("DDM cens+trunc", designDDM, pv, seed = 101,
                      TC = list(UC = 1.35, LC = 0.55, LT = 0.3, UT = 1.8))

  ## --- LBA, censored + truncated -------------------------------------------
  mf2 <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = mf2, model = LBA(posdrift = TRUE),
    formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ 1),
    constants = c(A = log(0.5)))
  pv <- sampled_pars(designLBA, doMap = FALSE)
  pv[c("B","v","v_lMTRUE","t0","sv")] <- c(log(1.5), .5, .5, log(0.2), log(1))
  out$LBA <- run_case("LBA cens+trunc", designLBA, pv, seed = 102,
                      TC = list(UC = 1.8, LC = .925, UT = 1.9, LT = .85))

  ## --- LBAIO (posdrift = FALSE), censored ----------------------------------
  designLBAIO <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = mf2, model = LBA(posdrift = FALSE),
    formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, sv ~ 1),
    constants = c(A = log(0.5)))
  pv <- sampled_pars(designLBAIO, doMap = FALSE)
  pv[c("B","v","v_lMTRUE","t0","sv")] <- c(log(1.5), .5, .5, log(0.2), log(1))
  out$LBAIO <- run_case("LBAIO cens", designLBAIO, pv, seed = 103,
                        TC = list(UC = 3, LC = 0.8))

  ## --- RDM + pContaminant, censored ----------------------------------------
  designRDMpc <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = mf2, model = RDM,
    formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM, pContaminant ~ 1),
    constants = c(s = log(1)))
  pv <- sampled_pars(designRDMpc, doMap = FALSE)
  pv[c("B","v","v_lMTRUE","A","t0","s_lMTRUE","pContaminant")] <-
    c(log(2), log(1), log(2), log(.5), log(0.2), log(.8), qnorm(0.15))
  out$RDMpC <- run_case("RDM pContaminant", designRDMpc, pv, seed = 104,
                        TC = list(UC = 1.775, LC = .825))

  ## --- RDMSWTN local_kill, sv > 0 (SWTN quadrature path) -------------------
  designSWTN <- design(
    factors = list(S = "Target", subjects = 1), Rlevels = c("Go"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(erlang_type = "local_kill", erlang_shape = 2), UC = 3)
  pv <- sampled_pars(designSWTN, doMap = FALSE)
  pv[c("v","B","A","t0","sv","mK")] <- c(log(3), log(1.1), log(0.3), log(0.2), log(0.5), log(2))
  out$RDMSWTN_local <- run_case("RDMSWTN local sv>0", designSWTN, pv, seed = 105)

  ## --- RDMSWTN global_kill (global-omission GL20 path) ---------------------
  designSWTNg <- design(
    factors = list(S = "Target", subjects = 1), Rlevels = c("Go"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
    constants = c(s = log(1), sv = log(0), A = log(0)),
    model = RDMSWTN(erlang_type = "global_kill", erlang_shape = 1), UC = 3)
  pv <- sampled_pars(designSWTNg, doMap = FALSE)
  pv[c("v","B","t0","mK")] <- c(log(3), log(1.1), log(0.2), log(2))
  out$RDMSWTN_global <- run_case("RDMSWTN global kill", designSWTNg, pv, seed = 106)

  ## --- BAwL local_kill ------------------------------------------------------
  designBAwL <- design(
    factors = list(S = "Target", subjects = 1), Rlevels = c("Go"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1, k ~ 1, mK ~ 1),
    constants = c(sv = log(1)),
    model = BAwL(erlang_type = "local_kill", erlang_shape = 1), UC = 3)
  pv <- sampled_pars(designBAwL, doMap = FALSE)
  pv[c("v","B","A","t0","k","mK")] <- c(1.5, log(1.1), log(0.3), log(0.2), log(0.5), log(2))
  out$BAwL <- run_case("BAwL local kill", designBAwL, pv, seed = 107)

  ## --- LogicalRules LBA (31-node GL pass) ----------------------------------
  funcs <- list(mismatch = function(d) ifelse(d$lM == TRUE, 0, 1),
                MatchY = function(d) ifelse(d$lM == TRUE & (d$lR == "A" | d$lR == "B"), 1, 0),
                MatchN = function(d) ifelse(d$lM == TRUE & (d$lR == "n_A" | d$lR == "n_B"), 1, 0),
                DT = function(d) ifelse(d$lM == TRUE & (d$lR == "A" | d$lR == "B") & d$S == "AB", 1, 0),
                Yes = function(d) ifelse(d$lR == "A" | d$lR == "B", 1, 0),
                No = function(d) ifelse(d$lR == "n_A" | d$lR == "n_B", 1, 0))
  mfLR <- function(d) dplyr::case_when(
    d$S == "NN" & d$lR == "n_A" ~ TRUE, d$S == "NN" & d$lR == "n_B" ~ TRUE,
    d$S == "AN" & d$lR == "A" ~ TRUE,  d$S == "AN" & d$lR == "n_B" ~ TRUE,
    d$S == "NB" & d$lR == "n_A" ~ TRUE, d$S == "NB" & d$lR == "B" ~ TRUE,
    d$S == "AB" & d$lR == "A" ~ TRUE,  d$S == "AB" & d$lR == "B" ~ TRUE,
    TRUE ~ FALSE)
  designLR <- design(
    fixed_accumulator_roles = factor(c("A","n_A","B","n_B"), levels = c("A","n_A","B","n_B")),
    matchfun = mfLR, model = LogicalRulesLBA, constants = c(sv = log(1)),
    factors = list(subjects = 1,
                   S = factor(rep(c("AB","AN","NB","NN"), 2), levels = c("AB","AN","NB","NN")),
                   LogicalRule = "OR"),
    Rlevels = c("yes","no"),
    formula = list(v ~ 0 + mismatch + MatchY + MatchN + DT, B ~ 0 + Yes + No, t0 ~ 1, A ~ 1, sv ~ 1),
    functions = funcs, UT = Inf)
  pv <- sampled_pars(designLR, doMap = FALSE)
  pv[c("v_mismatch","v_MatchY","v_MatchN","v_DT")] <- c(1, 2.5, 2.8, -0.25)
  pv[c("B_Yes","B_No","t0","A")] <- c(log(1), log(1.2), log(.2), log(.4))
  out$LogicalRules <- run_case("LogicalRules OR", designLR, pv, seed = 108, n_trials = 100)

  ## --- SSEXG no-ST, UC + LC censoring, staircase --------------------------
  lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)
  designSSE <- design(
    model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_go),
    covariates = "SSD",
    formula = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSE, doMap = FALSE)
  pv[c("mu_lMFALSE","mu_lMTRUE","sigma","tau","muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.65), log(0.50), log(0.05), log(0.15), log(0.20), log(0.04), log(0.08),
      qnorm(0.08), qnorm(0.08))
  out$SSEXG <- run_case("SSEXG stair UC+LC", designSSE, pv, seed = 109, n_trials = 400,
                        TC = list(UC = 1.2, LC = 0.3),
                        make_data_args = list(functions = list(SSD = make_ssd())))

  ## --- SSEXG with ST accumulator, UC --------------------------------------
  lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
  AccTypeFun <- function(d) {
    is_st <- as.character(d$lR) == "st"
    is_match <- !is_st & (as.character(d$S) == as.character(d$lR))
    factor(ifelse(is_st, "st", ifelse(is_match, "match", "mismatch")),
           levels = c("mismatch", "match", "st"))
  }
  mySSD <- function(d) SSD_function(d, SSD = c(0.20, 0.35), pSSD = c(0.25, 0.25))
  designSSEst <- design(
    model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right", "st"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_st, SSD = mySSD, AccType = AccTypeFun),
    formula = list(mu ~ 0 + AccType, sigma ~ 1, tau ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSEst, doMap = FALSE)
  pv[c("mu_AccTypemismatch","mu_AccTypematch","mu_AccTypest","sigma","tau",
       "muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.65), log(0.50), log(0.38), log(0.05), log(0.15),
      log(0.20), log(0.04), log(0.08), qnorm(0.05), qnorm(0.01))
  out$SSEXG_ST <- run_case("SSEXG ST UC", designSSEst, pv, seed = 110, n_trials = 400,
                           TC = list(UC = 1.2))

  ## --- SSRDEX no-ST, UC ----------------------------------------------------
  designSSR <- design(
    model = SSRDEX, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_go, SSD = mySSD),
    formula = list(v ~ 0 + lM, B ~ 1, A ~ 1, t0 ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSR, doMap = FALSE)
  pv[c("v_lMFALSE","v_lMTRUE","B","A","t0","muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.2), log(1.5), log(1), log(0.4), log(0.15),
      log(0.2), log(0.08), log(0.1), qnorm(0.08), qnorm(0.08))
  out$SSRDEX <- run_case("SSRDEX UC", designSSR, pv, seed = 111, n_trials = 400,
                         TC = list(UC = 1.5))

  ## --- density grids over the GL-touched quadratures -----------------------
  grid <- expand.grid(t = c(0.35, 0.6, 1.0, 1.8), mu = c(1, 3), A = c(0, 0.4),
                      sv = c(0, 0.6), lg = c(0, 0.8), lk = c(0, 0.5),
                      shape = c(1L, 2L), guess = c(FALSE, TRUE),
                      nodes = c(20L, 48L), posdrift = c(TRUE, FALSE))
  dens <- mapply(function(t, mu, A, sv, lg, lk, shape, guess, nodes, posdrift)
    EMC2:::drdmswtn(t, mu, 1.1, A, 1, 0.2, sv, lg, lk, nodes, FALSE, shape, guess, posdrift, 1.0),
    grid$t, grid$mu, grid$A, grid$sv, grid$lg, grid$lk, grid$shape, grid$guess,
    grid$nodes, grid$posdrift)
  cdfs <- mapply(function(t, mu, A, sv, lg, lk, shape, guess, nodes, posdrift)
    EMC2:::prdmswtn(t, mu, 1.1, A, 1, 0.2, sv, lg, lk, nodes, FALSE, shape, guess, posdrift, 1.0),
    grid$t, grid$mu, grid$A, grid$sv, grid$lg, grid$lk, grid$shape, grid$guess,
    grid$nodes, grid$posdrift)
  ggrid <- expand.grid(t = c(0.35, 0.6, 1.0, 1.8), mu = c(1, 3), A = c(0, 0.4),
                       lg = c(0, 0.8), lk = c(0, 0.5), shape = c(1L, 2L),
                       guess = c(FALSE, TRUE))
  gdens <- mapply(function(t, mu, A, lg, lk, shape, guess)
    EMC2:::dgbm(t, mu, 1.1, A, 1, 0.2, lg, lk, FALSE, shape, guess, 1.0),
    ggrid$t, ggrid$mu, ggrid$A, ggrid$lg, ggrid$lk, ggrid$shape, ggrid$guess)
  gcdf <- mapply(function(t, mu, A, lg, lk, shape, guess)
    EMC2:::pgbm(t, mu, 1.1, A, 1, 0.2, lg, lk, FALSE, shape, guess, 1.0),
    ggrid$t, ggrid$mu, ggrid$A, ggrid$lg, ggrid$lk, ggrid$shape, ggrid$guess)
  out$density_grids <- list(drdmswtn = dens, prdmswtn = cdfs, dgbm = gdens, pgbm = gcdf)
  cat(sprintf("%-22s drdmswtn=%d dgbm=%d values\n", "density grids", length(dens), length(gdens)))
  out
}

compare_all <- function(ref, new, tol = 1e-10) {
  stopifnot(identical(sort(names(ref)), sort(names(new))))
  worst <- 0; ok <- TRUE
  for (nm in names(ref)) {
    a <- if (nm == "density_grids") unlist(ref[[nm]]) else ref[[nm]]$lls
    b <- if (nm == "density_grids") unlist(new[[nm]]) else new[[nm]]$lls
    stopifnot(length(a) == length(b))
    same_nonfinite <- all((is.finite(a) == is.finite(b)) &
                          (!is.finite(a) | !is.finite(b) | TRUE)) &&
                      all(a[!is.finite(a)] == b[!is.finite(b)])
    d <- max(abs(a[is.finite(a) & is.finite(b)] - b[is.finite(a) & is.finite(b)]), 0)
    worst <- max(worst, d)
    status <- if (same_nonfinite && d <= tol) "OK" else {ok <- FALSE; "FAIL"}
    cat(sprintf("%-22s max|diff|=%.3e nonfinite-match=%s  %s\n", nm, d, same_nonfinite, status))
  }
  cat(sprintf("\nOverall: %s (worst %.3e, tol %.1e)\n", if (ok) "PASS" else "FAIL", worst, tol))
  invisible(ok)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  if (args[1] == "capture") {
    res <- capture_all()
    saveRDS(res, args[2])
    cat("saved", args[2], "\n")
  } else if (args[1] == "compare") {
    tol <- if (length(args) >= 4) as.numeric(args[4]) else 1e-10
    ok <- compare_all(readRDS(args[2]), readRDS(args[3]), tol)
    quit(status = if (ok) 0 else 1)
  }
}
