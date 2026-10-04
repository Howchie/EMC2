# Simulated fixtures with known truth for the burn benchmarks.
#   CASE=corr15|freeA15|wide15|lba71sim|jointsim Rscript benchmarks/mode-start/simulate_fixtures.R
# corr15/freeA15: the correlated RDM generator with group SDs x1.5 and 400
# trials per cell (same SD-to-subject-SE ratio as x3 at 100 trials).
# wide15: the 24-parameter wide RDM with group SD 0.3 and 600 trials per cell
# (same ratio as SD 0.5 at 200).
# lba71sim/jointsim: the lba71/joint designs (random 10/5-level factor F on
# the drift or location), simulated from fixed truths; F effects are
# subject-level noise around zero.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
case <- Sys.getenv("CASE")
file <- file.path(out, paste0("fixture_", case, ".rds"))
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR

if (case %in% c("corr15", "freeA15")) {
  env <- c(CORR_TRIALS_PER_CELL = "400", CORR_SD_SCALE = "1.5", CORR_STRENGTH = "1",
           CORR_FREE_A = if (case == "freeA15") "1" else "0",
           CORR_FIXTURE_OUT = file, EMC_LIB = Sys.getenv("EMC_LIB"))
  st <- system2("Rscript", "WorkingTests/generate_correlated_group_rdm_fixture.R",
                env = paste0(names(env), "=", env))
  if (st != 0) stop("fixture generation failed")
  quit(save = "no")
}

if (case == "wide15") {
  des <- q(design(factors = list(subjects = seq_len(30), E = c("speed", "neutral", "accuracy"),
                                 S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = matchfun, model = RDM,
    formula = list(v ~ lM * E, B ~ E * lR * S, t0 ~ E * S), constants = c(s = log(1), A = log(0))))
  gm <- sampled_pars(des); gm[] <- 0
  gm[c("v", "v_lMTRUE", "B", "t0")] <- c(log(2), log(1.5), log(1.1), log(0.25))
  cv <- diag(0.3^2, length(gm)); dimnames(cv) <- list(names(gm), names(gm))
  set.seed(20260929)
  se <- make_random_effects(des, gm, n_subj = 30, covariances = cv)
  dat <- make_data(se, des, n_trials = 600, verbose = FALSE)
  cat(case, ": accuracy", round(mean(dat$S == dat$R), 3), " median rt", round(median(dat$rt), 3), "\n")
  saveRDS(list(data = dat, group_means = gm, intended_covariance = cv, subject_effects = se,
               settings = list(formula = list(v ~ lM * E, B ~ E * lR * S, t0 ~ E * S),
                               constants = c(s = log(1), A = log(0)))), file)
  quit(save = "no")
}

n_subj <- 20L
n_trials <- 40L
fac <- function(n_F) list(subjects = seq_len(n_subj), S = c("left", "right"),
                          E = c("speed", "neutral", "accuracy"),
                          F = as.character(seq_len(n_F)))
# Group means and SDs for the non-F coefficients; F coefficients are 0 +- sd_F.
truth <- function(des, means, sds, sd_F) {
  pn <- names(sampled_pars(des))
  mu <- setNames(rep(0, length(pn)), pn)
  mu[names(means)] <- means
  sd <- setNames(ifelse(grepl("F", pn), sd_F, NA), pn)
  sd[names(sds)] <- sds
  stopifnot(!anyNA(sd))
  cov <- diag(sd^2, length(pn))
  dimnames(cov) <- list(pn, pn)
  list(mu = mu, cov = cov)
}
simulate <- function(des, tr, seed) {
  set.seed(seed)
  se <- make_random_effects(des, tr$mu, n_subj = n_subj, covariances = tr$cov)
  list(data = make_data(se, des, n_trials = n_trials, verbose = FALSE),
       subject_effects = se)
}
lba_design <- function(n_F) q(design(factors = fac(n_F), Rlevels = c("left", "right"),
  matchfun = matchfun, model = LBA, constants = c(sv = log(1)),
  formula = list(v ~ lM * E * F, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM)))
lba_truth <- function(des) truth(des,
  means = c(v = 1, v_lMTRUE = 1.5, v_Eneutral = 0.1, v_Eaccuracy = 0.2,
            B = log(1), B_Eneutral = 0.2, B_Eaccuracy = 0.4,
            A = log(0.5), t0 = log(0.2)),
  sds = c(v = 0.3, v_lMTRUE = 0.3, v_Eneutral = 0.2, v_Eaccuracy = 0.2,
          "v_lMTRUE:Eneutral" = 0.2, "v_lMTRUE:Eaccuracy" = 0.2,
          B = 0.2, B_Eneutral = 0.15, B_Eaccuracy = 0.15, B_lRright = 0.15,
          "B_Eneutral:lRright" = 0.1, "B_Eaccuracy:lRright" = 0.1,
          A = 0.2, t0 = 0.15, t0_Eneutral = 0.1, t0_Eaccuracy = 0.1, sv_lMTRUE = 0.1),
  sd_F = 0.3)

if (case == "lba71sim") {
  des <- lba_design(10L)
  tr <- lba_truth(des)
  sim <- simulate(des, tr, 71L)
  fx <- list(data = sim$data, subject_effects = sim$subject_effects, group_means = tr$mu,
             intended_covariance = tr$cov, model = "LBA")
} else if (case == "jointsim") {
  d1 <- lba_design(5L)
  d2 <- q(design(factors = fac(5L), Rlevels = c("left", "right"), matchfun = matchfun,
    model = LNR, formula = list(m ~ lM * E * F, s ~ lM, t0 ~ 1)))
  t1 <- lba_truth(d1)
  t2 <- truth(d2,
    means = c(m = -0.3, m_lMTRUE = -0.6, m_Eneutral = 0.1, m_Eaccuracy = 0.2,
              s = log(0.8), t0 = log(0.2)),
    sds = c(m = 0.15, m_lMTRUE = 0.15, m_Eneutral = 0.1, m_Eaccuracy = 0.1,
            "m_lMTRUE:Eneutral" = 0.1, "m_lMTRUE:Eaccuracy" = 0.1,
            s = 0.1, s_lMTRUE = 0.1, t0 = 0.15),
    sd_F = 0.15)
  s1 <- simulate(d1, t1, 72L)
  s2 <- simulate(d2, t2, 73L)
  fx <- list(data = list(s1$data, s2$data),
             subject_effects = list(s1$subject_effects, s2$subject_effects),
             group_means = list(t1$mu, t2$mu),
             intended_covariance = list(t1$cov, t2$cov), model = c("LBA", "LNR"))
} else stop("unknown CASE")

for (d in if (is.data.frame(fx$data)) list(fx$data) else fx$data) {
  err <- with(d, tapply(S != R, list(subjects, S, E, F), sum))
  cat(case, ": accuracy", round(mean(d$S == d$R), 3), " median rt", round(median(d$rt), 3),
      " cells with no errors", round(mean(err == 0), 2), "\n")
}
saveRDS(fx, file)
