# Identified free-A LBA fixture for the subject HMC panel.
#   CASE=lbaA TRIALS=1000 A_MEAN=1.0 EMC_LIB=<lib> Rscript benchmarks/subject-hmc/make_fixture.R
# The lba71sim LBA design (benchmarks/mode-start/simulate_fixtures.R) without
# the random factor F: v ~ lM * E, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM (17
# parameters), the same group means and SDs except the start-point range
# (A_MEAN, natural scale), TRIALS per S x E cell.
# A subject is kept only if its likelihood is identified at its own
# multi-start MLE: observed information positive definite with condition
# number below MAX_COND at two finite-difference steps. Otherwise it is
# redrawn from the same group distribution; the redraw counts are saved.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
ns <- asNamespace("EMC2")
case <- Sys.getenv("CASE", "lbaA")
n_subj <- as.integer(Sys.getenv("N_SUBJ", "20"))
n_trials <- as.integer(Sys.getenv("TRIALS", "1000"))
A_mean <- as.numeric(Sys.getenv("A_MEAN", "1.0"))
max_cond <- as.numeric(Sys.getenv("MAX_COND", "1e6"))
n_starts <- as.integer(Sys.getenv("N_STARTS", "8"))
seed <- as.integer(Sys.getenv("FIXTURE_SEED", "20261007"))
out <- Sys.getenv("OUT", file.path("benchmarks/subject-hmc/fixtures", paste0("fixture_", case, ".rds")))
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR
form <- list(v ~ lM * E, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM)
make_des <- function(subjects) q(design(
  factors = list(subjects = subjects, S = c("left", "right"), E = c("speed", "neutral", "accuracy")),
  Rlevels = c("left", "right"), matchfun = matchfun, model = LBA,
  constants = c(sv = log(1)), formula = form))
des <- make_des(seq_len(n_subj))
des1 <- make_des(1L)
pn <- names(sampled_pars(des))
mu <- setNames(rep(0, length(pn)), pn)
mu[c("v", "v_lMTRUE", "v_Eneutral", "v_Eaccuracy", "B", "B_Eneutral", "B_Eaccuracy", "A", "t0")] <-
  c(1, 1.5, 0.1, 0.2, log(1), 0.2, 0.4, log(A_mean), log(0.2))
sd <- c(v = 0.3, v_lMTRUE = 0.3, v_Eneutral = 0.2, v_Eaccuracy = 0.2,
        "v_lMTRUE:Eneutral" = 0.2, "v_lMTRUE:Eaccuracy" = 0.2,
        B = 0.2, B_Eneutral = 0.15, B_Eaccuracy = 0.15, B_lRright = 0.15,
        "B_Eneutral:lRright" = 0.1, "B_Eaccuracy:lRright" = 0.1,
        A = 0.2, t0 = 0.15, t0_Eneutral = 0.1, t0_Eaccuracy = 0.1, sv_lMTRUE = 0.1)
stopifnot(setequal(names(sd), pn))
cov <- diag(sd[pn]^2); dimnames(cov) <- list(pn, pn)

identified <- function(d, truth) {
  emc1 <- q(make_emc(d, des1, n_chains = 1, type = "single"))
  ch <- emc1[[1]]
  f <- function(x) { v <- ns$.emc_burn_ll_batch(matrix(x, 1), ch$data[[1]], ch$model, pn)
                     if (is.finite(v)) v else -1e10 }
  starts <- c(list(truth), lapply(seq_len(n_starts - 1L), function(k) mu + rnorm(length(pn), 0, 1.5 * sd[pn])))
  fits <- lapply(starts, function(x0) optim(x0, function(x) -f(x), method = "BFGS",
                                             control = list(maxit = 2000, reltol = 1e-12)))
  best <- fits[[which.max(-sapply(fits, `[[`, "value"))]]$par
  cond <- sapply(list(NULL, rep(3e-4, length(pn))), function(step) {
    H <- -optimHess(best, f, control = if (is.null(step)) list() else list(ndeps = step))
    ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) Inf else max(ev) / min(ev)
  })
  max(cond)
}

set.seed(seed)
se <- matrix(NA_real_, n_subj, length(pn), dimnames = list(as.character(seq_len(n_subj)), pn))
dat <- vector("list", n_subj)
redraws <- integer(n_subj)
conds <- numeric(n_subj)
for (s in seq_len(n_subj)) {
  repeat {
    redraws[s] <- redraws[s] + 1L
    if (redraws[s] > 50L) stop("subject ", s, ": no identified draw in 50 tries")
    p_vec <- q(make_random_effects(des1, mu, n_subj = 1, covariances = cov))[1, ]
    d <- q(make_data(matrix(p_vec, 1, dimnames = list("1", pn)), des1, n_trials = n_trials, verbose = FALSE))
    cnd <- identified(d, p_vec)
    if (!is.finite(cnd) || cnd > max_cond) next
    d$subjects <- s
    se[s, ] <- p_vec
    dat[[s]] <- d
    conds[s] <- cnd
    break
  }
}
data <- do.call(rbind, dat)
data$subjects <- factor(data$subjects, levels = seq_len(n_subj))
rownames(data) <- NULL
cat(sprintf("%s: %d subjects x %d trials per cell, %d parameters | redraws per kept subject: %s | max condition %.3g | accuracy %.3f, median RT %.3f\n",
            case, n_subj, n_trials, length(pn), paste(redraws, collapse = ","), max(conds),
            mean(data$S == data$R), median(data$rt)))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data, design = des, group_means = mu, intended_covariance = cov,
             subject_effects = se, redraws = redraws, condition = conds,
             settings = list(formula = form, constants = c(sv = log(1)), n_subj = n_subj,
                             n_trials = n_trials, seed = seed, A_mean = A_mean,
                             max_cond = max_cond, n_starts = n_starts)), out)
