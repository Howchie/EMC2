# RDM fixtures with RT-like subjects, regenerated from an existing fixture's
# design, group means and intended covariance.
#   CASE=corr15|freeA15|wide15 Rscript benchmarks/mode-start/plausible_fixtures.R
# Reads fixture_<case>.rds from MODE_START_SRC and writes the new one to
# MODE_START_OUT. A subject is kept only if its parameters and data are
# RT-like; otherwise it is redrawn from the same group distribution:
#   t0 in [0.15, 0.40] s, every cell's B in [0.5, 3.5], every accumulator's v
#   in [0.3, 6], every cell's median RT in [0.35, 1.5] s, and every RT, after
#   flooring to 1/60 s as make_emc() does, at least two frames above t0.
# The rejection truncates the group distribution slightly; the realised
# subject effects are saved with the fixture.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
src <- Sys.getenv("MODE_START_SRC", "benchmarks/mode-start/results")
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results_plausible")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
case <- Sys.getenv("CASE")
fx <- readRDS(file.path(src, paste0("fixture_", case, ".rds")))
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR
n_subj <- nlevels(fx$data$subjects)
n_trials <- as.integer(nrow(fx$data) / (n_subj * nlevels(fx$data$E) * nlevels(fx$data$S)))
des <- q(design(factors = list(subjects = seq_len(n_subj), E = levels(fx$data$E), S = levels(fx$data$S)),
                Rlevels = levels(fx$data$R), matchfun = matchfun, model = RDM,
                formula = fx$settings$formula, constants = fx$settings$constants))
# One-subject design for simulating a candidate subject.
des1 <- q(design(factors = list(subjects = 1L, E = levels(fx$data$E), S = levels(fx$data$S)),
                 Rlevels = levels(fx$data$R), matchfun = matchfun, model = RDM,
                 formula = fx$settings$formula, constants = fx$settings$constants))
gm <- fx$group_means
cv <- fx$intended_covariance

frame <- 1 / 60
plausible_pars <- function(p_vec) {
  m <- q(mapped_pars(des, p_vec))
  all(m$t0 >= 0.15 & m$t0 <= 0.40) && all(m$B >= 0.5 & m$B <= 3.5) &&
    all(m$v >= 0.3 & m$v <= 6)
}
plausible_data <- function(d, p_vec) {
  m <- q(mapped_pars(des, p_vec))
  t0 <- m$t0[1]
  med <- tapply(d$rt, list(d$E, d$S), median)
  all(med >= 0.35 & med <= 1.5) && min(floor(d$rt / frame) * frame) - t0 >= 2 * frame
}

set.seed(as.integer(Sys.getenv("FIXTURE_SEED", "20261007")))
se <- matrix(NA_real_, n_subj, length(gm), dimnames = list(NULL, names(gm)))
dat <- vector("list", n_subj)
redraws <- integer(n_subj)
for (s in seq_len(n_subj)) {
  repeat {
    p_vec <- q(make_random_effects(des1, gm, n_subj = 1, covariances = cv))[1, ]
    redraws[s] <- redraws[s] + 1L
    if (redraws[s] > 1000L) stop("subject ", s, ": no plausible draw in 1000 tries")
    if (!plausible_pars(p_vec)) next
    d <- q(make_data(matrix(p_vec, 1, dimnames = list("1", names(gm))), des1,
                     n_trials = n_trials, verbose = FALSE))
    d$subjects <- s
    if (!plausible_data(d, p_vec)) next
    se[s, ] <- p_vec
    dat[[s]] <- d
    break
  }
}
data <- do.call(rbind, dat)
data$subjects <- factor(data$subjects, levels = levels(fx$data$subjects))
rownames(data) <- NULL
cat(sprintf("%s: %d subjects, %d trials per cell | draws per accepted subject: median %d, max %d | accuracy %.3f, median RT %.3f\n",
            case, n_subj, n_trials, as.integer(median(redraws)), max(redraws),
            mean(data$S == data$R), median(data$rt)))
saveRDS(list(data = data, group_means = gm, intended_covariance = cv, subject_effects = se,
             settings = fx$settings, redraws = redraws,
             plausibility = "t0 [0.15, 0.40]; B [0.5, 3.5]; v [0.3, 6]; cell median RT [0.35, 1.5]; floored RT - t0 >= 2/60"),
        file.path(out, paste0("fixture_", case, ".rds")))
