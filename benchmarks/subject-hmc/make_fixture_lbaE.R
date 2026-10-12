# Practice-trend LBA fixture with a conventional design: speed/accuracy
# emphasis and response bias on the threshold gap, a stimulus effect on the
# rates, A and t0 as intercepts.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/make_fixture_lbaE.R
# v ~ lM + S with an exponential practice trend v + w exp(-d_ed trial),
# B ~ E + lR, A ~ 1, t0 ~ 1, sv ~ lM (11 parameters). Emphasis alternates in
# blocks so both levels span the practice curve. Speed emphasis leaves a
# small gap b - A, so start-point noise costs accuracy (about 0.87 -> 0.82,
# median RT 0.66 -> 0.38 s); that condition, with A shared across emphasis,
# is what identifies A. Every subject's information at its MLE is positive
# definite (screen_identification.R).
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
n_subj <- as.integer(Sys.getenv("N_SUBJ", "80"))
n_trials <- as.integer(Sys.getenv("TRIALS", "240"))
block <- as.integer(Sys.getenv("BLOCK", "30"))
seed <- as.integer(Sys.getenv("FIXTURE_SEED", "20261012"))
out <- Sys.getenv("OUT", "benchmarks/subject-hmc/fixtures/fixture_lbaE.rds")
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR
form <- list(v ~ lM + S, B ~ E + lR, A ~ 1, t0 ~ 1, sv ~ lM)
trend <- make_trend(par_names = "v", cov_names = "trial", kernels = "exp_decr", bases = "lin")
des <- q(design(factors = list(subjects = seq_len(n_subj), S = c("left", "right"),
                               E = c("accuracy", "speed")),
                Rlevels = c("left", "right"), matchfun = matchfun, model = LBA,
                constants = c(sv = log(1)), formula = form, trend = trend,
                covariates = "trial"))
pn <- names(sampled_pars(des))
mu <- c(v = 1, v_lMTRUE = 2, v_Sright = -0.3, B = log(0.8), B_Espeed = -1.8, B_lRright = 0.1,
        A = log(0.7), t0 = log(0.2), sv_lMTRUE = 0, v.w = -1, v.d_ed = log(0.02))
sd <- c(v = 0.3, v_lMTRUE = 0.3, v_Sright = 0.15, B = 0.2, B_Espeed = 0.15, B_lRright = 0.1,
        A = 0.2, t0 = 0.15, sv_lMTRUE = 0.1, v.w = 0.3, v.d_ed = 0.15)
stopifnot(setequal(names(mu), pn), setequal(names(sd), pn))
mu <- mu[pn]; sd <- sd[pn]
cov <- diag(sd^2); dimnames(cov) <- list(pn, pn)

set.seed(seed)
se <- q(make_random_effects(des, mu, n_subj = n_subj, covariances = cov))
se <- matrix(as.numeric(se), n_subj, length(pn), dimnames = list(as.character(seq_len(n_subj)), pn))
emph <- c("accuracy", "speed")[((seq_len(n_trials) - 1) %/% block) %% 2 + 1]
skeleton <- do.call(rbind, lapply(seq_len(n_subj), function(s)
  data.frame(subjects = s, S = sample(c("left", "right"), n_trials, TRUE), E = emph,
             trial = seq_len(n_trials), R = NA_character_, rt = NA_real_)))
skeleton$subjects <- factor(skeleton$subjects, levels = seq_len(n_subj))
skeleton$S <- factor(skeleton$S, levels = c("left", "right"))
skeleton$E <- factor(skeleton$E, levels = c("accuracy", "speed"))
skeleton$R <- factor(skeleton$R, levels = c("left", "right"))
data <- q(make_data(se, des, data = skeleton, verbose = FALSE))
rownames(data) <- NULL
acc <- tapply(data$S == data$R, data$E, mean)
mrt <- tapply(data$rt, data$E, median)
cat(sprintf("lbaE: %d subjects x %d trials, %d parameters | accuracy %s %.3f, %s %.3f | median RT %s %.3f, %s %.3f | P(R=right) %.3f\n",
            n_subj, n_trials, length(pn), names(acc)[1], acc[1], names(acc)[2], acc[2],
            names(mrt)[1], mrt[1], names(mrt)[2], mrt[2], mean(data$R == "right")))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data, design = des, group_means = mu, intended_covariance = cov,
             subject_effects = se,
             settings = list(model = "LBA", formula = form, constants = c(sv = log(1)),
                             trend = trend, covariates = "trial", n_subj = n_subj,
                             n_trials = n_trials, block = block, seed = seed)), out)
