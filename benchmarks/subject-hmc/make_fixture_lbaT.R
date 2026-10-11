# Practice-trend LBA fixture for the subject HMC panel: curved subject
# posteriors with a group-variance funnel.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/make_fixture_lbaT.R
# v ~ lM with an exponential practice trend v + w exp(-d_ed trial), B ~ 1,
# A ~ 1, t0 ~ 1, sv ~ lM (8 parameters). Few trials per subject leave each
# subject's (w, d_ed) pair on a banana-shaped ridge, so the group means are
# identified by pooling, and a small group SD of the log rate gives the
# funnel.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
n_subj <- as.integer(Sys.getenv("N_SUBJ", "80"))
n_trials <- as.integer(Sys.getenv("TRIALS", "150"))
seed <- as.integer(Sys.getenv("FIXTURE_SEED", "20261011"))
out <- Sys.getenv("OUT", "benchmarks/subject-hmc/fixtures/fixture_lbaT.rds")
q <- function(e) suppressMessages(e)
matchfun <- function(d) d$S == d$lR
form <- list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, sv ~ lM)
trend <- make_trend(par_names = "v", cov_names = "trial", kernels = "exp_decr", bases = "lin")
des <- q(design(factors = list(subjects = seq_len(n_subj), S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = matchfun, model = LBA,
                constants = c(sv = log(1)), formula = form, trend = trend,
                covariates = "trial"))
pn <- names(sampled_pars(des))
mu <- c(v = 1, v_lMTRUE = 1.5, B = log(1), A = log(0.5), t0 = log(0.2), sv_lMTRUE = 0,
        v.w = -1, v.d_ed = log(0.02))
sd <- c(v = 0.3, v_lMTRUE = 0.3, B = 0.2, A = 0.2, t0 = 0.15, sv_lMTRUE = 0.1,
        v.w = 0.3, v.d_ed = 0.15)
stopifnot(setequal(names(mu), pn), setequal(names(sd), pn))
mu <- mu[pn]; sd <- sd[pn]
cov <- diag(sd^2); dimnames(cov) <- list(pn, pn)

set.seed(seed)
se <- q(make_random_effects(des, mu, n_subj = n_subj, covariances = cov))
se <- matrix(as.numeric(se), n_subj, length(pn), dimnames = list(as.character(seq_len(n_subj)), pn))
skeleton <- do.call(rbind, lapply(seq_len(n_subj), function(s)
  data.frame(subjects = s, S = sample(c("left", "right"), n_trials, TRUE), trial = seq_len(n_trials),
             R = NA_character_, rt = NA_real_)))
skeleton$subjects <- factor(skeleton$subjects, levels = seq_len(n_subj))
skeleton$S <- factor(skeleton$S, levels = c("left", "right"))
skeleton$R <- factor(skeleton$R, levels = c("left", "right"))
data <- q(make_data(se, des, data = skeleton, verbose = FALSE))
rownames(data) <- NULL
early <- data$trial <= 30
cat(sprintf("lbaT: %d subjects x %d trials, %d parameters | accuracy early %.3f late %.3f | median RT early %.3f late %.3f\n",
            n_subj, n_trials, length(pn), mean((data$S == data$R)[early]), mean((data$S == data$R)[!early]),
            median(data$rt[early]), median(data$rt[!early])))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data, design = des, group_means = mu, intended_covariance = cov,
             subject_effects = se,
             settings = list(model = "LBA", formula = form, constants = c(sv = log(1)),
                             trend = trend, covariates = "trial", n_subj = n_subj,
                             n_trials = n_trials, seed = seed)), out)
