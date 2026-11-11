# Full SSEXG raw-likelihood benchmark. The default model's nonzero exg_lb and
# exgS_lb exercise both configurable lower bounds in the live trial route.
suppressPackageStartupMessages(
  library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
lI_go <- function(d) factor(rep(2L, nrow(d)), levels = 1:2)
des <- design(
  model = SSEXG(stop_method = "gl", stop_n_nodes = 64L),
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lI_go), covariates = "SSD",
  formula = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                 muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
p <- sampled_pars(des, doMap = FALSE)
p[c("mu_lMFALSE", "mu_lMTRUE", "sigma", "tau", "muS", "sigmaS", "tauS",
    "gf", "tf")] <- c(log(.65), log(.50), log(.05), log(.15), log(.20),
                        log(.04), log(.08), qnorm(.08), qnorm(.08))
dat <- make_data(p, des, n_trials = as.integer(Sys.getenv("N_TRIALS", "800")),
                 functions = list(SSD = make_ssd()))
emc <- make_emc(dat, des, type = "single", verbose = FALSE)
dadm <- emc[[1]]$data[[1]]
model <- emc[[1]]$model()
p_types <- names(model$p_types)
designs <- lapply(p_types, function(nm) {
  d <- attr(dadm, "designs")[[nm]]
  d[attr(d, "expand"), , drop = FALSE]
})
names(designs) <- p_types
constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
n_particles <- as.integer(Sys.getenv("N_PARTICLES", "32"))
set.seed(456)
particles <- matrix(rep(as.numeric(p), each = n_particles), nrow = n_particles) +
  matrix(rnorm(n_particles * length(p), sd = .15), nrow = n_particles)
colnames(particles) <- names(p)
fun <- get("calc_ll_oo", asNamespace("EMC2"))
call <- function() fun(particles, dadm, constants, designs, model$c_name,
                        model$bound, model$transform, model$pre_transform,
                        p_types, log(1e-10), model$trend)
result <- call()
n_rep <- as.integer(Sys.getenv("N_REP", "7"))
timing <- replicate(n_rep, system.time(invisible(call()))[["elapsed"]])
out <- Sys.getenv("OUT")
if (nzchar(out)) {
  data_file <- tempfile(); saveRDS(dadm, data_file)
  data_md5 <- unname(tools::md5sum(data_file)); unlink(data_file)
  saveRDS(list(result = result, timing = timing, data_md5 = data_md5,
               rows = nrow(dadm), particles = n_particles,
               session = sessionInfo()), out)
}
cat("unique rows", nrow(dadm), "particles", n_particles,
    "times", timing, "\n")
