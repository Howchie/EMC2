suppressPackageStartupMessages(library(EMC2))

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designRDM <- design(
  factors   = list(subjects = 1, S = c("left", "right")),
  Rlevels   = c("left", "right"),
  matchfun  = matchfun,
  model     = RDM,
  constants = c(s = log(1)),
  formula   = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1)
)

designRDMSWTN <- design(
  factors   = list(subjects = 1, S = c("left", "right")),
  Rlevels   = c("left", "right"),
  matchfun  = matchfun,
  model     = RDMSWTN,
  constants = c(s = log(1), sv = log(0), lambda = log(0)),
  formula   = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1)
)

p_vec <- sampled_pars(designRDM, doMap = FALSE)
p_vec[] <- c(log(1.5), log(0.5), log(1), log(0.2), log(0.4))
dat <- make_data(p_vec, designRDM, n_trials = 2000)

emcRDM  <- make_emc(dat, designRDM,     type = "single", compress = TRUE)
emcSWTN <- make_emc(dat, designRDMSWTN, type = "single", compress = TRUE)

# ---- extract LL ingredients ----
make_ll_args <- function(emc_obj) {
  s     <- emc_obj[[1]]
  dadm  <- EMC2:::.cache_ll_data_attrs(s$data[[1]])
  model <- s$model()                             # call the closure
  p_types <- names(model$p_types)
  designs <- list()
  for (p in p_types)
    designs[[p]] <- attr(dadm,"designs")[[p]][
      attr(attr(dadm,"designs")[[p]],"expand"), , drop=FALSE]
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  list(dadm=dadm, model=model, p_types=p_types,
       designs=designs, constants=constants)
}

args_rdm  <- make_ll_args(emcRDM)
args_swtn <- make_ll_args(emcSWTN)

N <- 500
pars_rdm  <- sampled_pars(designRDM,     doMap = FALSE)
pars_swtn <- sampled_pars(designRDMSWTN, doMap = FALSE)
pars_rdm[]  <- c(log(1.5), log(0.5), log(1), log(0.2), log(0.4))
pars_swtn[] <- c(log(1.5), log(0.5), log(1), log(0.2), log(0.4))

pm_rdm  <- matrix(rep(pars_rdm,  N), nrow=N, byrow=TRUE, dimnames=list(NULL,names(pars_rdm)))
pm_swtn <- matrix(rep(pars_swtn, N), nrow=N, byrow=TRUE, dimnames=list(NULL,names(pars_swtn)))

call_ll <- function(pm, a)
  EMC2:::calc_ll_oo(pm, a$dadm, constants=a$constants, designs=a$designs,
             type=a$model$c_name, a$model$bound, a$model$transform,
             a$model$pre_transform, p_types=a$p_types,
             min_ll=log(1e-10), NULL)

# warm up
for (i in 1:5) { call_ll(pm_rdm, args_rdm); call_ll(pm_swtn, args_swtn) }

R <- 100
t_rdm  <- system.time(for (i in seq_len(R)) call_ll(pm_rdm,  args_rdm))["elapsed"]
t_swtn <- system.time(for (i in seq_len(R)) call_ll(pm_swtn, args_swtn))["elapsed"]

cat(sprintf("\nTrials: %d   Particles: %d   Reps: %d\n", nrow(dat), N, R))
cat(sprintf("  RDM                    : %6.3f s  (%5.2f ms/call)\n",
            t_rdm,  1000*t_rdm/R))
cat(sprintf("  RDMSWTN (sv=0,lam=0)  : %6.3f s  (%5.2f ms/call)\n",
            t_swtn, 1000*t_swtn/R))
cat(sprintf("  Ratio                  : %.2fx  (RDMSWTN / RDM)\n\n",
            t_swtn / t_rdm))
