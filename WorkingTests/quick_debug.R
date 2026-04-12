
library(EMC2)
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
designLBA <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = LBA(posdrift=TRUE), 
  formula = list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM),
  constants = c(A = log(0.4))
)
p_vector <- sampled_pars(designLBA, doMap = FALSE)
p_vector[["B"]] <- log(1.2)
p_vector[["A"]] <- log(.4)
p_vector[["t0"]] <- log(0.15)
p_vector[["v_lMFALSE"]] <- .4
p_vector[["v_lMTRUE"]] <- 1.2
p_vector[["sv_lMFALSE"]] <- log(1.2)
p_vector[["sv_lMTRUE"]] <- log(1)

set.seed(123)
dat <- make_data(p_vector, designLBA, n_trials = 100, TC=list(UC = Inf))
emc <- make_emc(dat, designLBA, type = "single", rt_resolution = 1/60)
dadm <- emc[[1]]$data[[1]]
ll <- EMC2:::calc_ll_manager(t(as.matrix(p_vector)), dadm, designLBA$model)
cat("LL: ", ll, "\n")
