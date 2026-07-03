set.seed(42)
devtools::load_all()
library(microbenchmark)

des <- design(factors = list(subjects=1, S=c("s1","s2")),
              Rlevels = c("s1","s2"),
              model = DDM,
              formula = list(v~1, a~1, t0~1, sv~1, SZ~1))

pars <- c(v = 1, a = log(1.5), t0 = log(0.2), sv = log(1), SZ = log(0.1))
data <- make_data(pars, des, n_trials=1000)

emc <- make_emc(data, des, type="single")
pars_matrix <- matrix(pars, nrow = 1000, ncol = length(pars), byrow=TRUE)
colnames(pars_matrix) <- names(pars)

invisible(EMC2:::calc_ll_manager(pars_matrix, emc[[1]]$data[[1]], emc[[1]]$model))

res <- microbenchmark(
  ll = EMC2:::calc_ll_manager(pars_matrix, emc[[1]]$data[[1]], emc[[1]]$model),
  times = 5
)
print(res)
