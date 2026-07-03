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

# Generate extreme parameters for 100 particles
pars_matrix <- matrix(NA, nrow = 100, ncol = length(pars))
colnames(pars_matrix) <- names(pars)

for (i in 1:100) {
  pars_matrix[i, "v"] <- runif(1, -10, 10)
  pars_matrix[i, "a"] <- log(runif(1, 0.5, 5))
  pars_matrix[i, "t0"] <- log(runif(1, 0.05, 0.5))
  pars_matrix[i, "sv"] <- log(runif(1, 1e-4, 5))
  pars_matrix[i, "SZ"] <- log(runif(1, 0.1, 0.99)) # very high SZ
}

invisible(EMC2:::calc_ll_manager(pars_matrix, emc[[1]]$data[[1]], emc[[1]]$model))

res <- microbenchmark(
  ll = EMC2:::calc_ll_manager(pars_matrix, emc[[1]]$data[[1]], emc[[1]]$model),
  times = 5
)
print(res)
