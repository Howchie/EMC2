# Analytic two-mode conditional target for future cache validation.
# Posterior is an equally weighted mixture N(-3, .5^2), N(3, .5^2).
# Divide by the standard normal prior to define an unnormalised likelihood.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
log_target <- function(x) {
  a <- dnorm(x, -3, .5, log = TRUE); b <- dnorm(x, 3, .5, log = TRUE)
  m <- pmax(a, b); m + log((exp(a - m) + exp(b - m)) / 2)
}
ll_batch <- function(x) log_target(x[, 1]) - dnorm(x[, 1], log = TRUE)
objective <- function(x) -log_target(x)
starts <- c(-3.2, -2.8, 2.8, 3.2)
modes <- vapply(starts, function(x) optim(x, objective, method = 'BFGS')$par, 0)
stopifnot(all(abs(modes - c(-3, -3, 3, 3)) < 1e-4))
prior <- list(mu = 0, var = matrix(1))
set.seed(101)
checks <- lapply(modes, function(m) EMC2:::.emc_laplace_is_start(m, matrix(.25), prior, ll_batch))
result <- list(starts = starts, modes = modes, checks = checks,
               exact_mode_mass = c(left = .5, right = .5), final_seed = .Random.seed,
               session = sessionInfo())
saveRDS(result, Sys.getenv('OUT'))
print(result[c('starts', 'modes', 'exact_mode_mass')])
