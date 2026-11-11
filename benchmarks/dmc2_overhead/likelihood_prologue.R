# Measure the native parameter prologue against the complete likelihood call
# on the same EXG and LNR data used by paper_profile.R. This is a decision gate
# for a persistent native likelihood context, not a posterior benchmark.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))

case <- Sys.getenv("CASE", "EXG")
particle_grid <- sort(unique(as.integer(strsplit(
  Sys.getenv("NP_GRID", Sys.getenv("NP", "1,8,32,128,256")), ",",
  fixed = TRUE)[[1L]])))
n_batches <- as.integer(Sys.getenv("BATCHES", "5"))
full_reps <- as.integer(Sys.getenv("FULL_REPS", "5"))
prologue_reps <- as.integer(Sys.getenv("PROLOGUE_REPS", "100"))
stopifnot(case %in% c("EXG", "LNR"), length(particle_grid) > 1L,
          all(is.finite(particle_grid)), all(particle_grid > 0L),
          n_batches >= 3L, full_reps > 0L, prologue_reps > 0L)
n_particles <- max(particle_grid)

set.seed(123)
if (case == "EXG") {
  des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1,
                model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  centre <- log(c(mu = .4, sigma = .1, tau = .2))
  dat <- make_data(centre, des, n_trials = 1e4, rt_resolution = .001)
  dat <- make_missing(dat,
    UC = gamlss.dist::qexGAUS(.9, mu = .4, sigma = .1, nu = .2),
    rt_resolution = .001, verbose = FALSE)
} else {
  Dfun <- function(dadm) factor(abs(as.numeric(dadm$S) - as.numeric(dadm$lR)))
  des <- design(factors = list(subjects = 1, S = 1:3), Rlevels = 1:3,
                functions = list(D = Dfun),
                formula = list(m ~ D, s ~ 1, t0 ~ 1), model = LNR)
  centre <- sampled_pars(des)
  centre[] <- c(-.5, .5, 1, log(.75), log(.2))
  dat <- make_data(centre, des, n_trials = 1e4)
  dat <- make_missing(dat, LT = .41,
    UC = function(d) quantile(d$rt[is.finite(d$rt)], .9),
    rt_resolution = .001, verbose = FALSE)
}

emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                 rt_resolution = .001))
dadm <- emc[[1L]]$data[[1L]]
model_fn <- emc[[1L]]$model
model <- model_fn()
p_types <- names(model$p_types)
designs <- EMC2:::.oo_expanded_designs(dadm, expand = FALSE)
constants <- attr(dadm, "constants")
if (is.null(constants)) constants <- NA_real_

set.seed(101)
particles_all <- matrix(rnorm(n_particles * length(centre), sd = .01),
                        nrow = n_particles,
                        dimnames = list(NULL, names(centre)))
particles_all <- sweep(particles_all, 2L, as.numeric(centre), "+")
particles <- particles_all

native <- function() EMC2:::calc_ll_oo(
  particles, dadm, constants, designs, model$c_name, model$bound,
  model$transform, model$pre_transform, p_types, log(1e-10), model$trend)
manager <- function() EMC2:::calc_ll_manager(particles, dadm, model_fn, r_cores = 1L)
manager_validated <- function() EMC2:::calc_ll_manager(
  particles, dadm, model_fn, r_cores = 1L, ll_cache_validated = TRUE)
prologue <- function() EMC2:::pt_prologue_oo(
  particles, dadm, constants, designs, model$bound, model$transform,
  model$pre_transform, p_types, model$trend, return_pars = FALSE)

# Warm each path before recording batches. Batched timings keep sub-millisecond
# prologue calls above the system clock's resolution.
ll_native <- native()
ll_manager <- manager()
mapped <- prologue()
stopifnot(isTRUE(all.equal(ll_native, ll_manager, tolerance = 0)),
          length(mapped$n_ok) == n_particles)

measure <- function(fun, reps) {
  invisible(fun())
  vapply(seq_len(n_batches), function(i) {
    t0 <- proc.time()[["elapsed"]]
    for (j in seq_len(reps)) invisible(fun())
    1000 * (proc.time()[["elapsed"]] - t0) / reps
  }, numeric(1))
}

summarise <- function(x) c(median = median(x), q25 = unname(quantile(x, .25)),
                            q75 = unname(quantile(x, .75)))
manager_public_ms <- measure(manager, full_reps)
manager_block_ms <- measure(manager_validated, full_reps)
timings <- lapply(particle_grid, function(n) {
  particles <<- particles_all[seq_len(n), , drop = FALSE]
  list(n_particles = n,
       native_ms = measure(native, full_reps),
       prologue_ms = measure(prologue, prologue_reps))
})
names(timings) <- as.character(particle_grid)
summary <- lapply(timings, function(x) {
  c(n_particles = x$n_particles,
    native_ms = summarise(x$native_ms),
    prologue_ms = summarise(x$prologue_ms),
    prologue_share_pct = 100 * median(x$prologue_ms) / median(x$native_ms),
    native_remainder_ms = median(x$native_ms) - median(x$prologue_ms))
})
scaling <- data.frame(
  n_particles = particle_grid,
  native_ms = vapply(summary, function(x) x[["native_ms.median"]], numeric(1)),
  prologue_ms = vapply(summary, function(x) x[["prologue_ms.median"]], numeric(1)))
prologue_fit <- lm(prologue_ms ~ n_particles, data = scaling)

data_file <- tempfile(fileext = ".rds")
saveRDS(dadm, data_file)
data_md5 <- unname(tools::md5sum(data_file))
unlink(data_file)
out <- list(case = case, rows = nrow(dadm), particle_grid = particle_grid,
            data_md5 = data_md5, manager_public_ms = manager_public_ms,
            manager_block_ms = manager_block_ms,
            timings = timings, summary = summary, scaling = scaling,
            prologue_intercept_ms = unname(coef(prologue_fit)[[1L]]),
            prologue_slope_ms_per_particle = unname(coef(prologue_fit)[[2L]]),
            identical_manager_native = identical(ll_native, ll_manager),
            min_valid_rows = min(mapped$n_ok), max_valid_rows = max(mapped$n_ok),
            session = sessionInfo())
saveRDS(out, Sys.getenv("OUT"))

cat(sprintf("%s rows=%d data_md5=%s\n", case, nrow(dadm), data_md5))
cat(sprintf(paste0("manager public %.3f ms; block-validated manager %.3f ms; ",
                   "native bit-identical: %s\n"),
            median(manager_public_ms), median(manager_block_ms),
            identical(ll_native, ll_manager)))
for (i in seq_along(summary)) {
  x <- summary[[i]]
  cat(sprintf("N=%3d native %.3f ms; prologue %.3f ms (%.2f%%); remainder %.3f ms\n",
              x[["n_particles"]], x[["native_ms.median"]],
              x[["prologue_ms.median"]], x[["prologue_share_pct"]],
              x[["native_remainder_ms"]]))
}
cat(sprintf("linear prologue intercept %.3f ms, slope %.5f ms/particle\n",
            out$prologue_intercept_ms, out$prologue_slope_ms_per_particle))
