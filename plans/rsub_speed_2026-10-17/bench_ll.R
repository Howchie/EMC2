# Full compiled likelihood cost per accumulator-row: LBA, RDM, RSUB dense,
# RSUB packets, timed on the first forstmann subject with 100 particles.
# Run from the repo root:  taskset -c <core> Rscript plans/rsub_speed_2026-10-17/bench_ll.R
# EMC_LIB selects the library holding the build to time (default .lib/rigs).
lib <- Sys.getenv("EMC_LIB", ".lib/rigs")
suppressMessages(library(EMC2, lib.loc = lib))
matchfun <- function(d) d$S == d$lR
dat <- forstmann; set.seed(1)
models <- list(LBA = LBA, RDM = RDM, RSUB = RSUB, RSUB_PK = function() RSUB("packets"))
n_particles <- 100; n_rep <- 15

setup <- function(m, nm, A0) {
  form <- list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1)
  cons <- if (nm == "LBA") c(sv = log(1)) else c(s = log(1))
  form <- form[vapply(form, function(f) as.character(f[[2]]) %in% names(m()$p_types), TRUE)]
  if (A0) cons <- c(cons, A = log(0))
  des <- suppressMessages(design(data = dat, model = m, matchfun = matchfun,
                                 formula = form, constants = cons))
  invisible(utils::capture.output(p <- sampled_pars(des, doMap = FALSE)))
  p[] <- 0
  p["v"] <- if (nm == "LBA") log(3) else log(2); p["v_lMTRUE"] <- .5
  p["B"] <- log(1); if (!A0) p["A"] <- log(.4); p["t0"] <- log(.15)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  dadm <- emc[[1]]$data[[1]]; mod <- emc[[1]]$model
  pm <- matrix(p, n_particles, length(p), byrow = TRUE,
               dimnames = list(NULL, names(p)))
  list(f = function() EMC2:::calc_ll_manager(pm, dadm, mod), n = nrow(dadm))
}

for (A0 in c(FALSE, TRUE)) {
  S <- lapply(names(models), function(nm) setup(models[[nm]], nm, A0))
  names(S) <- names(models)
  for (s in S) s$f()
  # Interleave models within each repetition so machine load hits all alike.
  tm <- matrix(NA, n_rep, length(S), dimnames = list(NULL, names(S)))
  for (r in seq_len(n_rep)) for (nm in names(S)) {
    t0 <- proc.time()[3]; S[[nm]]$f(); tm[r, nm] <- proc.time()[3] - t0
  }
  med <- apply(tm, 2, median) / (n_particles * S[[1]]$n) * 1e9
  cat(if (A0) "A = 0" else "A > 0", " ns per accumulator-row:",
      paste(names(med), round(med), collapse = "  "),
      "\n   relative to LBA:", paste(names(med), round(med / med["LBA"], 2), collapse = "  "), "\n")
}
