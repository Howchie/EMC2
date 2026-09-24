# RNG isolation between test files.
#
# Many files call RNGkind("L'Ecuyer-CMRG") at top level, and until now none
# restored it.  Every later file in the same process then ran under L'Ecuyer,
# so a test that calls set.seed() without choosing a generator drew different
# numbers depending on which files ran before it.  With
# Config/testthat/parallel that also depends on how files are spread over the
# worker processes.
#
# local_rng_guard() records the generator kinds and .Random.seed and restores
# both when the calling frame exits.  Called at the top level of a test file,
# that is the end of the file, so whatever the file does to the generator
# stays inside it.  It changes nothing while the file runs.
local_rng_guard <- function(envir = parent.frame()) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  withr::defer({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, envir = envir)
  invisible(NULL)
}

# Pin the whole generator (kind, normal and sample kinds, and optionally the
# seed) for the calling frame, and restore it afterwards.  Use it in a test that
# draws random numbers, so the test does not depend on the generator it
# inherited.
local_test_rng <- function(seed = NULL, kind = "L'Ecuyer-CMRG",
                           envir = parent.frame()) {
  local_rng_guard(envir = envir)
  RNGkind(kind, "Inversion", "Rejection")
  if (!is.null(seed)) set.seed(seed)
  invisible(NULL)
}
