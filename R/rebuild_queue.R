# Rebuild jobs carry their own RNG stream or are deterministic. A shared
# counter balances uneven jobs without changing their result order. Forks
# inherit the inputs and exit after each rebuild, so model closures do not
# accumulate in the persistent sampling workers.
.emc_rebuild_queue <- new.env(parent = emptyenv())

.emc_rebuild_lapply <- function(X, FUN, n_cores, order = seq_along(X)) {
  n <- length(X)
  cores <- max(1L, min(n, as.integer(n_cores)))
  if (.Platform$OS.type != "unix" || cores <= 1L || n <= 2L * cores ||
      isTRUE(.emc_rebuild_queue$busy))
    return(auto_mclapply(X, FUN, mc.cores = cores, mc.set.seed = FALSE))
  stopifnot(length(order) == n, identical(sort(as.integer(order)), seq_len(n)))
  .emc_rebuild_queue$busy <- TRUE
  on.exit(.emc_rebuild_queue$busy <- FALSE, add = TRUE)
  map <- .emc_rebuild_queue$map
  if (is.null(map)) {
    path <- tempfile("emc-rebuild-claim-")
    on.exit(unlink(path), add = TRUE)
    map <- emc_claim_open(path)
    unlink(path)
    .emc_rebuild_queue$map <- map
  }
  emc_claim_reset(map, 1L)
  parts <- auto_mclapply(seq_len(cores), function(worker) {
    indices <- integer(n); values <- vector("list", n); used <- 0L
    repeat {
      claim <- emc_claim_next(map, 1L, n)
      if (claim < 0L) break
      i <- order[claim + 1L]
      used <- used + 1L
      indices[used] <- i
      values[used] <- list(FUN(X[[i]]))
    }
    list(indices = indices[seq_len(used)], values = values[seq_len(used)])
  }, mc.cores = cores, mc.set.seed = FALSE)
  out <- vector("list", n)
  for (part in parts) out[part$indices] <- part$values
  names(out) <- names(X)
  out
}
