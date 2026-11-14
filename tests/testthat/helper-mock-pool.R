# A pool of mock workers on real pipes, for tests of the chain's side of an
# exchange (src/pool_exchange.cpp). One forked process stands in for all the
# workers: it takes each worker's request in turn and `serve(w, msg)` says
# what that worker does, as a list with
#   records: the workers to announce on the completion channel (default `w`;
#            integer(0) for none),
#   reply:   the worker's reply; NULL sends none,
#   garbled: TRUE answers with bytes that are no message,
#   note:    anything the test wants back.
# `finish()` returns, per worker, the request it was sent and its note, and
# closes the pipes. NULL where FIFOs or forking are unavailable.
mock_worker_pool <- function(n_workers, serve, generation = 5L) {
  if (.Platform$OS.type == "windows" || !nzchar(Sys.which("mkfifo"))) return(NULL)
  dir <- tempfile("emc_mock_pool_")
  dir.create(dir)
  done <- EMC2:::.emc_wpool_done_open(dir)
  if (is.null(done)) return(NULL)
  pipe_pair <- function(name) {
    path <- file.path(dir, name)
    if (!identical(system2("mkfifo", shQuote(path), stdout = FALSE, stderr = FALSE), 0L))
      return(NULL)
    list(r = fifo(path, "rb", blocking = FALSE), w = fifo(path, "wb", blocking = FALSE))
  }
  req <- lapply(seq_len(n_workers), function(w) pipe_pair(sprintf("req%d", w)))
  ans <- lapply(seq_len(n_workers), function(w) pipe_pair(sprintf("ans%d", w)))
  if (any(vapply(c(req, ans), is.null, logical(1)))) return(NULL)
  server <- parallel::mcparallel({
    lapply(seq_len(n_workers), function(w) {
      msg <- EMC2:::.emc_wpool_recv(req[[w]]$r, deadline = 30)
      if (is.null(msg)) return(NULL)
      act <- serve(w, msg)
      records <- if (is.null(act$records)) w else act$records
      for (r in records) {
        writeBin(EMC2:::.emc_wpool_done_encode(r, generation, msg$request,
                                               EMC2:::.EMC_WPOOL_ST_OK), done$wc)
        flush(done$wc)
      }
      if (isTRUE(act$garbled)) {
        EMC2:::emc_wpool_send_cpp(ans[[w]]$w, as.raw(1:8), 30, 2^30, 4096)
      } else if (!is.null(act$reply)) {
        EMC2:::.emc_wpool_send(ans[[w]]$w, act$reply, deadline = 30)
      }
      list(msg = msg, note = act$note)
    })
  }, silent = TRUE)
  pool <- list(n = n_workers, dir = dir, jobs = vector("list", n_workers),
               wcs = lapply(req, `[[`, "w"), rcs = lapply(ans, `[[`, "r"),
               alive = TRUE, done = done, generation = generation)
  finish <- function() {
    log <- parallel::mccollect(server)[[1L]]
    for (cn in c(done[c("rc", "wc")], unlist(lapply(c(req, ans), function(p) p[c("r", "w")]),
                                              recursive = FALSE)))
      try(close(cn), silent = TRUE)
    unlink(dir, recursive = TRUE)
    log
  }
  list(pool = pool, finish = finish)
}
