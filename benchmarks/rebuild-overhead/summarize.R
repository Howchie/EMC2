# Export timing and memory summaries without copying any fixture data.
args <- commandArgs(TRUE)
root <- args[1L]; outdir <- args[2L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
profiles <- list()
for (file in list.files(file.path(root, "profile-nine"), pattern = "^profile-.+-(base|deliver)[.]rds$", full.names = TRUE)) {
  x <- readRDS(file)
  build <- if (grepl("-base[.]rds$", file)) "base" else "new"
  fixture <- sub("^profile-(.*)-(base|deliver)[.]rds$", "\\1", basename(file))
  p <- aggregate(x$profile[c("wall", "master_cpu", "worker_cpu", "idle_core_seconds")],
    list(component = x$profile$component), sum)
  p$fixture <- fixture; p$build <- build
  p$busy_cores <- with(p, ifelse(wall > 0, (master_cpu + worker_cpu) / wall, NA_real_))
  p$workers_used <- vapply(p$component, function(component) {
    j <- x$jobs$workers_used[x$jobs$component == component]
    if (length(j)) max(j) else 0L
  }, integer(1))
  profiles[[length(profiles) + 1L]] <- p
}
if (length(profiles)) {
  p <- do.call(rbind, profiles)
  write.csv(p, file.path(outdir, "components.csv"), row.names = FALSE)
  print(p, row.names = FALSE)
}
file <- file.path(root, "nine-core.rds")
if (file.exists(file)) {
  screen <- readRDS(file); x <- screen$timings
  # The first screen used an unwrapped auto_mclapply in its source overlay.
  # Its paired wall and identity results are valid; use installed profiles
  # for CPU attribution and idle core-seconds.
  columns <- if (isTRUE(screen$telemetry_complete))
    c("wall", "master_cpu", "worker_cpu", "idle_core_seconds") else "wall"
  p <- aggregate(x[columns], x[c("fixture", "component", "stage", "build")], mean)
  write.csv(p, file.path(outdir, "paired-components.csv"), row.names = FALSE)
  print(p, row.names = FALSE)
}
file <- file.path(root, "diagnostics.rds")
if (file.exists(file)) write.csv(readRDS(file), file.path(outdir, "diagnostics.csv"), row.names = FALSE)
file <- file.path(root, "memory-11-blocks.rds")
if (file.exists(file)) {
  x <- readRDS(file)
  write.csv(x$blocks, file.path(outdir, "memory-block-timings.csv"), row.names = FALSE)
  write.csv(x$memory, file.path(outdir, "memory-block-rss.csv"), row.names = FALSE)
  persistent <- split(x$memory, x$memory$pid)
  rows <- lapply(persistent, function(p) {
    p <- p[order(p$block), ]
    p <- p[is.finite(p$rss), ]
    if (!nrow(p)) return(NULL)
    data.frame(pid = p$pid[1L], blocks = nrow(p), first = p$rss[1L],
      last = tail(p$rss, 1L), peak = max(p$rss, na.rm = TRUE),
      range = diff(range(p$rss, na.rm = TRUE)))
  })
  write.csv(do.call(rbind, rows), file.path(outdir, "memory-summary.csv"), row.names = FALSE)
}

# Installed stopping checks and repeated dispatch lifetime.
rows <- list()
for (build in c("base", "deliver8")) {
  file <- file.path(root, paste0("diagnostics-", build, ".rds"))
  if (!file.exists(file)) next
  x <- readRDS(file)$timings
  x$build <- if (build == "base") "base" else "new"
  rows[[build]] <- x
}
if (length(rows)) write.csv(do.call(rbind, rows),
  file.path(outdir, "diagnostics-installed.csv"), row.names = FALSE)
file <- file.path(root, "diagnostics-memory.rds")
if (file.exists(file)) {
  x <- readRDS(file)
  write.csv(x$timings, file.path(outdir, "diagnostics-memory.csv"), row.names = FALSE)
  write.csv(x$dispatch, file.path(outdir, "diagnostics-worker-rss.csv"), row.names = FALSE)
}
if (length(profiles)) {
  p <- do.call(rbind, profiles)
  totals <- aggregate(p[c("wall", "master_cpu", "worker_cpu", "idle_core_seconds")],
    p[c("fixture", "build")], sum)
  write.csv(totals, file.path(outdir, "rebuild-totals.csv"), row.names = FALSE)
}
rows <- list()
for (rep in 1:2) for (build in c("base", "deliver8")) {
  file <- file.path(root, sprintf("forstmann-next-%d-%s.rds", rep, build))
  if (!file.exists(file)) next
  x <- readRDS(file)$profile
  rows[[length(rows) + 1L]] <- data.frame(fixture = "forstmann_sample_next_route",
    rep = rep, build = if (build == "base") "base" else "new", wall = sum(x$wall),
    master_cpu = sum(x$master_cpu), worker_cpu = sum(x$worker_cpu),
    idle_core_seconds = sum(x$idle_core_seconds))
}
if (length(rows)) write.csv(do.call(rbind, rows),
  file.path(outdir, "small-model-repeats.csv"), row.names = FALSE)
