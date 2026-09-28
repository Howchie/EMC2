# Aggregate completed suites; per-run metrics remain available in CSV files.
root <- "benchmarks/group-move"
files <- file.path(root,"results",c("development","holdout"),"metrics.csv")
x <- do.call(rbind,lapply(files[file.exists(files)],read.csv))
runs <- lapply(split(x,interaction(x$suite,x$p,x$obs_max,x$method,x$seed,drop=TRUE)),function(z) {
  data.frame(suite=z$suite[1],p=z$p[1],obs_max=z$obs_max[1],method=z$method[1],
    seed=z$seed[1],worst_rhat=max(z$rhat),min_bulk=min(z$bulk_ess),min_tail=min(z$tail_ess),
    variance_min=min(z$variance),variance_max=max(z$variance),max_abs_mean=max(abs(z$mean)),
    seconds=z$elapsed[1],ess_per_second=min(z$bulk_ess)/z$elapsed[1])
})
runs <- do.call(rbind,runs)
write.csv(runs,file.path(root,"results","run-summary.csv"),row.names=FALSE)
lines <- c("# Analytic benchmark results", "", "Each row summarizes three independent replicate runs, each with four chains.",
           "R-hat is the worst direction/replicate; ESS is the median across replicates of the minimum direction ESS.",
           "Variance ranges cover all directions and replicates; the analytic value is 1 after standardization.",
           "Seconds include adaptation, settling and retained sampling for all four chains on a shared machine.", "",
           "| Suite | d | Observation max eigenvalue | Method | Replicates | Worst R-hat | Median min bulk ESS | Median min tail ESS | Variance range | Median seconds | Median min ESS/s |",
           "|---|---:|---:|---|---:|---:|---:|---:|---|---:|---:|")
for(z in split(runs,interaction(runs$suite,runs$p,runs$obs_max,runs$method,drop=TRUE))) {
  lines <- c(lines,sprintf("| %s | %d | %g | %s | %d | %.3f | %.0f | %.0f | %.2f–%.2f | %.1f | %.1f |",
    z$suite[1],z$p[1],z$obs_max[1],z$method[1],nrow(z),max(z$worst_rhat),median(z$min_bulk),
    median(z$min_tail),min(z$variance_min),max(z$variance_max),median(z$seconds),median(z$ess_per_second)))
}
writeLines(lines,file.path(root,"results","SUMMARY.md"))
cat(paste(lines,collapse="\n"),"\n")
