# One-step stationary probe of the existing auxiliary-centred conditional-IS
# kernel, with an exactly known Gaussian posterior. No adaptation or fitting.
set.seed(20260928L)
probe <- function(d, posterior_sd, epsilon, repetitions = 2000L) {
  n <- round(50 * sqrt(d) * .4)
  moved <- esjd <- numeric(repetitions)
  for (it in seq_len(repetitions)) {
    x <- rnorm(d, sd = posterior_sd)
    z <- x + rnorm(d, sd = .5 * epsilon)
    cloud <- rbind(x, sweep(matrix(rnorm(n * d, sd = .5 * epsilon), n, d), 2, z, '+'))
    # N(z | x, S) / N(x | z, S) cancels exactly.
    logw <- -.5 * rowSums(cloud^2) / posterior_sd^2
    k <- sample.int(n + 1L, 1L, prob = exp(logw - max(logw)))
    moved[it] <- k != 1L
    esjd[it] <- sum((cloud[k, ] - x)^2) / (d * posterior_sd^2)
  }
  data.frame(d, posterior_sd, epsilon, moved = mean(moved), esjd = mean(esjd))
}
out <- do.call(rbind, lapply(c(13L, 24L), function(d)
  do.call(rbind, lapply(c(.1, .05, .02), function(sd)
    do.call(rbind, lapply(c(.4, .1, .04, .01), function(eps) probe(d, sd, eps)))))))
print(out, row.names = FALSE)
write.csv(out, 'benchmarks/subject-move/results/scale-probe.csv', row.names = FALSE)
