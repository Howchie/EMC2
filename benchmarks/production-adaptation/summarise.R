# Tables from the production-adaptation runs.
#   Rscript benchmarks/production-adaptation/summarise.R
suppressPackageStartupMessages(library(posterior))
out <- Sys.getenv("PADAPT_OUT", "benchmarks/production-adaptation/results")
files <- list.files(out, "^prod_.*\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
key <- function(r) paste(r$case, r$seed)
arms <- c("doubling", "frozen", "early", "extended")

# 0. Stopping under fit()'s sample rule: did the arm converge, after how many
# transitions, and at what CPU cost (production, and warm-up + production).
st <- do.call(rbind, lapply(res, function(r) if (!is.null(r$stop)) cbind(
  case = r$case, seed = r$seed, arm = r$arm, r$stop, prod_cpu = r$prod_cpu,
  total_cpu = r$warmup$burn_cpu + r$warmup$adapt_cpu + r$ext_cpu + r$prod_cpu)))
if (!is.null(st)) {
  cat("\n== Stopping per run\n")
  st$arm <- factor(st$arm, arms)
  print(st[order(st$case, st$seed, st$arm), ], digits = 3, row.names = FALSE)
  cat("\n== Stopping by arm: share converged, mean transitions and CPU\n")
  print(aggregate(cbind(converged, transitions, prod_cpu, total_cpu) ~ case + arm, st, mean),
        digits = 3, row.names = FALSE)
}

# 1. Movement by segment: the doubling arm against the frozen arm from the
# same warm-up. A rebuild that helps shows up as a ratio above 1 after it.
seg <- do.call(rbind, lapply(res, function(r) cbind(case = r$case, seed = r$seed, arm = r$arm, r$segments)))
w <- reshape(seg[seg$arm %in% c("doubling", "frozen"), c("case", "seed", "arm", "segment", "move_rate", "esjd_alpha", "esjd_mu")],
             idvar = c("case", "seed", "segment"), timevar = "arm", direction = "wide")
w$move_ratio <- w$move_rate.doubling / w$move_rate.frozen
w$esjd_alpha_ratio <- w$esjd_alpha.doubling / w$esjd_alpha.frozen
w$esjd_mu_ratio <- w$esjd_mu.doubling / w$esjd_mu.frozen
cat("\n== Segment movement, doubling / frozen (geometric mean over seeds)\n")
gm <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
segtab <- aggregate(cbind(move_ratio, esjd_alpha_ratio, esjd_mu_ratio) ~ case + segment, w, gm)
segtab$segment <- factor(segtab$segment, unique(seg$segment))
print(segtab[order(segtab$case, segtab$segment), ], digits = 3, row.names = FALSE)
cat("\n== Segment movement, frozen arm (absolute)\n")
fz <- aggregate(cbind(move_rate, move_rate_min, esjd_alpha, esjd_mu) ~ case + segment,
                seg[seg$arm == "frozen", ], mean)
fz$segment <- factor(fz$segment, unique(seg$segment))
print(fz[order(fz$case, fz$segment), ], digits = 3, row.names = FALSE)

# 2. Efficiency: worst-case bulk/tail ESS per CPU second. Production only,
# and total = warm-up + any extension + production.
eff <- do.call(rbind, lapply(res, function(r) {
  s <- r$summary
  total <- r$warmup$burn_cpu + r$warmup$adapt_cpu + r$ext_cpu + r$prod_cpu
  do.call(rbind, lapply(seq_len(nrow(s)), function(i) data.frame(
    case = r$case, seed = r$seed, arm = r$arm, what = s$what[i],
    ess_bulk_min = s$ess_bulk_min[i], ess_bulk_q05 = s$ess_bulk_q05[i],
    ess_tail_min = s$ess_tail_min[i], rhat_max = s$rhat_max[i],
    prod_cpu = r$prod_cpu, total_cpu = total,
    min_ess_per_ks_prod = 1000 * min(s$ess_bulk_min[i], s$ess_tail_min[i]) / r$prod_cpu,
    min_ess_per_ks_total = 1000 * min(s$ess_bulk_min[i], s$ess_tail_min[i]) / total)))
}))
cat("\n== Worst ESS (min of bulk, tail) per 1000 CPU s, ratio to doubling, geometric mean over seeds\n")
base <- eff[eff$arm == "doubling", c("case", "seed", "what", "min_ess_per_ks_prod", "min_ess_per_ks_total")]
names(base)[4:5] <- c("b_prod", "b_total")
m <- merge(eff, base, by = c("case", "seed", "what"))
m$ratio_prod <- m$min_ess_per_ks_prod / m$b_prod
m$ratio_total <- m$min_ess_per_ks_total / m$b_total
rt <- aggregate(cbind(ratio_prod, ratio_total) ~ case + what + arm, m, gm)
rt$arm <- factor(rt$arm, arms)
print(rt[order(rt$case, rt$what, rt$arm), ], digits = 3, row.names = FALSE)
cat("\n== Absolute: worst ESS, max R-hat, CPU (mean over seeds)\n")
ab <- aggregate(cbind(ess_bulk_min, ess_tail_min, rhat_max, prod_cpu, total_cpu) ~ case + what + arm, eff, mean)
ab$arm <- factor(ab$arm, arms)
print(ab[order(ab$case, ab$what, ab$arm), ], digits = 3, row.names = FALSE)

# 3. Do the arms agree? Posterior mean differences against the doubling arm,
# in Monte Carlo standard errors, for mu and log group SDs.
cat("\n== Agreement with doubling: share of |z| > 3 and max |z| (mu, logsd)\n")
pv <- do.call(rbind, lapply(res, function(r) cbind(case = r$case, seed = r$seed, arm = r$arm,
  r$per_var[r$per_var$what != "alpha", ])))
b <- pv[pv$arm == "doubling", c("case", "seed", "variable", "mean", "mcse_mean")]
names(b)[4:5] <- c("m0", "se0")
z <- merge(pv[pv$arm != "doubling", ], b, by = c("case", "seed", "variable"))
z$z <- (z$mean - z$m0) / sqrt(z$mcse_mean^2 + z$se0^2)
print(aggregate(z ~ case + arm, z, function(x) c(share_gt3 = mean(abs(x) > 3), max = max(abs(x)))),
      digits = 3, row.names = FALSE)
