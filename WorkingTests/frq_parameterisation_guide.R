# FRQ: the continuous relaxation, and where to start on the RT scale.
#
# Two things, both console-only (plus an optional PDF of the densities):
#
#   1. What the continuous relaxation actually is -- how a non-integer `alpha`
#      relates to the literal "N cues, wait for K" process, demonstrated
#      rather than asserted, and why both shapes are bounded below at 1.
#   2. A table of sensible starting parameterisations at human RT scale, with
#      the RT quantiles, omission rates and race accuracies they imply, so you
#      can pick a plausible neighbourhood before fitting anything.
#
# Run with:  Rscript WorkingTests/frq_parameterisation_guide.R

library(EMC2)

set.seed(20260816)

# ---------------------------------------------------------------------------
# Helpers.  Everything here is the *latent Beta* representation, which is what
# makes the relaxation work; the kernel's own inversion is EMC2:::frq_rate.
# ---------------------------------------------------------------------------

# Generative pair (p, lambda) implied by the estimation coordinates.
gen <- function(alpha, beta, h, tau) {
  pl <- EMC2:::frq_rate(alpha, beta, h, tau)
  c(p = unname(pl[1, "p"]), lambda = unname(pl[1, "lambda"]))
}

# Conditional quantile of decision time: F(x) = r * h.  Closed form, because
# q(x) = p(1 - exp(-lambda x)) is monotone with an explicit inverse.
qfrq <- function(r, alpha, beta, h, tau) {
  g <- gen(alpha, beta, h, tau)
  u <- qbeta(r * h, alpha, beta)
  -log1p(-u / g[["p"]]) / g[["lambda"]]
}

# Draw decision times from the latent-Beta representation directly (Inf when
# the quorum exceeds the attainable ceiling p).
rfrq_latent <- function(n, alpha, beta, h, tau) {
  g <- gen(alpha, beta, h, tau)
  u <- rbeta(n, alpha, beta)
  # Not ifelse(): that would evaluate log1p() on the u > p branch too and warn.
  out <- rep(Inf, n)
  fin <- u <= g[["p"]]
  out[fin] <- -log1p(-u[fin] / g[["p"]]) / g[["lambda"]]
  out
}

# The literal finite reservoir: N units, each available w.p. p, available ones
# registering at Exp(lambda); respond on the K-th registration.
rfrq_literal <- function(n, N, K, p, lambda) {
  avail <- matrix(runif(n * N) < p, nrow = n)
  lat <- matrix(rexp(n * N, lambda), nrow = n)
  lat[!avail] <- Inf
  apply(lat, 1, function(r) sort(r)[K])
}

rule <- function(txt) cat("\n", txt, "\n", strrep("-", nchar(txt)), "\n", sep = "")

# ===========================================================================
# PART 1.  Literal finite reservoir vs continuous relaxation
# ===========================================================================

rule("1a. Integer shapes ARE the literal finite-reservoir process")

# alpha = K, beta = N - K + 1.  Pick the generative pair first so the literal
# simulator and the closed form are being handed the same process.
N <- 7; K <- 3; p_lit <- 0.80; lam_lit <- 2.5
alpha <- K; beta <- N - K + 1
h_lit <- pbeta(p_lit, alpha, beta)
tau_lit <- -log1p(-qbeta(0.5 * h_lit, alpha, beta) / p_lit) / lam_lit

cat(sprintf("literal:  N = %d, K = %d, p = %.2f, lambda = %.2f\n",
            N, K, p_lit, lam_lit))
cat(sprintf("exposed:  alpha = %.0f, beta = %.0f, h = %.4f, tau = %.4f s\n",
            alpha, beta, h_lit, tau_lit))

Tk <- rfrq_literal(4e4, N, K, p_lit, lam_lit)
probe <- quantile(Tk[is.finite(Tk)], c(.1, .25, .5, .75, .9))
emp <- vapply(probe, function(x) mean(Tk <= x), numeric(1))
thy <- EMC2:::pfrq(as.numeric(probe), alpha, beta, h_lit, tau_lit)
cat(sprintf("max |F_empirical - F_closed_form|  = %.2e   (4e4 draws)\n",
            max(abs(emp - thy))))
cat(sprintf("omission rate: literal %.4f  vs  1 - h = %.4f\n",
            mean(is.infinite(Tk)), 1 - h_lit))
cat("=> at integer shapes the closed form is not an approximation to the\n",
    "   finite-cue process, it IS that process.\n", sep = "")

rule("1b. Non-integer shapes: no finite cue set, but an exact latent quorum")

# There is no (N, K) with K = 3.4, so `rfrq_literal` cannot be written at all.
# The latent representation U ~ Beta(alpha, beta) still holds exactly:
#   terminate when the cumulative registration state q(x) = p(1 - e^{-lambda x})
#   first reaches U;  never terminate if U > p.
alpha_c <- 3.4; beta_c <- 5.0; h_c <- 0.92; tau_c <- 0.30
Tc <- rfrq_latent(4e4, alpha_c, beta_c, h_c, tau_c)
probe <- quantile(Tc[is.finite(Tc)], c(.1, .25, .5, .75, .9))
emp <- vapply(probe, function(x) mean(Tc <= x), numeric(1))
thy <- EMC2:::pfrq(as.numeric(probe), alpha_c, beta_c, h_c, tau_c)
cat(sprintf("alpha = %.1f, beta = %.1f  ->  no integer (N, K) exists\n",
            alpha_c, beta_c))
cat(sprintf("max |F_latent_sim - F_closed_form| = %.2e   (4e4 draws)\n",
            max(abs(emp - thy))))
cat("=> the quorum U becomes a CONTINUOUS requirement on cumulative\n",
    "   registration probability rather than a count of cues.  The\n",
    "   distribution is exact; only the 'N discrete units' story is lost.\n",
    sep = "")

rule("1c. The relaxation interpolates the integer cases smoothly")

# Same h and tau throughout, so every row has the same completion probability
# and the same conditional median: alpha moves the SHAPE only.  Note that the
# reported N = alpha + beta - 1 is literal only at integer shapes.
cat(sprintf("%6s %6s %8s %8s %9s %9s %9s\n",
            "alpha", "beta", "N", "p", "q10", "q50", "q90"))
for (a in c(1, 1.5, 2, 2.5, 3, 4)) {
  b <- 3
  g <- gen(a, b, 0.95, 0.30)
  qs <- qfrq(c(.1, .5, .9), a, b, 0.95, 0.30)
  cat(sprintf("%6.1f %6.1f %8.1f %8.3f %9.3f %9.3f %9.3f\n",
              a, b, a + b - 1, g[["p"]], qs[1], qs[2], qs[3]))
}
cat("=> alpha is the leading-edge order: f(x) ~ x^(alpha-1) near zero, so\n",
    "   alpha = 1 leaves the density positive at the origin and alpha > 1\n",
    "   pushes it down.  Non-integer alpha sits genuinely between the\n",
    "   integer quorum depths, it is not a rounding of them.\n", sep = "")

# Confirm the leading-edge order numerically: the log-density slope at small x.
cat("\nleading-edge check, d log f / d log x as x -> 0 (should be alpha - 1):\n")
for (a in c(1.5, 2.5, 4)) {
  x <- c(1e-5, 2e-5)
  lf <- EMC2:::dfrq(x, a, 3, 0.95, 0.30, TRUE)
  cat(sprintf("  alpha = %.1f -> %.4f\n", a, diff(lf) / diff(log(x))))
}

rule("1d. Why the shapes are bounded below at 1")

cat("FRQ() alpha bounds: [",
    paste(FRQ()$bound$minmax[, "alpha"], collapse = ", "), "]\n")
cat("alpha < 1 is a mathematically valid transformed-Beta family, but the\n",
    "(h, tau) coordinates are not representable there in double precision,\n",
    "so it is not offered.  Two ways it fails:\n", sep = "")
cat(sprintf("  qbeta(0.99, 0.05, 0.05) = %.17g  -> p saturates at 1, so the\n",
            qbeta(0.99, 0.05, 0.05)))
cat("     fitted distribution is PROPER and h = 0.99 was silently ignored.\n")
cat(sprintf("  qbeta(0.10, 1e-4, 1e-4) = %.17g  -> p and u collapse together,\n",
            qbeta(0.1, 1e-4, 1e-4)))
cat("     so the kernel rejects an interior point: an artificial sampler cliff.\n")
cat("The dead zone reaches alpha = 0.5 at h = 1 - 1e-9, well inside anything\n",
    "a sampler would visit, which is why this is a bound and not a warning.\n",
    sep = "")

# ===========================================================================
# PART 2.  Sensible RT-scale starting parameterisations
# ===========================================================================

rule("2a. Single-accumulator starting points (seconds)")

# tau is the CONDITIONAL MEDIAN DECISION TIME, so the observed median RT of a
# lone accumulator is t0 + tau.  That is the whole trick to setting these: pick
# the RT median you want, subtract t0, and that is tau.
cfg <- list(
  "fast, sharp"        = c(alpha = 2.0, beta = 2.0, h = 0.99,  tau = 0.20, t0 = 0.15),
  "typical"            = c(alpha = 2.0, beta = 3.0, h = 0.98,  tau = 0.30, t0 = 0.20),
  "slow, skewed"       = c(alpha = 1.5, beta = 4.0, h = 0.97,  tau = 0.45, t0 = 0.22),
  "hard (low h)"       = c(alpha = 2.0, beta = 3.0, h = 0.80,  tau = 0.35, t0 = 0.20),
  "near non-defective" = c(alpha = 2.0, beta = 3.0, h = 0.999, tau = 0.30, t0 = 0.20)
)

cat(sprintf("%-20s %6s %6s %7s %7s %6s | %6s %7s | %6s %6s %6s %6s\n",
            "", "alpha", "beta", "h", "tau", "t0",
            "p", "lambda", "q10", "q50", "q90", "q99"))
for (nm in names(cfg)) {
  v <- cfg[[nm]]
  g <- gen(v[["alpha"]], v[["beta"]], v[["h"]], v[["tau"]])
  qs <- v[["t0"]] + qfrq(c(.1, .5, .9, .99), v[["alpha"]], v[["beta"]],
                         v[["h"]], v[["tau"]])
  cat(sprintf("%-20s %6.2f %6.2f %7.4f %7.2f %6.2f | %6.3f %7.2f | %6.3f %6.3f %6.3f %6.3f\n",
              nm, v[["alpha"]], v[["beta"]], v[["h"]], v[["tau"]], v[["t0"]],
              g[["p"]], g[["lambda"]], qs[1], qs[2], qs[3], qs[4]))
}
cat("\nRule of thumb: t0 + tau IS the conditional median RT, and t0 is the\n",
    "leading edge.  Set those two from the data first; alpha/beta then only\n",
    "trade off how heavy the right tail is against how sharp the rise is.\n",
    sep = "")

rule("2b. What h does at the RACE level")

# A two-accumulator race fails only if BOTH accumulators fail, so per-
# accumulator defect is far less visible than it looks.
cat(sprintf("%10s %10s %14s\n", "h_correct", "h_error", "race omission"))
for (hc in c(0.999, 0.99, 0.98, 0.95, 0.90, 0.80)) {
  he <- hc - 0.05
  cat(sprintf("%10.3f %10.3f %14.5f\n", hc, he, (1 - hc) * (1 - he)))
}
cat("=> if your data have ~0-1% omissions, h near 0.97-0.995 per accumulator\n",
    "   is already the right neighbourhood; do NOT read the omission rate as\n",
    "   1 - h directly.\n", sep = "")

rule("2c. Two-accumulator races: accuracy and RT together")

# Cross the design on h (attainability) and tau (speed given completion) --
# the recommended structure -- and read off what each buys you.
race <- function(n, cor, err) {
  Tc <- rfrq_latent(n, cor[["alpha"]], cor[["beta"]], cor[["h"]], cor[["tau"]])
  Te <- rfrq_latent(n, err[["alpha"]], err[["beta"]], err[["h"]], err[["tau"]])
  win <- ifelse(Tc < Te, "correct", "error")
  rt <- pmin(Tc, Te)
  om <- !is.finite(rt)
  ok <- !om
  list(omit = mean(om),
       acc = mean(win[ok] == "correct"),
       q = quantile(rt[ok], c(.1, .5, .9)),
       med_c = median(rt[ok][win[ok] == "correct"]),
       med_e = median(rt[ok][win[ok] == "error"]))
}

scen <- list(
  "tau only (timing)"    = list(cor = c(alpha=2, beta=3, h=0.98, tau=0.28),
                                err = c(alpha=2, beta=3, h=0.98, tau=0.55)),
  "h only (availability)"= list(cor = c(alpha=2, beta=3, h=0.98, tau=0.35),
                                err = c(alpha=2, beta=3, h=0.55, tau=0.35)),
  "dual (both)"          = list(cor = c(alpha=2, beta=3, h=0.98, tau=0.28),
                                err = c(alpha=2, beta=3, h=0.75, tau=0.45))
)
t0 <- 0.20
cat(sprintf("%-22s %7s %9s %7s %7s %7s | %8s %8s\n",
            "", "acc", "omission", "q10", "q50", "q90", "med(cor)", "med(err)"))
for (nm in names(scen)) {
  r <- race(4e4, scen[[nm]]$cor, scen[[nm]]$err)
  cat(sprintf("%-22s %7.3f %9.5f %7.3f %7.3f %7.3f | %8.3f %8.3f\n",
              nm, r$acc, r$omit, t0 + r$q[1], t0 + r$q[2], t0 + r$q[3],
              t0 + r$med_c, t0 + r$med_e))
}
cat("=> all three land near the same accuracy on completely different RT:\n",
    "   * tau-only gives SLOW errors (the classic race signature): an error\n",
    "     only wins when a genuinely slow correct draw loses to it.\n",
    "   * h-only gives FAST errors, the opposite sign.  Holding tau fixed\n",
    "     while lowering h lowers the attainable ceiling p and RAISES lambda\n",
    "     to keep the conditional median at tau, so the weak accumulator is\n",
    "     not slower -- it is quicker but usually absent.  It therefore only\n",
    "     ever wins early.  It is also slower overall, because a failed\n",
    "     competitor removes the faster of the two finishing times.\n",
    "   * dual cancels the two and gives no correct/error difference at all.\n",
    "   * the omission column separates them independently -- h-only leaks an\n",
    "     order of magnitude more no-responses at the same accuracy.\n",
    "   So h and tau are not two routes to the same fit: choose between them\n",
    "   on the error/correct RT ordering, then let the sampler tune.\n", sep = "")

rule("2d. The same thing as an EMC2 design, on the sampled (transformed) scale")

matchfun <- function(d) d$S == d$lR
ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
dat <- forstmann
des <- design(data = dat, model = FRQ, matchfun = matchfun,
              formula = list(alpha ~ E, beta ~ 1, h ~ lM, tau ~ lM, t0 ~ 1),
              contrasts = list(h = list(lM = ADmat), tau = list(lM = ADmat)),
              report_p_vector = FALSE,LT=.25,UT=1.5)



# Sampled scale: alpha/beta/tau/t0 are log, h is probit.  These are the
# "typical" row above, with the correct/error split from scenario 3.
pv <- c(alpha = log(2), beta = log(3),
        h = qnorm(0.93), h_lMd = 1.2,          # h: 0.83 error -> 0.99 correct
        tau = log(0.35), tau_lMd = -0.45,      # tau: 0.44 error -> 0.28 correct
        t0 = log(0.20))
pv <- pv[names(sampled_pars(des))]
prior_FRQ <- prior(des, mu_mean = pv)
cat("sampled_pars scale:\n"); print(round(pv, 3))
cat("\nmapped to the natural + generative scale:\n")
print(mapped_pars(des, pv, digits = 4))

sim <- make_data(pv, design = des, n_trials = 4000)
cat(sprintf("\nsimulated %d trials | omissions %.4f | accuracy %.3f\n",
            nrow(sim), mean(!is.finite(sim$rt)),
            mean(sim$R[is.finite(sim$rt)] == sim$S[is.finite(sim$rt)])))
ok <- is.finite(sim$rt)
cat("RT quantiles by accuracy:\n")
print(round(do.call(rbind, tapply(sim$rt[ok], sim$R[ok] == sim$S[ok],
                                  quantile, c(.1, .5, .9))), 3))

cat("\nStart here: alpha ~ 2, beta ~ 3, t0 ~ 0.2, tau ~ 0.3 with the design on\n",
    "h and tau, and let the shapes stay at ~1 df each.  See ?FRQ on why\n",
    "crossing alpha/beta with the same factors as tau gives a ridge.\n", sep = "")
emc = make_emc(dat,des,prior=prior_FRQ)
emc = fit(emc,cores_per_chain = 10, cores_for_chains = 3)
pred = predict(emc)
plot_cdf(dat,pred,defective_factor = "R",factors = c("S","E"))
