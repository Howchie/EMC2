# Fixed-threshold two-counter race with common-CV shared gamma input draws.

# E[nu1^r nu2^q exp(-s.nu)] under the shared-gamma construction, which is the
# moment the joint density below is built from.  c = 0 is the deterministic-rate
# limit, where the expectation is just the integrand.
.pcorr_m <- function(r, q, s, nu, c, rho) {
  if (c == 0) return(nu[1]^r * nu[2]^q * exp(-sum(s * nu)))
  a  <- 1 / c^2          # gamma shape implied by the common CV
  a0 <- rho * a          # shape carried by the SHARED component
  b  <- a / nu           # rate of each racer's own gamma
  # Rising factorial x(x+1)...(x+n-1); n = 0 is an empty product, and x = 0
  # with n > 0 kills the term (a shape of zero contributes nothing).
  z <- function(x, n) {
    if (!n) 1 else if (!x) 0 else exp(lgamma(x + n) - lgamma(x))
  }
  o <- 0
  for (u in 0:r) for (v in 0:q) {
    o <- o +
      choose(r, u) * choose(q, v) *
      z(a0, u + v) * (1 + s[1] / b[1] + s[2] / b[2])^(-a0 - u - v) *
      z(a - a0, r - u) * (1 + s[1] / b[1])^(-a + a0 - r + u) *
      z(a - a0, q - v) * (1 + s[2] / b[2])^(-a + a0 - q + v)
  }
  o / b[1]^r / b[2]^q
}

# The correlated path is a genuine two-racer joint law, so it needs exactly two
# accumulators, a c and rho shared between them, and the self-excitation and
# overshoot terms switched off.
.pcorr_check <- function(p) {
  n <- c("nu", "c", "gamma", "k", "omega", "t0", "rho")
  if (!is.matrix(p) || nrow(p) != 2 || !all(n %in% colnames(p))) {
    stop("PCOUNTERcorr requires exactly two racers and columns ",
         paste(n, collapse = ", "))
  }
  if (any(!is.finite(p)) ||
      any(p[, "nu"] <= 0) |
      any(p[, "c"] < 0) |
      any(p[, "gamma"] != 0) |
      any(p[, "omega"] != 0) |
      any(p[, "rho"] < 0 | p[, "rho"] > 1) |
      diff(p[, "c"]) != 0 |
      diff(p[, "rho"]) != 0 |
      any(p[, "rho"] > 0 & p[, "c"] == 0)) {
    stop("PCOUNTERcorr requires common c/rho, c>0 if rho>0, gamma=omega=0")
  }
  p
}

#' Correlated gamma-rate Poisson counter race
#' @export
PCOUNTERcorr <- function() {
  req <- c("nu", "c", "gamma", "k", "omega", "t0", "rho")
  list(
    type = "RACE",
    c_name = "PCOUNTERcorr",
    correlated = TRUE,
    correlation_type = "pcounter_shared_gamma",
    p_types = c(nu = log(10), c = log(.2), gamma = log(0), k = log(3),
                omega = log(0), t0 = log(0), rho = qnorm(.5)),
    p_types_canonical = req,
    transform = list(
      func = c(nu = "exp", c = "exp", gamma = "exp", k = "exp",
               omega = "exp", t0 = "exp", rho = "pnorm"),
      lower = c(rho = 0),
      upper = c(rho = 1)
    ),
    bound = list(
      minmax = cbind(nu = c(1e-6, Inf), c = c(1e-4, Inf), gamma = c(1e-4, Inf),
                     k = c(0, Inf), omega = c(1e-4, Inf), t0 = c(0, Inf),
                     rho = c(.001, .999)),
      exception = c(gamma = 0, k = 0, omega = 0, t0 = 0, rho = 0)
    ),
    Ttransform = function(p, d) p,
    rfun = function(data = NULL, pars) {
      lr <- data$lR
      if (length(levels(lr)) != 2 || nrow(pars) %% 2) {
        stop("PCOUNTERcorr requires two racers")
      }
      z <- sapply(seq_len(nrow(pars) / 2), function(j) {
        p <- .pcorr_check(pars[(2 * j - 1):(2 * j), , drop = FALSE])
        a <- 1 / p[1, "c"]^2
        # Shared component first, then the two private ones: the draw ORDER is
        # part of the reproducible stream, so it must not be rearranged.
        u <- if (!is.finite(a)) {
          p[, "nu"]
        } else {
          p[, "nu"] / a *
            ((if (p[1, "rho"] == 0) 0 else rgamma(1, p[1, "rho"] * a)) +
             (if (p[1, "rho"] == 1) c(0, 0) else rgamma(2, (1 - p[1, "rho"]) * a)))
        }
        k <- 2 + floor(p[, "k"] + .5)
        p[, "t0"] + rgamma(2, k, u)
      })
      w <- max.col(-t(z))
      data.frame(R = factor(levels(lr)[w], levels = levels(lr)),
                 rt = z[cbind(w, seq_len(ncol(z)))])
    },
    race_density = dPCOUNTERcorr,
    dfun = function(...) stop("joint density: use dPCOUNTERcorr"),
    pfun = function(...) stop("joint model"),
    log_likelihood = function(...) stop("compiled PCOUNTERcorr path required")
  )
}

#' @export
dPCOUNTERcorr <- function(rt, pars, winner = 1L) {
  p <- .pcorr_check(pars); w <- as.integer(winner); l <- 3L - w
  k <- 2L + floor(p[, "k"] + .5)
  vapply(rt, function(t) {
    x <- t - p[, "t0"]
    if (x[w] <= 0) return(0)
    front <- x[w]^(k[w] - 1L) / factorial(k[w] - 1L)
    # Loser has not started: it contributes no counts, so the sum collapses.
    if (x[l] <= 0) {
      return(front * .pcorr_m(k[w], 0L, c(x[w], 0), p[, "nu"],
                              p[1, "c"], p[1, "rho"]))
    }
    # Otherwise marginalise the loser's count n over the states in which it has
    # not yet reached its own threshold.
    sum(vapply(0:(k[l] - 1L), function(n) {
      r <- if (w == 1L) k[1L] else n
      q <- if (w == 1L) n else k[2L]
      front * x[l]^n / factorial(n) *
        .pcorr_m(r, q, x, p[, "nu"], p[1, "c"], p[1, "rho"])
    }, numeric(1)))
  }, numeric(1))
}
