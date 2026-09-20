\name{PCOUNTERcorr}
\alias{PCOUNTERcorr}
\alias{dPCOUNTERcorr}
\title{Correlated-rate Poisson counter race}
\description{
  A two-racer Poisson counter race with gamma-distributed trialwise event
  rates.  The two rates have common coefficient of variation and a shared
  gamma component, giving correlation \code{rho}; the compiled likelihood
  evaluates the finite fixed-threshold race density analytically.
}
\usage{
PCOUNTERcorr()
dPCOUNTERcorr(rt, pars, winner = 1L)
}
\arguments{
  \item{rt}{Response times on the observed scale.}
  \item{pars}{A two-row matrix with columns \code{nu}, \code{c},
    \code{gamma}, \code{k}, \code{omega}, \code{t0}, and \code{rho}.}
  \item{winner}{The winning racer, either 1 or 2.}
}
\details{
  The initial implementation is deliberately the fixed-threshold,
  non-reinforcing member: \code{gamma = 0} and \code{omega = 0}.  Both racers
  must share \code{c} and \code{rho}; \code{k} is converted to the integer
  threshold \code{2 + floor(k + 0.5)}.  The shared-rate construction has
  marginal mean \code{nu}, coefficient of variation \code{c}, and rate
  correlation \code{rho}.
}
\value{
  \code{PCOUNTERcorr} returns an EMC2 race-model specification.
  \code{dPCOUNTERcorr} returns the winning-racer density.
}
\seealso{\code{\link{PCOUNTER}}}
\keyword{models}
