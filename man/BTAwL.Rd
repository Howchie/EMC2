% Generated for the BTAwL ballistic transient accumulator models.
\name{BTAwL}
\alias{BTAwL}
\alias{BTAwL_mixed}
\alias{BTAwLMix}
\title{Ballistic Transient Accumulator with Leak}
\description{
The ballistic analogue of Smith's transient accumulator with leak.  The
transient-only member uses a trialwise launch strength and an onset-evoked
pulse.  \code{BTAwL_mixed} adds a shared trialwise strength for Smith's
sustained and transient channels, with mixture weight \code{pi}.
}
\usage{
BTAwL(posdrift = TRUE, drift_distribution = c("normal", "lognormal"),
      mixture = FALSE, chart = c("endpoint", "rate"))

BTAwL_mixed(posdrift = TRUE,
            drift_distribution = c("normal", "lognormal"),
            chart = c("endpoint", "rate"))

BTAwLMix(posdrift = TRUE,
         drift_distribution = c("normal", "lognormal"),
         chart = c("endpoint", "rate"))
}
\arguments{
\item{posdrift}{Logical; truncate normal launch strengths below zero.}
\item{drift_distribution}{Launch distribution, either \code{"normal"} or
\code{"lognormal"}.}
\item{mixture}{If \code{TRUE}, \code{"shared"}, or \code{"mixed"}, return
the shared-strength sustained/transient extension.}

\item{chart}{Clearance chart. \code{"endpoint"} samples the observable hard
transient endpoint \code{Ttrans}; \code{"rate"} samples the intrinsic time
constant \code{tau}.}
}
\details{
The transient-only model solves
\deqn{dX/du = V (u/tau) exp(-u/tau) - k X}{dX/du = V (u/tau) exp(-u/tau) - kX}
with a uniform launch point.  The mixed model uses
\deqn{V_S = pi V, V_T = (1-pi)V}{V_S = pi V, V_T = (1-pi)V}.
At \code{pi = 0} the mixed kernel is exactly the transient-only kernel;
\code{pi = 1} is the sustained-only member.
}
\value{A race-model specification suitable for \code{design}.}
\seealso{\code{BAwL}, \code{BAwD}}
