% Generated for the BTAwL ballistic transient/sustained local race.
\name{BTAwL}
\alias{BTAwL}
\title{Ballistic transient/sustained accumulator race with leak}
\usage{
BTAwL(posdrift = TRUE, drift_distribution = c("normal", "lognormal"),
      chart = c("endpoint", "rate"))
}
\arguments{
\item{posdrift}{Logical; truncate normal launch strengths below zero.}
\item{drift_distribution}{Launch distribution, either \code{"normal"} or
\code{"lognormal"}.}
\item{chart}{Clearance chart. \code{"endpoint"} samples the transient
endpoint \code{Ttrans}; \code{"rate"} samples its time constant
\code{tau_t}.}
}
\description{
The full BTAwL model is a local race between independent sustained and
transient processes. Each process has its own start point and launch draw.
The single base launch distribution is allocated by \code{pi}: the normal
launch location and scale are both multiplied by \code{pi} or \code{1 - pi},
whereas the lognormal log-location shifts by the corresponding log weight.
}
\details{
The transient process follows
\deqn{dX_T/du = V_T (u/\tau_t) \exp(-u/\tau_t) - k X_T}{dX_T/du = V_T (u/\tau_t) exp(-u/\tau_t) - k X_T}
and the sustained process follows
\deqn{dX_S/du = V_S (1 - \exp(-u/\tau_s)) - k X_S}{dX_S/du = V_S (1 - exp(-u/\tau_s)) - k X_S}.
The response is the first local crossing of either process. The likelihood and
simulator use the same independent-launch local-race composition. At
\code{pi = 0} the model is transient-only; at \code{pi = 1} it is sustained-only.
}
\value{A race-model specification suitable for \code{design}.}
\seealso{\code{BTAwLTransient}, \code{BTAwLSustained}, \code{BAwL}}
