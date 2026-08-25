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

The model uses the following parameter matrix. \code{B} is the distance from
the upper end of the start-point range to the threshold, so the kernel receives
\code{b = B + A}.
\tabular{llllll}{
   \strong{Parameter} \tab \strong{Transform} \tab \strong{Natural scale} \tab \strong{Default} \tab \strong{Mapping} \tab \strong{Interpretation} \cr
   \emph{mu} \tab identity \tab [-Inf, Inf] \tab 0 \tab  \tab Mean of log launch strength for a lognormal launch. \cr
   \emph{sigma} \tab log \tab [0, Inf] \tab log(1) \tab  \tab SD of log launch strength for a lognormal launch. \cr
   \emph{B} \tab log \tab [0, Inf] \tab log(1) \tab \emph{b} = \emph{B} + \emph{A} \tab Distance from the upper start-point range to the threshold. \cr
   \emph{A} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Start-point range. \cr
   \emph{t0} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Non-decision time. \cr
   \emph{k} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Leak rate. \cr
   \emph{tau_s} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Sustained-drive time constant. \cr
   \emph{Ttrans} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Transient-drive endpoint parameter. \cr
   \emph{pi} \tab probit \tab [0, 1] \tab qnorm(.5) \tab  \tab Probability allocated to the sustained process. \cr
}

With \code{drift_distribution = "normal"}, \code{mu} and \code{sigma} are
replaced by \code{v} and \code{sv}; the launch is truncated positive when
\code{posdrift = TRUE}. \code{pContaminant} and \code{pGuess} are optional
nuisance parameters appended by the data pipeline.

The evidence scale is not intrinsic: multiplying the launch strength, threshold
distance, and start-point range by the same positive constant leaves response
times unchanged. Fix one scale parameter (normally \code{sv} for a normal
launch, or one intercept among \code{mu}, \code{B}, and \code{A} for a
lognormal launch) when specifying a design.
}
\value{A race-model specification suitable for \code{design}.}
\seealso{\code{BTAwLTransient}, \code{BTAwLSustained}, \code{BAwL}}
