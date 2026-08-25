% Generated for the pure sustained member of BTAwL.
\name{BTAwLSustained}
\alias{BTAwLSustained}
\title{Pure sustained BTAwL}
\usage{
BTAwLSustained(posdrift = TRUE, drift_distribution = c("normal", "lognormal"))
}
\arguments{
\item{posdrift}{Logical; truncate normal launch strengths below zero.}
\item{drift_distribution}{Launch distribution, either \code{"normal"} or
\code{"lognormal"}.}
}
\description{A pure sustained accumulator wrapper around the BTAwL sustained
kernel.}
\details{
The sustained wrapper uses the following parameter matrix:
\tabular{llllll}{
   \strong{Parameter} \tab \strong{Transform} \tab \strong{Natural scale} \tab \strong{Default} \tab \strong{Mapping} \tab \strong{Interpretation} \cr
   \emph{mu} \tab identity \tab [-Inf, Inf] \tab 0 \tab  \tab Mean of log launch strength. \cr
   \emph{sigma} \tab log \tab [0, Inf] \tab log(1) \tab  \tab SD of log launch strength. \cr
   \emph{B} \tab log \tab [0, Inf] \tab log(1) \tab \emph{b} = \emph{B} + \emph{A} \tab Distance from the upper start-point range to the threshold. \cr
   \emph{A} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Start-point range. \cr
   \emph{t0} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Non-decision time. \cr
   \emph{k} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Leak rate. \cr
   \emph{tau_s} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Sustained-drive time constant. \cr
}
With \code{drift_distribution = "normal"}, \code{mu} and \code{sigma} are
replaced by \code{v} and \code{sv}; the launch is truncated positive when
\code{posdrift = TRUE}. \code{pContaminant} and \code{pGuess} are optional
nuisance parameters appended by the data pipeline.
}
\value{A race-model specification suitable for \code{design}.}
