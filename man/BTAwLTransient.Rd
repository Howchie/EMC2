% Generated for the pure transient member of BTAwL.
\name{BTAwLTransient}
\alias{BTAwLTransient}
\title{Pure transient BTAwL}
\usage{
BTAwLTransient(posdrift = TRUE, drift_distribution = c("normal", "lognormal"),
               chart = c("endpoint", "rate"))
}
\arguments{
\item{posdrift}{Logical; truncate normal launch strengths below zero.}
\item{drift_distribution}{Launch distribution, either \code{"normal"} or
\code{"lognormal"}.}
\item{chart}{Clearance chart. \code{"endpoint"} samples \code{Ttrans};
\code{"rate"} samples \code{tau}.}
}
\description{A pure transient accumulator wrapper around the BTAwL transient
kernel.}
\details{
The transient wrapper uses the following parameter matrix:
\tabular{llllll}{
   \strong{Parameter} \tab \strong{Transform} \tab \strong{Natural scale} \tab \strong{Default} \tab \strong{Mapping} \tab \strong{Interpretation} \cr
   \emph{mu} \tab identity \tab [-Inf, Inf] \tab 0 \tab  \tab Mean of log launch strength. \cr
   \emph{sigma} \tab log \tab [0, Inf] \tab log(1) \tab  \tab SD of log launch strength. \cr
   \emph{B} \tab log \tab [0, Inf] \tab log(1) \tab \emph{b} = \emph{B} + \emph{A} \tab Distance from the upper start-point range to the threshold. \cr
   \emph{A} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Start-point range. \cr
   \emph{t0} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Non-decision time. \cr
   \emph{k} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Leak rate. \cr
   \emph{Ttrans} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Transient-drive endpoint parameter. \cr
}
With \code{drift_distribution = "normal"}, \code{mu} and \code{sigma} are
replaced by \code{v} and \code{sv}; the launch is truncated positive when
\code{posdrift = TRUE}. \code{pContaminant} and \code{pGuess} are optional
nuisance parameters appended by the data pipeline.
}
\value{A race-model specification suitable for \code{design}.}
