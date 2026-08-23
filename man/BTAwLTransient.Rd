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
\value{A race-model specification suitable for \code{design}.}
