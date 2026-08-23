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
\value{A race-model specification suitable for \code{design}.}
