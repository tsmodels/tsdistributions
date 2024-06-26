% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spd.R
\name{estimate.tsdistribution.spdspec}
\alias{estimate.tsdistribution.spdspec}
\title{Estimates the parameters of a semi-parametric distribution.}
\usage{
\method{estimate}{tsdistribution.spdspec}(object, method = "pwm", ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsdistribution.spdspec}.}

\item{method}{a choice of \dQuote{Grimshaw}, \dQuote{obre} or \dQuote{nlm} from 
\code{\link[mev]{fit.gpd}} or \dQuote{pwm} for the probability weighted moments estimator.}

\item{...}{additional parameters passed to the gpd estimation function.}
}
\value{
An object of class \dQuote{tsdistribution.spdestimate} with slots for
the upper, lower and interior kernel fitted values.
}
\description{
Estimates the parameters of a semi-parametric distribution.
}
\details{
The estimation defaults to the Probability Weighted Moments (pwm) of 
Hosking (1985), and alternative methods are provided via the \dQuote{mev} package. 
For the interior of the distribution, the \code{\link[KernSmooth]{bkde}} function is used
to calculate the kernel density.
}
\references{
\insertRef{Hosking1985}{tsdistributions}
}
