% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spd.R
\name{spd_modelspec}
\alias{spd_modelspec}
\title{Specification of a semi-parametric distribution model}
\usage{
spd_modelspec(
  y,
  lower = 0.1,
  upper = 0.9,
  kernel_type = c("normal", "box", "epanech", "biweight", "triweight"),
  ...
)
}
\arguments{
\item{y}{a numeric vector}

\item{lower}{the probability for the lower GPD tail.}

\item{upper}{the probability for the upper GPD tail.}

\item{kernel_type}{the choice of the kernel to use from the \code{\link[KernSmooth]{bkde}} 
function.}

\item{...}{not currently used}
}
\value{
An object of class \dQuote{tsdistribution.spd_spec}.
}
\description{
Specification of a semi-parametric distribution model
}
\examples{
spec <- spd_modelspec(rnorm(1000))
}
