#' Summary of estimated distribution
#'
#' @param object an object of class tsdistribution.estimate.
#' @param digits the number of significant digits to use when printing,.
#' @param ... additional parameters passed to the summary method.
#' @return The function computes and returns a list of summary statistics of
#' the fitted model given in object.
#' @method summary tsdistribution.estimate
#' @aliases summary
#' @rdname summary
#' @export
#'
#'

summary.tsdistribution.estimate <- function(object, digits = 4, ...)
{
 value <- NULL
 estimate <- NULL
 pmatrix <- object$spec$parmatrix
 pmatrix[estimate ==  1, value := object$pars]
 printout <- data.frame("Parameter" = pmatrix[estimate == 1]$parameter, "Est[Value]" = pmatrix[estimate == 1]$value, check.names = FALSE)
 if (!is.null(object$hessian)) {
  S <- try(suppressWarnings(.make_standard_errors(pmatrix, object$hessian)), silent = TRUE)
  if (!inherits(S,'try-error')) {
   printout <- cbind(printout, S)
  }
 }
 print(kable(printout, right = FALSE, digits = digits, row.names = FALSE, format = "simple"))
 return(invisible(printout))
}

#' Extract Model Coefficients
#'
#' @param object an object of class tsdistribution.estimate.
#' @param ... other arguments.
#' @return The function computes and returns a list of summary statistics of
#' the fitted model given in object.
#' @method coef tsdistribution.estimate
#' @aliases coef
#' @rdname coef
#' @export
#'
#'
coef.tsdistribution.estimate <- function(object, ...)
{
 return(object$pars)
}

#' Extract the moments of an estimated distribution
#'
#' @param object an object of class tsdistribution.estimate.
#' @param ... other arguments.
#' @return The function computes and returns a vector of the first four 
#' moments of the distribution based on the estimated parameters. The kurtosis
#' always represents the value in excesss of 3.
#' @method tsmoments tsdistribution.estimate
#' @aliases tsmoments
#' @rdname tsmoments
#' @export
#'
#'
tsmoments.tsdistribution.estimate <- function(object, ...)
{
  p <- c(coef(object),0,0,0)
  s <- dskewness(object$spec$distribution, skew = p[3], shape = p[4], lambda = p[5])
  k <- dkurtosis(object$spec$distribution, skew = p[3], shape = p[4], lambda = p[5])
  return(c(mu = as.numeric(p[1]), sigma = as.numeric(p[2]), skewness = s, kurtosis = k))
}

#' Extract Log-Likelihood
#'
#' @param object an object of class tsdistribution.estimate.
#' @param ... other arguments.
#' @return Returns an object of class logLik. This is a number with at least 
#' one attribute, "df" (degrees of freedom), giving the number of 
#' (estimated) parameters in the model.
#' @method logLik tsdistribution.estimate
#' @aliases logLik
#' @rdname logLik
#' @export
#'
#'
logLik.tsdistribution.estimate <- function(object, ...)
{
  np <- length(object$pars)
  structure(-object$llh, df = np, class = "logLik")
}
