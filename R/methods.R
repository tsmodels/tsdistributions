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