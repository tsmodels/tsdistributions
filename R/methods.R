#' Summary of estimated distribution
#'
#' @param object an object of class tsdistribution.estimate.
#' @param digits the number of significant digits to use when printing,.
#' @param vcov_type the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).
#' @param ... additional parameters passed to the summary method.
#' @return The function computes and returns a list of summary statistics of
#' the fitted model given in object.
#' @method summary tsdistribution.estimate
#' @aliases summary
#' @rdname summary
#' @export
#'
#'

summary.tsdistribution.estimate <- function(object, digits = 4, vcov_type = "H", ...)
{
    value <- NULL
    estimate <- NULL
    n_obs <- object$nobs
    n_parameters <- length(coef(object))
    V <- try(vcov(object, type = vcov_type), silent = TRUE)
    est <- object$parmatrix[estimate == 1]$value
    if (inherits(V, 'try-error')) {
        V <- matrix(NaN, ncol = n_parameters, nrow = n_parameters)
        se <- rep(NaN, n_parameters)
        tval <- rep(NaN, n_parameters)
        pval <- rep(NaN, n_parameters)
    } else {
        se <- sqrt(diag(V))
        tval <- est / se
        pval <- 2*(1 - pnorm(abs(tval)))
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    coefficients <- as.data.frame(cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = pval))
    rownames(coefficients) <- par_names
    loglik <- -object$loglik
    distribution <- object$spec$distribution
    coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
    setnames(coefficients, "rn","term")
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = loglik, n_obs = n_obs, n_parameters = n_parameters,
                AIC = AIC(object),
                BIC = BIC(object))
    class(out) <- "summary.tsdistribution"
    return(out)
}


#' Model Estimation Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsdistribution}
#' @param x an object of class \dQuote{summary.tsdistribution}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param table.caption an optional string for the table caption.
#' @param ... not currently used.
#' @return Console output of the object summary.
#' @aliases print.summary.tsdistribution
#' @method print summary.tsdistribution
#' @rdname print
#' @export
#'
#'
print.summary.tsdistribution <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), 
                                  table.caption = paste0(toupper(x$distribution)," Model Summary\n"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters
    coefs <- copy(x$coefficients)
    coef_names <- coefs$term
    coefs <- coefs[,term := NULL]
    coefs <- as.data.frame(coefs)
    rownames(coefs) <- coef_names
    if (!is.null(table.caption)) cat(table.caption)
    cat("\nCoefficients:\n")
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    cat("\nN:", as.integer(x$n_obs))
    cat(",  ")
    cat("dof:", as.integer(x$n_parameters + 1))
    cat(",  ")
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = digits)))
    cat("\n")
    invisible(x)
}

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsdistribution.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return A numeric value.
#' @aliases AIC
#' @method AIC tsdistribution.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsdistribution.estimate <- function(object, ..., k = 2)
{
    out <- ( -2.0 * as.numeric(object$loglik) + k * object$df)
    return(out)
}

#' Bayesian Information Criterion
#'
#' @description Extract the BIC from an estimated model.
#' @param object an object of class \dQuote{tsdistribution.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases BIC
#' @method BIC tsdistribution.estimate
#' @rdname BIC
#' @export
#'
#'
BIC.tsdistribution.estimate <- function(object, ...)
{
    out <- -2 * as.numeric(object$loglik) + object$df * log(object$nobs)
    return(out)
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
    out <- -object$loglik
    attr(out,"nobs") <- object$nobs
    attr(out,"df") <- object$df
    class(out) <- "logLik"
    return(out)
}


#' The Covariance Matrix of the Estimated Parameters
#'
#' @param object an object of class tsdistribution.estimate
#' @param adjust logical. Should a finite sample adjustment be made? This amounts
#' to multiplication with n/(n-k) where n is the number of observations and k
#' the number of estimated parameters.
#' @param type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @return The variance-covariance matrix of the estimated parameters.
#' @method vcov tsdistribution.estimate
#' @aliases vcov
#' @rdname vcov
#' @export
#'
vcov.tsdistribution.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    type <- match.arg(type[1],c("H","OP","QMLE","NW"))
    N <- nrow(estfun(object))
    if (type == "H") {
        V <- solve(bread(object))
    } else if (type == "QMLE") {
        bread. <- bread(object)
        meat. <- meat_tsdistribution(object, adjust = adjust)
        bread. <- solve(bread.)
        V <- bread. %*% meat. %*% bread.
    } else if (type == "OP") {
        V <- vcovOPG(object, adjust = adjust)
    } else if (type == "NW") {
        bread. <- bread(object)
        meat. <- meatHAC_tsdistribution(object, adjust = adjust, ...)
        bread. <- solve(bread.)
        V <- bread. %*% meat. %*% bread.
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    colnames(V) <- rownames(V) <- par_names
    return(V)
}
