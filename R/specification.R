#' Specification of distribution model
#'
#' @param y a numeric vector
#' @param distribution the type of distribution. Valid choices are norm (Normal),
#' snorm (Skew Normal), std (Student), sstd (Skew Student), ged (Generalized Error),
#' sged (Skew Generalized Error), nig (Normal Inverse Gaussian), gh (Generalized Hyperbolic),
#' ghst (Generalized Hyperbolic Skew Student) and jsu (Johnson's SU).
#' @details
#' All distributions are parameterized in terms of their mean (\sQuote{mu}), standard deviation
#' \sQuote{sigma}, skew \sQuote{skew} and shape \sQuote{shape} parameters. Additionally,
#' for the Generalized Hyperbolic distribution, there is an extra shape parameter
#' \dQuote{lambda} arising from the GIG mixing distribution.
#' Parameters can be fixed post initialization by setting setting specific values
#' to the \sQuote{value} column in the parmatrix table and setting the \sQuote{estimate}
#' variable to 0 (instead of 1).
#' @param ... not currently used
#' @returns An object of class \dQuote{tsdistribution.spec}
#' @examples
#' spec <- distribution_modelspec(rnorm(1000), distribution = "gh")
#' # fix lambda and shape
#' spec$parmatrix[parameter == 'lambda', value := 30]
#' spec$parmatrix[parameter == 'lambda', estimate := 0]
#' @export distribution_modelspec
#' @export
#'
#'
#'
distribution_modelspec <- function(y, distribution = "norm", ...)
{
    distribution <- match.arg(distribution[1], valid_distributions())
    setup <- list()
    setup$target$y <- as.numeric(y)
    setup$distribution <- distribution
    setup$parmatrix <- .distribution_spec(distribution, y, model = "distribution", ...)
    class(setup) <- "tsdistribution.spec"
    return(setup)
}

.distribution_spec <- function(distribution, y, model = "arima")
{
    sig <- sd(y, na.rm = TRUE)
    mu <- mean(y, na.rm = TRUE)
    if (model == "distribution") {
        tmp <- rbind(
            data.table(parameter = "mu", value = mu, lower = -Inf, upper = Inf, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
            data.table(parameter = "sigma", value = sig, lower = 1e-14, upper = 10 * sig, include = 1, estimate = 1, equation = "distribution", group = "distribution"))
    } else if (model == "arima") {
        tmp <- data.table(parameter = "sigma", value = sig, lower = 1e-14, upper = 10 * sig, include = 1, estimate = 1, equation = "distribution", group = "distribution")
    } else {
        tmp <- NULL
    }
    if (distribution == "norm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, include = 0, estimate = 0, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, include = 0, estimate = 0, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "ged") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, include = 0, estimate = 0, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 2, lower = 0.1, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "std") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, include = 0, estimate = 0, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "snorm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.5, lower = 0.1, upper = 10, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, include = 0, estimate = 0, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "sged") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 2, lower = 0.1, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "sstd") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "nig") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 0.4, lower = 0.01, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "gh") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                    data.table(parameter = "shape", value = 2, lower = 0.25, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                    data.table(parameter = "lambda", value = -0.5, lower = -30, upper = 30, include = 1, estimate = 1, equation = "distribution", group = "distribution"))
        return(tmp)
    }
    if (distribution == "jsu") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = -20, upper = 20, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 1, lower = 0.1, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
    } 
    if (distribution == "ghst") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -80, upper = 80, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "shape", value = 8.2, lower = 4.01, upper = 100, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, include = 0, estimate = 0, equation = "distribution", group = "distribution"))
        return(tmp)
 }
}
