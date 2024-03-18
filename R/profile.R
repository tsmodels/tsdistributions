#' Model Parameter Profiling
#'
#' @description Profiles the model parameters under the specified distribution.
#' @details The function profiles the parameters of a model by simulating and
#' then estimating multiple paths from the assumed distribution. This makes it possible
#' to obtain a better understanding of the convergence properties (RMSE) of each
#' parameter under different data sizes.
#' @param object an object of class \dQuote{tsdistribution.spec} with pre-set parameters.
#' @param nsim the number of paths to generate.
#' @param sizes a vector of data sizes for which to simulate and estimate.
#' @param seed an object specifying if and how the random number generator
#' should be initialized. See the simulate documentation for more details.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... not currently used.
#' @note The function can use parallel functionality as long as the user has set
#' up a \code{\link[future]{plan}} using the future package.
#' @returns An object of class \dQuote{tsdistribution.profile}.
#' @aliases tsprofile
#' @method tsprofile tsdistribution.spec
#' @rdname tsprofile.tsdistribution.spec
#' @export
#'
#'
tsprofile.tsdistribution.spec <- function(object, nsim = 100, sizes = c(800, 1000, 1500, 2000, 3000), seed = NULL, trace = FALSE, ...)
{
    ":=" <- NULL
    parameter <- NULL
    actual <- error <- value <- NULL
    .error <- .absolute_error <- .squared_error <- .absolute_percent_error <- NULL
    spec <- copy(object)
    par_names <- object$parmatrix[estimate == 1]$parameter
    true_pars <- object$parmatrix[estimate == 1]$value
    true_skewness <- dskewness(object$distribution, skew = object$parmatrix[parameter == "skew"]$value,
                               shape = object$parmatrix[parameter == "shape"]$value,
                               lambda = object$parmatrix[parameter == "lambda"]$value)
    true_kurtosis <- dkurtosis(object$distribution, skew = object$parmatrix[parameter == "skew"]$value,
                               shape = object$parmatrix[parameter == "shape"]$value,
                               lambda = object$parmatrix[parameter == "lambda"]$value)
    # dkurtosis returns the excess kurtosis...add back 3
    true_kurtosis <- true_kurtosis + 3
    true_pars <- c(true_pars, true_skewness, true_kurtosis)
    par_names <- c(par_names, "skewness","kurtosis")
    true_table <- data.table(parameter = par_names, actual = true_pars)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # sort the sizes
    sizes <- sort(sizes)
    
    sim <- rdist(object$distribution, n = nsim * max(sizes), mu = object$parmatrix[parameter == "mu"]$value, 
                 sigma = object$parmatrix[parameter == "sigma"]$value, skew = object$parmatrix[parameter == "skew"]$value,
                 shape = object$parmatrix[parameter == "shape"]$value, lambda = object$parmatrix[parameter == "lambda"]$value)
    sim <- matrix(sim, ncol = nsim, nrow = max(sizes))
    sim_out <- list()
    if (trace) {
        prog_trace <- progressor(nsim)
    }
    start_time <- Sys.time()
    size_results <- lapply(1:nsim, function(i){
        if (trace) prog_trace()
        sim_out %<-% future_lapply(1:length(sizes), function(j){
            # set the threads to 1 for each fork
            setDTthreads(1)
            spec_new <- distribution_modelspec(sim[1:sizes[j], i], distribution = object$distribution)
            # account for any fixed parameters
            if (any(object$parmatrix$estimate == 0)) {
                fixed_pars <- object$parmatrix[estimate == 0]$parameter
                fixed_pars_value <- object$parmatrix[estimate == 0]$value
                spec_new$parmatrix[parameter %in% fixed_pars, value := fixed_pars_value]
                spec_new$parmatrix[parameter %in% fixed_pars, estimate := 0]
            }
            mod <- try(estimate(spec_new), silent = TRUE)
            if (inherits(mod, 'try-error')) {
                coeff <- spec_new$parmatrix[estimate == 1]$value * as.numeric(NA)
                skewness <- kurtosis <- as.numeric(NA)
            } else{
                # check results and catch errors
                coeff <- unname(coef(mod))
                skewness <- dskewness(object$distribution, skew = mod$parmatrix[parameter == "skew"]$value,
                                      shape = mod$parmatrix[parameter == "shape"]$value,
                                      lambda = mod$parmatrix[parameter == "lambda"]$value)
                kurtosis <- dkurtosis(object$distribution, skew = mod$parmatrix[parameter == "skew"]$value,
                                      shape = mod$parmatrix[parameter == "shape"]$value,
                                      lambda = mod$parmatrix[parameter == "lambda"]$value)
                # dkurtosis returns the excess kurtosis...add back 3
                kurtosis <- kurtosis + 3
            }
            par_names <- spec_new$parmatrix[estimate == 1]$parameter
            return(data.table(parameter = c(par_names, "skewness", "kurtosis"), value = c(coeff, skewness, kurtosis),
                              draw = i, size = sizes[j]))
        },future.packages = c("tsmethods", "tsdistributions", "data.table"))
        sim_out <- eval(sim_out)
        sim_out <- rbindlist(sim_out)
        return(sim_out)
    })
    size_results <- rbindlist(size_results)
    size_results <- merge(size_results, true_table, by = "parameter", all.x = TRUE)
    size_results[, .error := actual - value]
    size_results[, .squared_error := .error^2]
    size_results[, .absolute_error := abs(.error)]
    size_results[, .absolute_percent_error := .absolute_error/abs(actual)]
    summary_table <- size_results[,list(RMSE = sqrt(mean(.squared_error, na.rm = TRUE)), MAE = mean(.absolute_error, na.rm = TRUE), MAPE = mean(.absolute_percent_error, na.rm = TRUE)), by = c('parameter','size')]
    sol <- list()
    sol$profile <- size_results
    sol$summary <- summary_table
    sol$elapsed <- Sys.time() - start_time
    sol$model <- list(distribution = object$distribution, parameters = par_names)
    class(sol) <- 'tsdistribution.profile'
    return(sol)
}



#' Distribution Profile Summary
#'
#' @description Summary method for class \dQuote{tsdistribution.profile}
#' @param object an object of class \dQuote{tsdistribution.profile}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param measure either one of the 3 included measure in the summary slot of
#' the returned object \dQuote{RMSE}, \dQuote{MAE} or \dQuote{MAPE}, else any
#' other user calculated measure which has been generated in the summary table
#' post processing.
#' @param ... not currently used.
#' @returns A list with summary information of class \dQuote{summary.tsdistribution.profile},
#' including a table with each actual parameter against the measure chosen across each
#' size in the profile.
#' @method summary tsdistribution.profile
#' @rdname summary.tsdistribution.profile
#' @export
#'
#'
summary.tsdistribution.profile <- function(object, digits = 4, measure = "RMSE", ...)
{
    actual <- NULL
    valid_measures <- colnames(object$summary)[!colnames(object$summary) %in% c("parameter","size")]
    if (!measure %in% valid_measures) stop("\nmeasure not in summary table!")
    tab <- dcast(object$summary, parameter ~ size, value.var = measure)
    actuals <- object$profile[,list(actual = mean(actual, na.rm = TRUE)), by = "parameter"]
    # re-order rows
    tab <- merge(actuals, tab, by = "parameter")
    tab <- tab[object$model$parameters]
    colnames(tab)[1] <- "parameter/size"
    out <- list(table = tab, distribution = object$model$distribution, measure = measure)
    class(out) <- "summary.tsdistribution.profile"
    return(out)
}

#' Profile Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsdistribution.profile}
#' @param x an object of class \dQuote{summary.tsdistribution.profile}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @returns Invisibly returns the original summary object and prints out to the console.
#' @aliases print.summary.tsdistribution.profile
#' @method print summary.tsdistribution.profile
#' @rdname print.tsdistribution.profile
#' @export
#'
#'
print.summary.tsdistribution.profile <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    distribution <- paste0(toupper(x$distribution))
    cat(paste0("\n",x$measure," Profile Summary : ",distribution))
    cat("\n\n")
    df <- as.data.frame(x$table)
    r_names <- df[,1]
    df <- df[,-1]
    rownames(df) <- r_names
    print(df, digits = digits)
    invisible(x)
}
