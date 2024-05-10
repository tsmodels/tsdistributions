#' Specification of a semi-parametric distribution model
#'
#' @param y a numeric vector
#' @param lower the probability for the lower GPD tail.
#' @param upper the probability for the upper GPD tail.
#' @param kernel_type the choice of the kernel to use from the \code{\link{bkde}}
#' function.
#' @param ... not currently used
#' @returns An object of class \dQuote{tsdistribution.spd_spec}.
#' @examples
#' spec <- spd_modelspec(rnorm(1000))
#' @export spd_modelspec
#' @export
#'
#'
spd_modelspec <- function(y, lower = 0.1, upper = 0.9, kernel_type = c("normal","box","epanech","biweight","triweight"), ...)
{
    setup <- list()
    setup$target$y <- as.numeric(y)
    setup$distribution <- "spd"
    setup$gpd$lower <- lower
    setup$gpd$upper <- upper
    setup$kernel$type <- kernel_type[1]
    n <- length(y)
    x <- sort(as.numeric(y))
    lower_point <- trunc(n * lower)
    upper_point <- trunc(n * (1 - upper))
    upper_quantile <- x[n - upper_point]
    lower_quantile <- x[lower_point]
    upper_init_pars <- .spd_init_parameters(x, upper_quantile)
    lower_init_pars <- .spd_init_parameters(-1.0 * x, -1.0 * lower_quantile)
    paramatrix <- rbind(
        data.table(parameter = "lower_scale", value = lower_init_pars[1], lower = 1e-12, upper = 1000, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
        data.table(parameter = "lower_shape", value = lower_init_pars[2], lower = -0.5, upper = 0.5, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
        data.table(parameter = "upper_scale", value = upper_init_pars[1], lower = 1e-12, upper = 1000, include = 1, estimate = 1, equation = "distribution", group = "distribution"),
        data.table(parameter = "upper_shape", value = upper_init_pars[2], lower = -0.5, upper = 0.5, include = 1, estimate = 1, equation = "distribution", group = "distribution"))
    setup$parmatrix <- paramatrix
    setup$gpd$lower_quantile <- lower_quantile
    setup$gpd$upper_quantile <- upper_quantile
    class(setup) <- "tsdistribution.spdspec"
    return(setup)
}

#' Estimates the parameters of a semi-parametric distribution.
#'
#' @param object an object of class \dQuote{tsdistribution.spdspec}.
#' @param method a choice of \dQuote{Grimshaw}, \dQuote{obre} or \dQuote{nlm} from 
#' \code{\link[mev]{fit.gpd}} or \dQuote{pwm} for the probability weighted moments estimator.
#' @param ... additional parameters passed to the gpd estimation function.
#' @returns An object of class \dQuote{tsdistribution.spdestimate} with slots for
#' the upper, lower and interior kernel fitted values.
#' @details The estimation defaults to the Probability Weighted Moments (pwm) of 
#' Hosking (1985), and alternative methods are provided via the \dQuote{mev} package. 
#' For the interior of the distribution, the \code{\link[KernSmooth]{bkde}} function is used
#' to calculate the kernel density.
#' @method estimate tsdistribution.spdspec
#' @aliases estimate.tsdistribution.spdspec
#' @references 
#' \insertRef{Hosking1985}{tsdistributions}
#' @rdname estimate.tsdistribution.spdspec
#' @export
#'
#'
estimate.tsdistribution.spdspec <- function(object, method = "pwm", ...)
{
    value <- parameter <- NULL
    # first the gpd fit
    if (method != "pwm") {
        # Numerical optimization of the generalized Pareto distribution for data **exceeding** threshold
        # so we reflect for the lower tail
        lt <- fit.gpd(-object$target$y, threshold = -object$gpd$lower_quantile, method = method, ...)
        ut <- fit.gpd(object$target$y, threshold = object$gpd$upper_quantile, method = method, ...)
        lower_tail <- list(pars = lt$param, hessian = solve(lt$vcov), 
                           threshold = object$gpd$lower_quantile, 
                           solution = lt)
        upper_tail <- list(pars = ut$param, hessian = solve(ut$vcov), 
                           threshold = object$gpd$upper_quantile, 
                           solution = ut)
    } else {
        lower_tail <- .gpd_pwm(object, type = "lower")
        upper_tail <- .gpd_pwm(object, type = "upper")
    }
    # then the kernel fit
    x <- sort(object$target$y)
    n <- length(x)
    kernelfit <- bkde(x, object$kernel$type, gridsize = as.integer(n), range.x = c(1.5 * min(x), 1.5 * max(x)))
    sol <- list()
    sol$gpd$lower <- lower_tail
    sol$gpd$lower$prob <- object$gpd$lower
    sol$gpd$upper <- upper_tail
    sol$gpd$upper$prob <- object$gpd$upper
    sol$kernel$type <- object$kernel$type
    sol$kernel$fit <- kernelfit
    sol$parmatrix <- copy(object$parmatrix)
    sol$parmatrix[parameter == "lower_scale", value := lower_tail$pars[1]]
    sol$parmatrix[parameter == "lower_shape", value := lower_tail$pars[2]]
    sol$parmatrix[parameter == "upper_scale", value := upper_tail$pars[1]]
    sol$parmatrix[parameter == "upper_shape", value := upper_tail$pars[2]]
    sol$target$y <- object$target$y
    sol$nobs <- length(x)
    sol$df <- sum(sol$parmatrix$estimate)
    class(sol) <- "tsdistribution.spdestimate"
    sol$loglik <- sum(-dspd(sol$target$y, sol, log = TRUE))
    return(sol)
}


#' @method coef tsdistribution.spdestimate
#' @aliases coef
#' @rdname coef.tsdistribution.estimate
#' @export
#'
#'
coef.tsdistribution.spdestimate <- function(object, ...)
{
    estimate <- NULL
    cf <- c(object$parmatrix[estimate == 1]$value)
    names(cf) <- object$parmatrix[estimate == 1]$parameter
    return(cf)
}

#' Summary of estimated SPD distribution
#'
#' @param object an object of class \dQuote{tsdistribution.spdestimate}.
#' @param ... additional parameters passed to the summary method.
#' @returns A list of summary statistics of the fitted model given in object.
#' @details
#' The standard errors assume a blog diagonal covariance structure between
#' the upper and lower Generalized Pareto Tails.
#' @method summary tsdistribution.spdestimate
#' @aliases summary.tsdistribution.spdestimate
#' @rdname summary.tsdistribution.spdestimate
#' @export
#'
#'
summary.tsdistribution.spdestimate <- function(object, ...)
{
    value <- NULL
    estimate <- NULL
    nobs <- object$nobs
    df <- object$df
    V <- try(vcov(object), silent = TRUE)
    est <- object$parmatrix[estimate == 1]$value
    if (inherits(V, 'try-error')) {
        V <- matrix(NaN, ncol = df, nrow = df)
        se <- rep(NaN, df)
        tval <- rep(NaN, df)
        pval <- rep(NaN, df)
    } else {
        se <- sqrt(diag(V))
        tval <- est / se
        pval <- 2*(1 - pnorm(abs(tval)))
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    coefficients <- as.data.frame(cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = pval))
    par_names <- gsub("lower_","[lower tail] ", par_names)
    par_names <- gsub("upper_","[upper tail] ", par_names)
    rownames(coefficients) <- par_names
    distribution <- "spd"
    coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
    setnames(coefficients, "rn","term")
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = -object$loglik, n_obs = nobs, n_parameters = df,
                AIC = AIC(object),
                BIC = BIC(object))
    class(out) <- "summary.spd"
    return(out)
}

#' @aliases print.summary.tsdistribution
#' @method print summary.spd
#' @rdname print.summary.tsdistribution
#' @export
#'
print.summary.spd <- function(x, digits = max(3L, getOption("digits") - 3L), 
                              signif.stars = getOption("show.signif.stars"),  
                              table.caption = paste0(toupper(x$distribution)," Model Summary\n"), ...)
{
    term <- NULL
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
    cat("dof:", as.integer(x$n_parameters))
    cat(",  ")
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = digits)))
    cat("\n")
    invisible(x)
}

#' @aliases AIC
#' @method AIC tsdistribution.spdestimate
#' @rdname AIC.tsdistribution.estimate
#' @export
#'
#'
AIC.tsdistribution.spdestimate <- function(object, ..., k = 2)
{
    out <- ( -2.0 * as.numeric(object$loglik) + k * object$df)
    return(out)
}

#' @aliases BIC
#' @method BIC tsdistribution.spdestimate
#' @rdname BIC.tsdistribution.estimate
#' @export
#'
#'
BIC.tsdistribution.spdestimate <- function(object, ...)
{
    out <- -2 * as.numeric(object$loglik) + object$df * log(object$nobs)
    return(out)
}

#' @method logLik tsdistribution.spdestimate
#' @aliases logLik
#' @rdname logLik.tsdistribution.estimate
#' @export
#'
#'
logLik.tsdistribution.spdestimate <- function(object, ...)
{
    out <- -object$loglik
    attr(out,"nobs") <- object$nobs
    attr(out,"df") <- object$df
    class(out) <- "logLik"
    return(out)
}



#' @method vcov tsdistribution.spdestimate
#' @aliases vcov
#' @rdname vcov.tsdistribution.estimate
#' @export
#'
vcov.tsdistribution.spdestimate <- function(object, ...)
{
    V <- solve(bread(object))
    par_names <- object$parmatrix[estimate == 1]$parameter
    colnames(V) <- rownames(V) <- par_names
    return(V)
}

.spd_init_parameters <- function(x, threshold)
{
    exceedances <- x[x > threshold]
    excess <- exceedances - threshold
    xbar <- mean(excess)
    sigma_sqr <- var(excess)
    shape_0 <- -0.5 * (((xbar * xbar)/sigma_sqr) - 1)
    scale_0 <- 0.5 * xbar * (((xbar * xbar)/sigma_sqr) + 1)
    theta <- c(scale_0, shape_0)
    return(theta)
}


.gpd_pwm <- function(object, type = "lower")
{
    if (type == "lower") {
        x <- -object$target$y
        threshold <- -object$gpd$lower_quantile
        original_threshold <- object$gpd$lower_quantile
    } else {
        x <- object$target$y
        threshold <- object$gpd$upper_quantile
        original_threshold <- object$gpd$upper_quantile
    }
    x <- sort(x)
    exceedances <- x[x > threshold]
    excess <- exceedances - threshold
    n <- length(excess)
    mu <- mean(excess)
    gamma <- -0.35
    delta <- 0
    pvec <- ((1:n) + gamma)/(n + delta)
    a1 <- mean(sort(excess) * (1 - pvec))
    shape <- 2 - mu/(mu - 2 * a1)
    scale <- (2 * mu * a1)/(mu - 2 * a1)
    pars <- c(scale, shape)
    names(pars) <- c("scale", "shape")
    denom <- n * (1 - 2 * shape) * (3 - 2 * shape)
    if (shape > 0.5) {
        denom <- NA
        warning("Asymptotic Standard Errors not available for PWM when shape>0.5.")
    }
    one_one <- (7 - 18 * shape + 11 * shape^2 - 2 * shape^3) * scale^2
    two_two <- (1 - shape) * (1 - shape + 2 * shape^2) * (2 - shape)^2
    cdiag <- scale * (2 - shape) * (2 - 6 * shape + 7 * shape^2 - 2 * shape^3)
    hessian <- solve(matrix(c(one_one, cdiag, cdiag, two_two), 2)/denom)
    return(list(pars = pars, hessian = hessian, threshold = original_threshold))
}

#' Semi-Parametric Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the semi parametric distribution (spd) which has generalized
#' Pareto tails and kernel fitted interior.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param object an object of class \dQuote{tsdistribution.spdestimate} returned
#' from calling \code{\link{estimate.tsdistribution.spdspec}}.
#' @param linear logical, if TRUE (default) interior smoothing function uses linear 
#' interpolation rather than constant.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator.
#' @rdname spd
#' @export
#'
#'
#'
dspd <- function(x, object, linear = TRUE, log = FALSE)
{
    parameter <- NULL
    if (!is(object,"tsdistribution.spdestimate")) stop("\nobject must be of class tsdistribution.spdestimate")
    y <- sort(x)
    n <- length(x)
    values <- vector(length = n, mode = "numeric")
    orig_data <- sort(object$target$y)
    n_orig <- length(orig_data)
    lower_threshold <-  object$gpd$lower$threshold
    upper_threshold <- object$gpd$upper$threshold
    upper_pars <- object$parmatrix[parameter %in% c("upper_scale","upper_shape")]$value
    lower_pars <- object$parmatrix[parameter %in% c("lower_scale","lower_shape")]$value
    # this is the estimate of CDF in the interval
    kernfit <- object$kernel$fit
    kernel_x <- kernfit$x
    kernel_y <- kernfit$y
    kernel_z <- cumsum(kernel_y/sum(kernel_y))
    
    if (any(x > max(kernel_x)) && any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(x), max(x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    
    if (any(x > max(kernel_x)) && !any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(kernel_x), max(x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    
    if (!any(x > max(kernel_x)) && any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(x), max(kernel_x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    
    mid_x <- lower_threshold <= kernel_x & kernel_x <= upper_threshold
    xm <- as.numeric(kernel_x[mid_x])
    nm <- length(xm)
    # This is the count of the number of kernel fitted occurrences < Lower Threshold
    mid_start <- sum(kernel_x < lower_threshold)
    
    # This is the count of the number of kernel fitted occurrences <= Upper Threshold
    mid_end <- sum(kernel_x <= upper_threshold) + 1
    # this is the estimate of PDF in the interval:
    check <- try(.interp_pdf(xm, kernel_x, kernel_z), silent = TRUE)
    if (!inherits(check, "try-error"))
    {
        values[mid_x] <- check
    }
    upper_tail_x <- kernel_x > upper_threshold
    lower_tail_x <- kernel_x < lower_threshold
    if (length(kernel_x[kernel_x > upper_threshold]) != 0) {
        # this is the estimate of F at upper threshold:
        start <- values[mid_x][length(values[mid_x])]
        zz <- dgp(kernel_x[kernel_x > upper_threshold] - upper_threshold, loc = 0, scale = upper_pars[1], shape = upper_pars[2])
        xd <- function(x) (zz/sum(zz) * ((1 - x) * sum(kernel_y)))[1] - start
        ff = uniroot(f = xd, interval = c(0.5, 1))
        values[upper_tail_x] <- zz/sum(zz) * ((1 - ff$root) * sum(kernel_y))
    }
    if (length(kernel_x[kernel_x < lower_threshold]) != 0) {
        # this is the estimate of F at lower threshold:
        zz = dgp((lower_threshold - kernel_x[kernel_x < lower_threshold]), loc = 0, scale = lower_pars[1], shape = lower_pars[2])
        start <- values[mid_x][1]
        xd <- function(x) (zz/sum(zz) * (x * sum(kernel_y)))[length(zz)] - start
        ff <- uniroot(f = xd, interval = c(0, 0.5))
        values[lower_tail_x] <- zz/sum(zz) * (ff$root * sum(kernel_y))
    }
    final_values <- approx(x = kernel_x, y = values, xout = sort(x), method = "linear", ties = "ordered")$y
    retval <- final_values
    retval[sort.list(x)] <- final_values
    if (log) retval <- log(retval)
    return(retval)
}

#'
#' @rdname spd
#' @export
#'
pspd <- function(q, object, linear = TRUE, lower_tail = TRUE)
{
    parameter <- NULL
    if (!is(object, "tsdistribution.spdestimate")) stop("\nobject must be of class tsdistribution.spdestimate")
    x <- sort(q)
    n <- length(x)
    values <- vector(length = n, mode = "numeric")
    orig_data <- sort(as.matrix(object$target$y))
    n_orig <- length(orig_data)
    # define regions of parametric and kernel interaction
    lower_threshold <-  object$gpd$lower$threshold
    upper_threshold <- object$gpd$upper$threshold
    upper_pars <- object$parmatrix[parameter %in% c("upper_scale","upper_shape")]$value
    lower_pars <- object$parmatrix[parameter %in% c("lower_scale","lower_shape")]$value
    # interior coordinates: estimate of CDF in the interval
    kernfit <- object$kernel$fit
    kernel_x <- kernfit$x
    kernel_y <- kernfit$y
    kernel_z <- cumsum(kernel_y/sum(kernel_y))
    if (any(x > max(kernel_x)) && any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(x), max(x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    
    if (any(x > max(kernel_x)) && !any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(kernel_x), max(x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    
    if (!any(x > max(kernel_x)) && any(x < min(kernel_x))) {
        kfit <- bkde(orig_data, object$kernel$type, gridsize = as.integer(length(orig_data)), range.x = c(min(x), max(kernel_x)))
        kernel_x <- kfit$x
        kernel_y <- kfit$y
        kernel_z <- cumsum(kernel_y/sum(kernel_y))
    }
    # interior values
    mid_x <- lower_threshold <= x & x <= upper_threshold
    xm <- as.numeric(x[mid_x])
    nm <- length(xm)
    # The kernel estimation is done on all the data points so we must truncate...
    # This is the count of the number of kernel fitted occurrences < Lower Threshold
    mid_start <- sum(kernel_x < lower_threshold)
    # This is the count of the number of kernel fitted occurrences <= Upper Threshold
    mid_end <- sum(kernel_x <= upper_threshold) + 1
    # This is the kernel fitted data in interior
    mid_pts <- as.numeric(kernel_x[mid_start:mid_end])
    # This is the number of kernel fitted points in between
    n_mid_pts <- length(mid_pts)
    old_ind <- vector(length = nm, mode = "numeric")
    old_ind[1:nm] <- -1
    # Now we find the location of xm (midpoints of p input data) in the midpoints region of the
    # kernel fitted data
    options(show.error.messages = FALSE)
    mid_index <- approx(x = sort(mid_pts), y = 1:n_mid_pts, xout = xm, method = "linear", ties = "ordered")$y + mid_start
    options(show.error.messages = TRUE)
    # kernel_x data corresponding to the p data
    mid_val <- kernel_x[mid_index]
    mid_val[is.na(mid_val)] <- 0
    mid_index_L <- mid_index - 1
    mid_index_L[mid_index_L == 0] <- NA
    # shifted [backwards] data corresponding to the p data for linear interpolation
    mid_val_S <- kernel_x[mid_index_L]
    mid_val_S[is.na(mid_val_S)] <- (x[mid_x])[is.na(mid_val_S)]
    # corresponding cumulative probability kernel_z
    p_mid_val <- kernel_z[mid_index]
    p_mid_val[is.na(p_mid_val)] <- 0
    
    p_mid_val_S <- kernel_z[mid_index_L]
    p_mid_val_S[is.na(p_mid_val_S)] <- 0

    if (linear) {
        values[mid_x] <- p_mid_val_S + (p_mid_val - p_mid_val_S) * (xm - mid_val_S) / (mid_val - mid_val_S)
    } else {
        values[mid_x] <- p_mid_val_S
    }
    # this is the estimate of CDF at upper threshold:
    p_upper <- kernel_z[mid_end - 1] + (kernel_z[mid_end] - kernel_z[mid_end - 1]) * (upper_threshold - kernel_x[mid_end - 1]) / (kernel_x[mid_end] - kernel_x[mid_end - 1])
    p_lower <- kernel_z[mid_start] + (kernel_z[mid_start + 1] - kernel_z[mid_start]) * (lower_threshold - kernel_x[mid_start]) / (kernel_x[mid_start + 1] - kernel_x[mid_start])
    upper_tail_x <- x > upper_threshold
    lower_tail_x <- x < lower_threshold
    upper_scale <- upper_pars[1]
    upper_shape <- upper_pars[2]
    b <- upper_threshold
    if (length(x[x > upper_threshold]) != 0) {
        values[upper_tail_x] <- 1 - (1 - p_upper) * pgp(x[x > upper_threshold] - upper_threshold, loc = 0, scale = upper_scale, shape = upper_shape, lower.tail = FALSE)
        if (upper_shape < 0 & sum(x > b - upper_scale/upper_shape) > 0) {
            values[x > b - upper_scale/upper_shape]  <- 1.0
        }
    }
    # this is the estimate of CDF at lower threshold:
    if (length(x[x < lower_threshold]) != 0)  {
        lower_scale <- lower_pars[1]
        lower_shape <- lower_pars[2]
        b <- lower_threshold
        if (length(x[x < lower_threshold]) == 0) xlin <- 0 else xlin <- x[x < lower_threshold]
        values[lower_tail_x] <- p_lower * pgp(lower_threshold - xlin, loc = 0, scale = lower_scale, shape = lower_shape, lower.tail = FALSE)
        if (lower_shape < 0 & sum(x < b + lower_scale/lower_shape) > 0) {
            values[x < b + lower_scale/lower_shape] <- 0
        }
    }
    out <- values
    out[sort.list(q)] <- out
    if (!lower_tail) out <- 1 - out
    return(out)
}

#'
#' @rdname spd
#' @export
#'
qspd <- function(p, object, linear = TRUE, lower_tail = TRUE)
{
    if (!lower_tail) {
        p <- 1 - p
    }
    parameter <- NULL
    if (!is(object, "tsdistribution.spdestimate")) stop("\nobject must be of class tsdistribution.spdestimate")
    x <- sort(p)
    n <- length(x)
    orig_data <- sort(object$target$y)
    n_orig <- length(orig_data)
    lower_threshold <-  object$gpd$lower$threshold
    upper_threshold <- object$gpd$upper$threshold
    upper_pars <- object$parmatrix[parameter %in% c("upper_scale","upper_shape")]$value
    lower_pars <- object$parmatrix[parameter %in% c("lower_scale","lower_shape")]$value
    # this is the estimate of F in the interval:
    values <- vector(length = n, mode = "numeric")
    valid_points <- x >= 0 & x <= 1
    values[!valid_points] <- NA
    mid_start <- sum(orig_data < lower_threshold)
    mid_end <- sum(orig_data <= upper_threshold) + 1
    prob_points <- ppoints(orig_data)
    p_upper <- prob_points[mid_end - 1] + (prob_points[mid_end] - prob_points[mid_end - 1]) * (upper_threshold - orig_data[mid_end]) / (orig_data[mid_end] - orig_data[mid_end - 1])
    p_lower <- prob_points[mid_start] + (prob_points[mid_start + 1] - prob_points[mid_start]) * (lower_threshold - orig_data[mid_start]) / (orig_data[mid_start + 1] - orig_data[mid_start])
    mid_x <-  p_lower <= x & x <= p_upper
    xm <- as.double(x[mid_x])
    nm <- as.integer(length(xm))
    mid_points <- as.numeric(prob_points[mid_start:mid_end])
    n_mid_points <- length(mid_points)
    oldind <- vector(length = nm, mode = "numeric")
    oldind[1:nm] <- -1
    mid_index <- approx(y = 1:n_mid_points, x = mid_points, xout = xm, method = "constant", ties = "ordered")$y + mid_start
    mid_val <- prob_points[mid_index]
    mid_index_L <- mid_index - 1
    mid_index[mid_index == 0] <- NA
    mid_index[is.na(mid_index)] <- 0
    
    mid_val_S <- prob_points[mid_index_L]
    mid_val_S[is.na(mid_val_S)] <- 0
    
    p_mid_val <- orig_data[mid_index]
    p_mid_val[is.na(p_mid_val)] <- 0
    
    p_mid_val_S <- orig_data[mid_index_L]
    p_mid_val_S[is.na(p_mid_val_S)] <- 0
    
    if (linear) {
        values[mid_x] <- p_mid_val_S + (p_mid_val - p_mid_val_S) * (xm - mid_val_S)/(mid_val - mid_val_S)
    } else {
        values[mid_x] <- p_mid_val_S
    }
    # this is the estimate of F at upper threshold:
    upper_tail_x <- x > p_upper & valid_points
    lower_tail_x <- x < p_lower & valid_points
    if (length(values[upper_tail_x]) != 0) {
        xu <- x[x > p_upper & valid_points]
        values[upper_tail_x] <- upper_threshold + qgp((xu - p_upper)/(1 - p_upper), loc = 0, scale = upper_pars[1], shape = upper_pars[2])
    }
    # this is the estimate of F at lower threshold:
    if (length(values[lower_tail_x]) != 0) {
        xl <- x[x < p_lower & valid_points]
        values[lower_tail_x] <- lower_threshold - qgp(1 - xl/p_lower, loc = 0, scale = lower_pars[1], shape = lower_pars[2])
    }
    retval <- values
    retval[sort.list(p)] <- values
    return(retval)
}

#'
#' @rdname spd
#' @export
#'
rspd <- function(n, object, linear = TRUE)
{
    # inverse transform sampling method
    if (!is(object, "tsdistribution.spdestimate")) stop("\nobject must be of class tsdistribution.spdestimate")
    u <- runif(n)
    r <- qspd(u, object, linear = linear)
    return(r)
}
