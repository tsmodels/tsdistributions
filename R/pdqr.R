# This file contains the location-scale invariant parameterization of the
# distributions used in the tsdistributions package.
# ------------------------------------------------------------------------------
# Skew Generalized Hyberolic Student's T
# alpha = abs(beta)+1e-12, lambda = -nu/2
# ------------------------------------------------------------------------------
# Location-Scale Invariant Parametrization
.paramGHST <- function(betabar, nu){
    # Alexios Galanos 2012
    # betabar is the skew parameter = beta*delta (parametrization 4 in Prause)
    # nu is the shape parameter
    # details can be found in the vignette
    delta <- (((2 * betabar^2) / ((nu - 2) * (nu - 2) * (nu - 4))) + (1/(nu - 2)))^(-0.5)
    beta <- betabar / delta
    mu <- -((beta * (delta^2)) / (nu - 2))
    return(c(mu, delta, beta, nu))
}

#' Generalized Hyperbolic Skewed Student Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the generalized hyperbolic skew student distribution parameterized in 
#' terms of mean, standard deviation, skew and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname ghst
#' @export
#'
#'
#'
dghst <- function(x, mu = 0, sigma = 1, skew = 1, shape = 8, log = FALSE)
{
    if (any(abs(skew) < 1e-12)) skew[which(abs(skew) < 1e-12)] <- 1e-12
    val_length <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
    max_n = max(val_length)
    if (val_length[1] != max_n) x <- rep(x[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    if (val_length[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dghst", x = as.double(x), mu = as.double(mu),
          sigma = as.double(sigma), skew = as.double(skew),
          shape = as.double(shape), ans = ans, n = as.integer(max_n),
          logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#'
#' @rdname ghst
#' @export
#'
rghst <- function(n, mu = 0, sigma = 1, skew = 1, shape = 8)
{
    if (any(abs(skew) < 1e-12)) skew[which(abs(skew) < 1e-12)] <- 1e-12
    val_length = c(length(mu), length(sigma), length(skew), length(shape))
    if (val_length[1] != n) mu <- rep(mu[1], n)
    if (val_length[2] != n) sigma <- rep(sigma[1], n)
    if (val_length[3] != n) skew <- rep(skew[1], n)
    if (val_length[4] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rghst", n = as.integer(n), mu = as.double(mu),
              sigma = as.double(sigma), skew = as.double(skew),
              shape = as.double(shape), ans = ans,
              PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#'
#' @rdname ghst
#' @export
#'
pghst <- function(q, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE)
{
    if (any(abs(skew) < 1e-12)) skew[which(abs(skew) < 1e-12)] <- 1e-12
    val_length <- c(length(q), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(val_length)
    if (val_length[1] != max_n) q <- rep(q[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    if (val_length[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    for (i in 1:max_n) {
        ans[i] <- .pghst(q[i], mu = mu[i], sigma = sigma[i], skew = skew[i],
                          shape = shape[i], lower_tail = lower_tail, log = log)
    }
    return(ans)
}

.pghst <- function(q, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE)
{
    param <- .paramGHST(skew, shape)
    # scale the parameters
    mux <- param[1] * sigma + mu
    delta <- param[2] * sigma
    beta <- param[3] / sigma
    nu <- param[4]
    ans <- pskewhyp(q, mu = mux, delta = delta, beta = beta, nu = nu, log.p = log, lower.tail = lower_tail)
    return(ans)
}

#'
#' @rdname ghst
#' @export
#'
qghst <- function(p, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE) {
    if (any(abs(skew) < 1e-12)) skew[which(abs(skew) < 1e-12)] <- 1e-12
    val_length <- c(length(p), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(val_length)
    if (val_length[1] != max_n) p <- rep(p[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    if (val_length[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    for (i in 1:max_n) {
        ans[i] <- .qghst(p[i], mu = mu[i], sigma = sigma[i], skew = skew[i], shape = shape[i], lower_tail = lower_tail, log = log)
    }
    return(ans)
}

.qghst <- function(p, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) {
        p <- 1 - p
        lower_tail <- TRUE
    }
    param <- .paramGHST(skew, shape)
    # scale the parameters
    mu <- param[1] * sigma + mu
    delta <- param[2] * sigma
    beta <- param[3] / sigma
    nu <- param[4]
    ans <- qskewhyp(p, mu = mu, delta = delta, beta = beta, nu = nu,
                    lower.tail = lower_tail, log.p = log, method = c("spline","integrate")[1], 
                    nInterpol = 1000)
    return(ans)
}

#' Skewed Normal Distribution of Fernandez and Steel
#'
#' @description Density, distribution, quantile function and random number
#' generation for the skewed normal distribution parameterized in 
#' terms of mean, standard deviation and skew parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname snorm
#' @export
#'
#'
#'
dsnorm <- function(x, mu = 0, sigma = 1, skew = 1.5, log = FALSE)
{
    val_length <- c(length(x), length(mu), length(sigma), length(skew))
    max_n <- max(val_length)
    if (val_length[1] != max_n) x <- rep(x[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dsnorm", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew), ans = ans,
                  n = as.integer(max_n), logr = as.integer(log),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname snorm
#' @export
psnorm <- function(q, mu = 0, sigma = 1, skew = 1.5, lower_tail = TRUE, log = FALSE)
{
    val_length = c(length(q), length(mu), length(sigma), length(skew))
    max_n = max(val_length)
    if (val_length[1] != max_n) q <- rep(q[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_psnorm", q = as.double(q), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew), ans = ans,
                  n = as.integer(max_n), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname snorm
#' @export
qsnorm <- function(p, mu = 0, sigma = 1, skew = 1.5, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) {
        p <- 1 - p
    }
    val_length <- c(length(p), length(mu), length(sigma), length(skew))
    max_n <- max(val_length)
    if (val_length[1] != max_n) p <- rep(p[1], max_n)
    if (val_length[2] != max_n) mu <- rep(mu[1], max_n)
    if (val_length[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (val_length[4] != max_n) skew <- rep(skew[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_qsnorm", p = as.double(p), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew), ans = ans,
                  n = as.integer(max_n), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname snorm
#' @export
rsnorm <- function(n, mu = 0, sigma = 1, skew = 1.5)
{
    val_length = c(length(mu), length(sigma), length(skew))
    if (val_length[1] != n) mu <- rep(mu[1], n)
    if (val_length[2] != n) sigma <- rep(sigma[1], n)
    if (val_length[3] != n) skew <- rep(skew[1], n)
    ans <- double(n)
    sol <- try(.C("c_rsnorm", n = as.integer(n), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew), ans = ans,
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' Generalized Error Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the generalized error distribution parameterized in 
#' terms of mean, standard deviation and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname ged
#' @export
#'
#'
#'
dged <- function(x, mu = 0, sigma = 1, shape = 2, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dged", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n), logr = as.integer(log),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname ged
#' @export
pged <- function(q, mu = 0, sigma = 1, shape = 2, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_pged", q = as.double(q), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname ged
#' @export
qged <- function(p, mu = 0, sigma = 1, shape = 2, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_qged", p = as.double(p), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname ged
#' @export
rged <- function(n, mu = 0, sigma = 1, shape = 2)
{
    value_len <- c(length(mu), length(sigma), length(shape))
    if (value_len[1] != n) mu <- rep(mu[1], n)
    if (value_len[2] != n) sigma <- rep(sigma[1], n)
    if (value_len[3] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rged", n = as.integer(n), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' Skewed Generalized Error Distribution of Fernandez and Steel
#'
#' @description Density, distribution, quantile function and random number
#' generation for the skewed generalized error distribution parameterized in 
#' terms of mean, standard deviation, skew and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname sged
#' @export
#'
#'
#'
dsged <- function(x, mu = 0, sigma = 1, skew = 1.5, shape = 2, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dsged", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape), ans = ans, n = as.integer(max_n),
                  logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname sged
#' @export
psged <- function(q, mu = 0, sigma = 1, skew = 1.5, shape = 2, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_psged", q = as.double(q), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape), ans = ans,
                  n = as.integer(max_n), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname sged
#' @export
qsged <- function(p, mu = 0, sigma = 1, skew = 1.5, shape = 2, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol = try(.C("c_qsged", p = as.double(p), mu = as.double(mu),
                 sigma = as.double(sigma), skew = as.double(skew),
                 shape = as.double(shape), ans = ans, n = as.integer(max_n),
                 PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname sged
#' @export
rsged <- function(n, mu = 0, sigma = 1, skew = 1.5, shape = 2)
{
    value_len <- c(length(mu), length(sigma), length(skew), length(shape))
    if (value_len[1] != n) mu <- rep(mu[1], n)
    if (value_len[2] != n) sigma <- rep(sigma[1], n)
    if (value_len[3] != n) skew <- rep(skew[1], n)
    if (value_len[4] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rsged", n = as.integer(n), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape),
                  ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

Heaviside <- function(x, a = 0)
{
    result <- (sign(x - a) + 1)/2
    return(result)
}

#' Student Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the student distribution parameterized in terms of mean,
#' standard deviation and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname std
#' @export
#'
#'
#'
dstd <- function(x, mu = 0, sigma = 1, shape = 5, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dstd", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n), logr = as.integer(log),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname std
#' @export
pstd <- function(q, mu = 0, sigma = 1, shape = 5, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_pstd", q = as.double(q), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname std
#' @export
qstd <- function(p, mu = 0, sigma = 1, shape = 5, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(mu), length(sigma), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_qstd", p = as.double(p), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, n = as.integer(max_n),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname std
#' @export
rstd <- function(n, mu = 0, sigma = 1, shape = 5)
{
    value_len <- c(length(mu), length(sigma), length(shape))
    if (value_len[1] != n) mu <- rep(mu[1], n)
    if (value_len[2] != n) sigma <- rep(sigma[1], n)
    if (value_len[3] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rstd", n = as.integer(n), mu = as.double(mu),
                  sigma = as.double(sigma), shape = as.double(shape),
                  ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' Skewed Student Distribution of Fernandez and Steel
#'
#' @description Density, distribution, quantile function and random number
#' generation for the skewed student distribution parameterized in 
#' terms of mean, standard deviation, skew and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param mu mean.
#' @param n number of observations.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname sstd
#' @export
#'
#'
#'
dsstd <- function(x, mu = 0, sigma = 1, skew = 1.5, shape = 5, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dsstd", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape), ans = ans, n = as.integer(max_n),
                  logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname sstd
#' @export
psstd <- function(q, mu = 0, sigma = 1, skew = 1.5, shape = 5, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol = try(.C("c_psstd", q = as.double(q), mu = as.double(mu),
                 sigma = as.double(sigma), skew = as.double(skew),
                 shape = as.double(shape), ans = ans, n = as.integer(max_n),
                 PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname sstd
#' @export
qsstd <- function(p, mu = 0, sigma = 1, skew = 1.5, shape = 5, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_qsstd", p = as.double(p), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape), ans = ans, n = as.integer(max_n),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname sstd
#' @export
rsstd <- function(n, mu = 0, sigma = 1, skew = 1.5, shape = 5)
{
    value_len <- c(length(mu), length(sigma), length(skew), length(shape))
    if (value_len[1] != n) mu <- rep(mu[1], n)
    if (value_len[2] != n) sigma <- rep(sigma[1], n)
    if (value_len[3] != n) skew <- rep(skew[1], n)
    if (value_len[4] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rsstd", n = as.integer(n), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew),
                  shape = as.double(shape),
                  ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
      return(sol)
    } else {
      return(sol$ans)
  }
}
# ------------------------------------------------------------------------------
# Generalized Hyperbolic Distribution (standard representation)
# ------------------------------------------------------------------------------
.unitroot <- function(f, interval, lower = min(interval), upper = max(interval), tol = .Machine$double.eps^0.25, ...)
{

    if (is.null(args(f))) {
    if (f(lower) * f(upper) >= 0) return(NA)
    } else {
        if (f(lower, ...) * f(upper, ...) >= 0) return(NA)
    }
    ans <- uniroot(f = f, interval = interval, lower = lower, upper = upper, tol = tol, ...)
    return(ans$root)
}

.kappaGH <- function(x, lambda = 1)
{
    # A function implemented by Diethelm Wuertz
    #   Returns modified Bessel function ratio
    stopifnot(x >= 0)
    stopifnot(length(lambda) == 1)
    if (lambda == -0.5) {
    # NIG:
        kappa <- 1/x
    } else {
    # GH:
        kappa <- (besselK(x, lambda + 1, expon.scaled = TRUE)/besselK(x, lambda, expon.scaled = TRUE)) / x
    }
    return(kappa)
}

.deltaKappaGH <- function(x, lambda = 1)
{
 return(.kappaGH(x, lambda + 1) - .kappaGH(x, lambda))
}

.paramGH <- function(rho = 0, zeta = 1, lambda = 1)
{
    #   Change parameterizations to alpha(zeta, rho, lambda)
    Rho2 <- 1 - rho^2
    alpha <- zeta^2 * .kappaGH(zeta, lambda) / Rho2
    alpha <- alpha * (1 + rho^2 * zeta^2 * .deltaKappaGH(zeta, lambda) / Rho2)
    alpha <- sqrt(alpha)
    beta <- alpha * rho
    delta <- zeta / (alpha * sqrt(Rho2))
    mu <- -beta * delta^2 * .kappaGH(zeta, lambda)
    return(c(alpha = alpha, beta = beta, delta = delta, mu = mu))
}

dghyp <- function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
{
    value_len <- c(length(x), length(alpha), length(beta), length(delta), length(mu), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) alpha <- rep(alpha[1], max_n)
    if (value_len[3] != max_n) beta <- rep(beta[1], max_n)
    if (value_len[4] != max_n) delta <- rep(delta[1], max_n)
    if (value_len[5] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    if (any(alpha <= 0)) stop("alpha must be greater than zero")
    if (any(delta <= 0)) stop("delta must be greater than zero")
    if (any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
    ans <- double(max_n)
    sol <- try(.C("c_dghyp", x = as.double(x), alpha = as.double(alpha),
                  beta = as.double(beta), delta = as.double(delta),
                  mu = as.double(mu), lambda = as.double(lambda),
                  ans = ans, n = as.integer(max_n), logr = as.integer(log),
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

pghyp <- function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(alpha), length(beta), length(delta), length(mu), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) alpha <- rep(alpha[1], max_n)
    if (value_len[3] != max_n) beta <- rep(beta[1], max_n)
    if (value_len[4] != max_n) delta <- rep(delta[1], max_n)
    if (value_len[5] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    if (any(alpha <= 0)) stop("alpha must be greater than zero")
    if (any(delta <= 0)) stop("delta must be greater than zero")
    if (any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
    ans <- rep(NA, max_n)
    for (i in 1:max_n) {
        Integral = integrate(dghyp, -Inf, q[i], stop.on.error = FALSE, alpha = alpha[i], 
                             beta = beta[i], delta = delta[i], mu = mu[i], 
                             lambda = lambda[i])
        ans[i] = as.numeric(unlist(Integral)[1])
    }
    if (!lower_tail) ans <- 1 - ans
    if (log) ans <- log(ans)
    return(ans)
}

qghyp <- function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(alpha), length(beta), length(delta), length(mu), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) alpha <- rep(alpha[1], max_n)
    if (value_len[3] != max_n) beta <- rep(beta[1], max_n)
    if (value_len[4] != max_n) delta <- rep(delta[1], max_n)
    if (value_len[5] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    if (any(alpha <= 0)) stop("alpha must be greater than zero")
    if (any(delta <= 0)) stop("delta must be greater than zero")
    if (any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
    # Internal Function:
    .froot <- function(x, alpha, beta, delta, mu, lambda, p)
    {
        pghyp(q = x, alpha = alpha, beta = beta, delta = delta, mu = mu, lambda = lambda) - p
    }
    # Quantiles:
    ans <- rep(NA, max_n)
    for (i in 1:max_n) {
        lower <- -1
        upper <-  1
        counter <- 0
        iteration <- NA
        while (is.na(iteration)) {
            iteration = .unitroot(f = .froot, interval = c(lower,upper), alpha = alpha[i], 
                                  beta = beta[i], delta = delta[i], mu = mu[i], 
                                  lambda = lambda[i], p = p[i])
            counter <- counter + 1
            lower <- lower - 2^counter
            upper <- upper + 2^counter
        }
        ans[i] <- iteration
    }
    if (log) ans <- log(ans)
    return(ans)
}


#' Generalized Hyperbolic Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the generalized hyperbolic distribution parameterized in 
#' terms of mean, standard deviation, skew and two shape parameters (shape and
#' lambda)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param lambda additional shape parameter determining subfamilies of this
#' distributions.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @return d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname gh
#' @export
#'
#'
#'
dgh <- function(x, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    ans <- double(max_n)
    sol = try(.C("c_dgh", x = as.double(x), mu = as.double(mu), 
                 sigma = as.double(sigma), skew = as.double(skew), 
                 shape = as.double(shape), lambda = as.double(lambda), 
                 ans = ans, n = as.integer(max_n), logr = as.integer(log), 
                 PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname gh
#' @export
pgh <- function(q, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(skew), length(shape), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    tmp <- t(apply(cbind(skew, shape, lambda), 1, function(x) .paramGH(x[1], x[2], x[3])))
    ans <- pghyp((q - mu)/sigma, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda, lower_tail = lower_tail, log = log)
    # equivalent: pghyp(q, alpha/sigma, beta/sigma, delta*sigma, mu*sigma+mu, lambda)
    return(as.numeric(ans))
}

#' @rdname gh
#' @export
qgh <- function(p, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(p), length(mu), length(sigma), length(skew), length(shape), length(lambda))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    if (value_len[6] != max_n) lambda <- rep(lambda[1], max_n)
    tmp <- t(apply(cbind(skew, shape, lambda), 1, function(x) .paramGH(x[1], x[2], x[3])))
    ans <- mu + sigma * as.numeric(qghyp(p, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda, lower_tail = lower_tail, log = log))
    return(ans)
}

#' @rdname gh
#' @export
rgh <- function(n, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1)
{
    params <- .paramGH(rho = skew, zeta = shape, lambda = lambda)
    out <- rghyp(n = n, mu = params[4], delta = params[3], alpha = params[1], beta = params[2], lambda = lambda)
    out <- mu + sigma * out
    return(out)
}

# ------------------------------------------------------------------------------
# Normal Inverse Gaussian (NIG) Distribution
# ------------------------------------------------------------------------------

#' Normal Inverse Gaussian Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for the normal inverse gaussian distribution generalized parameterized in 
#' terms of mean, standard deviation, skew and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname nig
#' @export
#'
#'
#'
dnig <- function(x, mu = 0, sigma = 1, skew = 0, shape = 1, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_dnig", x = as.double(x), mu = as.double(mu),
                  sigma = as.double(sigma), skew = as.double(skew), 
                  shape = as.double(shape), ans = ans, n = as.integer(max_n), 
                  logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname nig
#' @export
pnig <- function(q, mu = 0, sigma = 1, skew = 0, shape = 1, lower_tail = TRUE, log = FALSE)
{
    return(pgh(q, mu, sigma, skew, shape, lambda = -0.5, lower_tail = lower_tail, log = log))
}

#' @rdname nig
#' @export
qnig <- function(p, mu = 0, sigma = 1, skew = 0, shape = 1, lower_tail = TRUE, log = FALSE)
{
    return(qgh(p, mu, sigma, skew, shape, lambda = -0.5, lower_tail = lower_tail, log = log))
}

#' @rdname nig
#' @export
rnig <- function(n, mu = 0, sigma = 1, skew = 0, shape = 1)
{
    return(rgh(n, mu, sigma, skew, shape, lambda = -0.5))
}


#' Johnson's SU Distribution
#'
#' @description Density, distribution, quantile function and random number
#' generation for Johnson's SU distribution parameterized in 
#' terms of mean, standard deviation, skew and shape parameters.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname jsu
#' @export
#'
#'
#'
djsu <- function(x, mu = 0, sigma = 1, skew = 1, shape = 0.5, log = FALSE)
{
    value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) x <- rep(x[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_djsu", x = as.double(x), mu = as.double(mu), 
                  sigma = as.double(sigma), skew = as.double(skew), 
                  shape = as.double(shape), ans = ans, n = as.integer(max_n), 
                  logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

#' @rdname jsu
#' @export
pjsu <- function(q, mu = 0, sigma = 1, skew = 1, shape = 0.5, lower_tail = TRUE, log = FALSE)
{
    value_len <- c(length(q), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) q <- rep(q[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_pjsu", q = as.double(q), mu = as.double(mu), 
                  sigma = as.double(sigma), skew = as.double(skew), 
                  shape = as.double(shape), ans = ans, n = as.integer(max_n), 
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (!lower_tail) sol$ans <- 1 - sol$ans
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname jsu
#' @export
qjsu <- function(p, mu = 0, sigma = 1, skew = 1, shape = 0.5, lower_tail = TRUE, log = FALSE)
{
    if (!lower_tail) p <- 1 - p
    value_len <- c(length(p), length(mu), length(sigma), length(skew), length(shape))
    max_n <- max(value_len)
    if (value_len[1] != max_n) p <- rep(p[1], max_n)
    if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
    if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
    if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
    if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
    ans <- double(max_n)
    sol <- try(.C("c_qjsu", p = as.double(p), mu = as.double(mu), 
                  sigma = as.double(sigma), skew = as.double(skew), 
                  shape = as.double(shape), ans = ans, n = as.integer(max_n), 
                  PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        if (log) sol$ans <- log(sol$ans)
        return(sol$ans)
    }
}

#' @rdname jsu
#' @export
rjsu <- function(n, mu = 0, sigma = 1, skew = 1, shape = 0.5)
{
    value_len <- c(length(mu), length(sigma), length(skew), length(shape))
    if (value_len[1] != n) mu <- rep(mu[1], n)
    if (value_len[2] != n) sigma <- rep(sigma[1], n)
    if (value_len[3] != n) skew <- rep(skew[1], n)
    if (value_len[4] != n) shape <- rep(shape[1], n)
    ans <- double(n)
    sol <- try(.C("c_rjsu", n = as.integer(n), mu = as.double(mu), 
                  sigma = as.double(sigma), skew = as.double(skew), 
                  shape = as.double(shape), 
                  ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
    if (inherits(sol, 'try-error')) {
        return(sol)
    } else {
        return(sol$ans)
    }
}

# ------------------------------------------------------------------------------
# Distribution Wrapper Functions

.nigtransform <- function(rho, zeta)
{
    nigpars <- t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = -0.5)))
    colnames(nigpars) <- c("alpha", "beta", "delta", "mu")
    return(nigpars)
}

.ghyptransform <- function(rho, zeta, lambda)
{
    n <- length(zeta)
    ghyppars <- t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
    ghyppars <- cbind(ghyppars, rep(lambda, n))
    colnames(ghyppars) <- c("alpha", "beta", "delta", "mu", "lambda")
    return(ghyppars)
}

.nigscale <- function(mu, sigma, skew, shape)
{
    nigpars <- t(apply(cbind(skew, shape), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = -0.5)))
    xdensity <- matrix(0, ncol = 4, nrow = length(sigma))
    # alpha, beta, delta, mu
    xdensity[,4] <- nigpars[,1]/sigma
    xdensity[,3] <- nigpars[,2]/sigma
    xdensity[,2] <- nigpars[,3]*sigma
    xdensity[,1] <- nigpars[,4]*sigma + mu
    # technically: mu, delta, beta, alpha
    colnames(xdensity) <- c("mu", "delta", "beta", "alpha")
    return(xdensity)
}

.ghstscale <- function(mu, sigma, skew, shape)
{
    ghpars <- t(apply(cbind(skew, shape), 1, FUN = function(x) .paramGHST(betabar = x[1], nu = x[2])))
    xdensity <- matrix(0, ncol = 5, nrow = length(sigma))
    # alpha, beta, delta, mu
    xdensity[,5] <- -shape/2
    xdensity[,4] <- (abs(ghpars[,3]) + 1e-12)/sigma
    xdensity[,3] <- ghpars[,3]/sigma
    xdensity[,2] <- ghpars[,2]*sigma
    xdensity[,1] <- ghpars[,1]*sigma + mu
    # technically: mu, delta, beta, alpha
    colnames(xdensity) <- c("mu", "delta", "beta", "alpha", "lambda")
    return(xdensity)
}

.ghypscale <- function(mu, sigma, skew, shape, lambda)
{
    if (length(lambda) > 1) { 
        ghpars <- t(apply(cbind(skew, shape, lambda), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = x[3])))
    } else { 
        ghpars <- t(apply(cbind(skew, shape), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
    }
    xdensity <- matrix(0, ncol = 4, nrow = length(sigma))
    # alpha, beta, delta, mu
    xdensity[,4] <- ghpars[,1]/sigma
    xdensity[,3] <- ghpars[,2]/sigma
    xdensity[,2] <- ghpars[,3]*sigma
    xdensity[,1] <- ghpars[,4]*sigma + mu
    # technically: mu, delta, beta, alpha
    colnames(xdensity) <- c("mu", "delta", "beta", "alpha")
    return(xdensity)
}


#' Parameter Transformation
#' @description Transforms parameters from standardized representation to distribution
#' specific representation for the nig and gh distributions.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param lambda additional shape parameter for the Generalized Hyperbolic
#' distribution.
#' @returns The (alpha, beta, delta, mu) representation.
#' @rdname ghyptransform
#' @export
#'
#'
#'
nigtransform <- function(mu = 0, sigma = 1,  skew = 0, shape = 3)
{
    return(.nigscale(mu = mu, sigma = sigma, skew = skew, shape = shape))
}

#' @rdname ghyptransform
#' @export
ghyptransform <- function(mu = 0, sigma = 1,  skew = 0, shape = 3, lambda = -0.5)
{
    return(.ghypscale(mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda))
}


#' Distributions pqdr wrapper
#'
#' @description Density, distribution, quantile function and random number
#' generation for all the distributions in the package.
#' @param distribution a valid distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape  shape parameter.
#' @param lambda additional shape parameter for the Generalized Hyperbolic
#' distribution.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @returns d gives the density, p gives the distribution function, q gives the quantile function
#' and r generates random deviates. Output depends on x or q length, or n for the random number
#' generator
#' @rdname ddist
#' @export
#'
ddist <- function(distribution = "norm", x, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, log = FALSE)
{
    distribution <- match.arg(distribution[1], valid_distributions())
    ans <- switch(distribution,
                  "norm" = dnorm(x, mean = mu, sd = sigma, log = log),
                  "snorm" = dsnorm(x, mu = mu, sigma = sigma, skew = skew, log = log), 
                  "std" = dstd(x, mu = mu, sigma = sigma, shape = shape, log = log), 
                  "sstd" = dsstd(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log), 
                  "ged" = dged(x, mu = mu, sigma = sigma, shape = shape, log = log), 
                  "sged" = dsged(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log),
                  "nig" = dnig(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log), 
                  "gh" = dgh(x, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda, log = log), 
                  "jsu" = djsu(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log), 
                  "ghst" = dghst(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log)
                  )
    return(ans)
}

#' @rdname ddist
#' @export
pdist <- function(distribution = "norm", q, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, lower_tail = TRUE, log = FALSE)
{
    distribution <- match.arg(distribution[1], valid_distributions())
    ans <- switch(distribution, 
                  "norm" = pnorm(q, mean = mu, sd = sigma, lower.tail = lower_tail, log.p = log), 
                  "snorm" = psnorm(q, mu = mu, sigma = sigma, skew = skew, lower_tail = lower_tail, log = log), 
                  "std" = pstd(q, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log), 
                  "sstd" = psstd(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "ged" = pged(q, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log), 
                  "sged" = psged(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "nig" = pnig(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "gh" = pgh(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda, lower_tail = lower_tail, log = log), 
                  "jsu" = pjsu(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "ghst" = pghst(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log)
                  )
    return(ans)
}

#' @rdname ddist
#' @export
qdist <- function(distribution = "norm", p, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, lower_tail = TRUE, log = FALSE)
{
    distribution <- match.arg(distribution[1], valid_distributions())
    ans <- switch(distribution, 
                  "norm" = qnorm(p, mean = mu, sd = sigma, lower.tail = lower_tail, log.p = log), 
                  "snorm" = qsnorm(p, mu = mu, sigma = sigma, skew = skew, lower_tail = lower_tail, log = log), 
                  "std" = qstd(p, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log), 
                  "sstd" = qsstd(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "ged" = qged(p, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log), 
                  "sged" = qsged(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
                  "nig" = qnig(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "gh" = qgh(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda, lower_tail = lower_tail, log = log), 
                  "jsu" = qjsu(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log), 
                  "ghst" = qghst(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
                  )
    return(ans)
}

#' @rdname ddist
#' @export
rdist <- function(distribution = "norm", n, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5)
{
    distribution <- match.arg(distribution[1], valid_distributions())
    ans <- switch(distribution, 
                  "norm" = rnorm(n, mean = mu, sd = sigma), 
                  "snorm" = rsnorm(n, mu = mu, sigma = sigma, skew = skew), 
                  "std" = rstd(n, mu = mu, sigma = sigma, shape = shape), 
                  "sstd" = rsstd(n, mu = mu, sigma = sigma, skew = skew, shape = shape), 
                  "ged" = rged(n, mu = mu, sigma = sigma, shape = shape), 
                  "sged" = rsged(n, mu = mu, sigma = sigma, skew = skew, shape = shape), 
                  "nig" =  rnig(n, mu = mu, sigma = sigma, skew = skew, shape = shape), 
                  "gh" = rgh(n, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda), 
                  "jsu" = rjsu(n, mu = mu, sigma = sigma, skew = skew, shape = shape), 
                  "ghst" = rghst(n, mu = mu, sigma = sigma, skew = skew, shape = shape)
                  )
    return(ans)
}
