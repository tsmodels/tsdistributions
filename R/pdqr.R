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
 delta <- ( ((2 * betabar^2)/((nu - 2)*(nu - 2)*(nu - 4))) + (1/(nu - 2)) )^(-0.5)
 beta <- betabar/delta
 mu <- -( (beta * (delta^2))/(nu - 2))
 return( c(mu, delta, beta, nu) )
}

#' Skewed Generalized Hyperbolic Student Distribution
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of the 
#' generalized hyperbolic student distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
#' and r generates random deviates. Output depends on x or q length, or n for the random number 
#' generator
#' @rdname ghst
#' @export
#'
#'
#'
dsghst <- function(x, mu = 0, sigma = 1, skew = 1, shape = 8, log = FALSE)
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
rsghst <- function(n, mu = 0, sigma = 1, skew = 1, shape = 8)
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
psghst <- function(q, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE)
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
  ans[i] <- .psghst(q[i], mu = mu[i], sigma = sigma[i], skew = skew[i], 
                   shape = shape[i], lower_tail = lower_tail, log = log)
 }
 return(ans)
}

.psghst <- function(q, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE) 
{
 param <- .paramGHST(skew, shape)
 # scale the parameters
 mux <- param[1]*sigma + mu
 delta <- param[2]*sigma
 beta <- param[3]/sigma
 nu <- param[4]
 ans <- pskewhyp(q, mu = mux, delta = delta, beta = beta, nu = nu, log.p = log, lower.tail = lower_tail)
 return(ans)
}

#'
#' @rdname ghst
#' @export
#'
qsghst <- function(p, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE) {
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
  ans[i] <- .qsghst(p[i], mu = mu[i], sigma = sigma[i], skew = skew[i], shape = shape[i], lower_tail = lower_tail, log = log)
 }
 return( ans )
}

.qsghst <- function(p, mu = 0, sigma = 1, skew = 1, shape = 8, lower_tail = TRUE, log = FALSE)
{
 if (!lower_tail) {
  p <- 1 - p
  lower_tail <- TRUE
 }
 param <- .paramGHST(skew, shape)
 # scale the parameters
 mu <- param[1]*sigma + mu
 delta <- param[2]*sigma
 beta <- param[3]/sigma
 nu <- param[4]
 ans <- qskewhyp(p, mu = mu, delta = delta, beta = beta, nu = nu, 
                 lower.tail = lower_tail, log.p = log, method = c("spline","integrate")[1])
 return(ans)
}

#' Skewed Normal Distribution of Fernandez and Steel
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of the 
#' skewed normal distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
#' generation for the location scale invariant parameterization of the 
#' generalized error distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
pged <- function(q, mu = 0, sigma = 1, shape = 2, lower_tail = TRUE, log = TRUE)
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
#' generation for the location scale invariant parameterization of the 
#' skewed generalized error distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
################################################################################
Heaviside <- function(x, a = 0) 
{   
 result <- (sign(x - a) + 1)/2
 return(result)
}
# ------------------------------------------------------------------------------
# Student Distribution
# ------------------------------------------------------------------------------

#' Student Distribution
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of the 
#' student distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
#' generation for the location scale invariant parameterization of the 
#' skewed student distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param mu mean.
#' @param n number of observations.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
  kappa <- (besselK(x, lambda + 1, expon.scaled = TRUE)/besselK(x, lambda, expon.scaled = TRUE) ) / x
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
 alpha <- alpha * ( 1 + rho^2 * zeta^2 * .deltaKappaGH(zeta, lambda) / Rho2)
 alpha <- sqrt(alpha)  
 beta <- alpha * rho
 delta <- zeta / ( alpha * sqrt(Rho2) )
 mu <- -beta * delta^2 * .kappaGH(zeta, lambda)
 return(c(alpha = alpha, beta = beta, delta = delta, mu = mu))
}

.rgigjd <- function(n, theta)
{	
 #	Original Version by David Scott
 lambda <- theta[1]
 chi <- theta[2]
 psi <- theta[3]
 if (chi < 0) stop("chi can not be negative")
 if (psi < 0) stop("psi can not be negative")
 if ((lambda >= 0) & (psi == 0)) stop("When lambda >= 0, psi must be > 0")
 if ((lambda <= 0) & (chi == 0)) stop("When lambda <= 0, chi must be > 0")
 if (chi == 0) stop("chi = 0, use rgamma")
 if (psi == 0) stop("algorithm only valid for psi > 0")	
 alpha <- sqrt(psi/chi)
 beta <- sqrt(psi*chi)
 m <- (lambda - 1 + sqrt((lambda - 1)^2 + beta^2))/beta
 g <- function(y) {
  0.5 * beta * y^3 - y^2 * (0.5 * beta * m + lambda + 1) + y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m
 }
 upper <- m
 while (g(upper) <= 0) upper = 2 * upper
 yM <- uniroot(g, interval = c(0, m))$root
 yP <- uniroot(g, interval = c(m, upper))$root
 a <- (yP - m) * (yP/m)^(0.5 * (lambda - 1)) * exp(-0.25 * beta * (yP + 1/yP - m - 1/m))
 b <- (yM - m) * (yM/m)^(0.5 * (lambda - 1)) * exp(-0.25 * beta * (yM + 1/yM - m - 1/m))
 c <- -0.25 * beta * (m + 1/m) + 0.5 * (lambda - 1) * log(m)	
 output <- numeric(n)
 for (i in 1:n) {
  need.value <- TRUE
  while (need.value == TRUE) {
   R1 <- runif(1)
   R2 <- runif(1)
   Y <- m + a * R2/R1 + b * (1 - R2)/R1
   if (Y > 0) {
    if (-log(R1) >= -0.5 * (lambda - 1) * log(Y) + 0.25 * beta * (Y + 1/Y) + c) {
     need.value <- FALSE
    }
   }
  }
  output[i] <- Y
 }
 return(output/alpha)
}

.rgigjd1 <- function(n, theta)
{	
 #	Original Version by David Scott
 if (length(theta) == 2) theta <- c(1, theta)
 lambda <- 1
 chi <- theta[2]
 psi <- theta[3]
 if (chi < 0) stop("chi can not be negative") 
 if (psi < 0) stop("psi can not be negative")	
 if (chi == 0) stop("chi = 0, use rgamma")
 if (psi == 0) stop("When lambda >= 0, psi must be > 0")	
 alpha <- sqrt(psi/chi)
 beta <- sqrt(psi*chi)
 m <- abs(beta)/beta
 g <- function(y) {
  0.5 * beta * y^3 - y^2 * (0.5 * beta * m + lambda + 1) + y * (-0.5 * beta) + 0.5 * beta * m
 }
 upper <- m
 while (g(upper) <= 0) upper <- 2 * upper
 yM <- uniroot(g,interval = c(0, m))$root
 yP <- uniroot(g,interval = c(m, upper))$root
 a <- (yP - m) * exp(-0.25 * beta * (yP + 1/yP - m - 1/m))
 b <- (yM - m) * exp(-0.25 * beta * (yM + 1/yM - m - 1/m))
 c <- -0.25 * beta * (m + 1/m)
 output <- numeric(n)
 for (i in 1:n) {
  need.value <- TRUE
  while (need.value == TRUE) {
   R1 <- runif(1)
   R2 <- runif(1)
   Y <- m + a * R2/R1 + b * (1 - R2)/R1
   if (Y > 0) {
    if (-log(R1) >= 0.25 * beta * (Y + 1/Y) + c) {
     need.value <- FALSE
    }
   }
  }
  output[i] <- Y
 }	
 return(output/alpha)
}

dgh <- function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
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
 sol <- try(.C("c_dgh", x = as.double(x), alpha = as.double(alpha), 
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

pgh <- function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, lower_tail = TRUE, log = FALSE)
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
  Integral = integrate(dgh, -Inf, q[i], stop.on.error = FALSE, alpha = alpha[i], 
                       beta = beta[i], delta = delta[i], mu = mu[i],
                       lambda = lambda[i])
  ans[i] = as.numeric(unlist(Integral)[1])
 }
 if (!lower_tail) ans <- 1 - ans
 if (log) ans <- log(ans)
 return( ans )
}

qgh <- function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, lower_tail = TRUE, log = FALSE)
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
  pgh(q = x, alpha = alpha, beta = beta, delta = delta,
      mu = mu, lambda = lambda) - p
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


rgh <- function(n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
 value_len <- c(length(alpha), length(beta), length(delta), length(mu), length(lambda))
 max_n <- n
 if (value_len[1] != max_n) alpha <- rep(alpha[1], max_n)
 if (value_len[2] != max_n) beta <- rep(beta[1], max_n)
 if (value_len[3] != max_n) delta <- rep(delta[1], max_n)
 if (value_len[4] != max_n) mu <- rep(mu[1], max_n)
 if (value_len[5] != max_n) lambda <- rep(lambda[1], max_n)
 chi <- delta^2
 psi <- alpha^2 - beta^2
 if (any(alpha <= 0)) stop("alpha must be greater than zero")
 if (any(delta <= 0)) stop("delta must be greater than zero")
 if (any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
 V <- cbind(lambda, chi, psi)
 X <- apply(V, 1, function(x){
  ifelse(x[1] == 1, 
         .rgigjd1(1, c(x[1], x[2], x[3])), 
         .rgigjd(1,  c(x[1], x[2], x[3])))
 })
 sigma <- sqrt(as.numeric(X))
 Z <- rnorm(n)
 Y <- mu + beta * sigma^2 + sigma*Z
 return(Y)
}

#' (Standardized) Generalized Hyperbolic Distribution
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of the 
#' standardized generalized hyperbolic distribution.
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
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
#' and r generates random deviates. Output depends on x or q length, or n for the random number 
#' generator
#' @rdname sgh
#' @export
#'
#'
#'
dsgh <- function(x, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, log = FALSE)
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
 sol = try(.C("c_dghyp", x = as.double(x), mu = as.double(mu), 
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

#' @rdname sgh
#' @export
psgh <- function(q, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, lower_tail = TRUE, log = FALSE) 
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
 ans <- pgh((q - mu)/sigma, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda, lower_tail = lower_tail, log = log)
 # equivalent: pgh(q, alpha/sigma, beta/sigma, delta*sigma, mu*sigma+mu, lambda)
 return(as.numeric(ans))
}

#' @rdname sgh
#' @export
qsgh <- function(p, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, lower_tail = TRUE, log = FALSE) 
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
 ans <- mu + sigma*as.numeric( qgh(p, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda, lower_tail = lower_tail, log = log) )
 return( ans )
}

#' @rdname sgh
#' @export
rsgh <- function(n, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1)
{
 value_len <- c(length(mu), length(sigma), length(skew), length(shape), length(lambda))
 if (value_len[1] != n) mu <- rep(mu[1], n)
 if (value_len[2] != n) sigma <- rep(sigma[1], n)
 if (value_len[3] != n) skew <- rep(skew[1], n)
 if (value_len[4] != n) shape <- rep(shape[1], n)
 if (value_len[5] != n) lambda <- rep(lambda[1], n)
 ans <- double(n)
 sol <- try(.C("c_rghyp", n = as.integer(n), mu = as.double(mu), 
              sigma = as.double(sigma), skew = as.double(skew), 
              shape = as.double(shape), lambda = as.double(lambda),
              ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
 if (inherits(sol, 'try-error')) {
  return(sol)
 } else {
  return(sol$ans)
 }
}
# ------------------------------------------------------------------------------
# Normal Inverse Gaussian (NIG) Distribution
# ------------------------------------------------------------------------------
.qnigC <- function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lower_tail = TRUE, log = FALSE)
{   
 if (alpha <= 0) stop("Invalid parameters: alpha <= 0.\n")
 if (alpha^2 <= beta^2) stop("Invalid parameters: alpha^2 <= beta^2.\n")
 if (delta <= 0) stop("Invalid parameters: delta <= 0.\n")
 if ((sum(is.na(p)) > 0)) {
   stop("Invalid probabilities:\n", p,"\n")
 } else {
  if (sum(p < 0) + sum(p > 1) > 0) stop("Invalid probabilities:\n", p,"\n")
 }
 if (!lower_tail) p <- 1 - p
 n <- length(p)
 q <- rep(0, n)
 ans <- .C("qNIG",
                 p = as.double(.CArrange(p,1,1,n)),
                 i_mu = as.double(mu),
                 i_delta = as.double(delta),
                 i_alpha = as.double(alpha),
                 i_beta = as.double(beta),
                 i_n = as.integer(n),
                 q = as.double(.CArrange(q, 1, 1, n)), PACKAGE = "tsdistributions")
 ans <- ans[[7]]
 ans[ans <= -1.78e+308] <- -Inf
 ans[ans >= 1.78e+308] <- Inf
 if (log) ans <- log(ans)
 return(ans)
}


.CArrange <- function(obj, i, j, n)
{
 # Description:
 #   Arrange input matrices and vectors in a suitable way for the C program
 #   Matrices are transposed because the C program stores matrices by row 
 #   while R stores matrices by column
 # Arguments:
 #   i - length of first dimension
 #   j - length of second dimension
 #   n - length of third dimension
 # Value:
 #   out - transformed data set
 # Author: 
 #   Daniel Berg <daniel at nr.no> (Kjersti Aas <Kjersti.Aas at nr.no>)
 #   Date: 12 May 2005
 #   Version: 1.0.2	
 if (is.null(obj)) stop("Missing data")
 if (is.vector(obj)) {
  if (i == 1 & j == 1 & length(obj) == n) out <- as.double(obj)
  else stop("Unexpected length of vector")
 } else if (is.matrix(obj)) {
  if (nrow(obj) == i && ncol(obj) == j) out <- as.double(rep(t(obj), n))
  else stop("Unexpected dimensions of matrix")
 } else {
  stop("Unexpected object")
 }	
 return(out) 
}

dnig <- function(x, alpha = 1, beta = 0, delta = 1, mu = 0, log = FALSE)
{   
 return( dgh(x = x, alpha = alpha, beta = beta, delta = delta, mu = mu, 
             lambda = -0.5, log = log) )
}
pnig <- function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lower_tail = TRUE, log = FALSE)
{   
 return( pgh(q = q, alpha = alpha, beta = beta, delta = delta, mu = mu, 
             lambda = -0.5, lower_tail = lower_tail, log = log) )
}
qnig <- function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lower_tail = TRUE, log = FALSE)
{
 value_len <- c(length(p), length(alpha), length(beta), length(delta), length(mu))
 max_n <- max(value_len)
 if (value_len[1] != max_n) p <- rep(p[1], max_n)
 if (value_len[2] != max_n) alpha <- rep(alpha[1], max_n)
 if (value_len[3] != max_n) beta <- rep(beta[1], max_n)
 if (value_len[4] != max_n) delta <- rep(delta[1], max_n)
 if (value_len[5] != max_n) mu <- rep(mu[1], max_n)
 ans <- rep(NA, max_n)
 for (i in 1:max_n) {
  ans[i] <- .qnigC(p = p[i], alpha = alpha[i], beta = beta[i], delta = delta[i], mu = mu[i], lower_tail = lower_tail, log = log)
 }
 return(ans)
}

rnig <- function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
 return(rgh(n, alpha = alpha, beta = beta, delta = delta, mu = mu, lambda = -0.5))
}
# ------------------------------------------------------------------------------
# Standardized Normal Inverse Gaussian (NIG) Distribution
# ------------------------------------------------------------------------------
.qsnigC <- function(p, rho = 0, zeta = 1, lower_tail = TRUE, log = FALSE) 
{
 param <- .paramGH(rho, zeta, lambda = -0.5)
 return( .qnigC(p, param[1], param[2], param[3], param[4], lower_tail = lower_tail, log = log) )
}

#' (Standardized) Normal Inverse Gaussian Distribution
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of the 
#' standardized normal inverse gaussian distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
#' and r generates random deviates. Output depends on x or q length, or n for the random number 
#' generator
#' @rdname snig
#' @export
#'
#'
#'
dsnig <- function(x, mu = 0, sigma = 1, skew = 0, shape = 1, log = FALSE) 
{
 value_len <- c(length(x), length(mu), length(sigma), length(skew), length(shape))
 max_n <- max(value_len)
 if (value_len[1] != max_n) x <- rep(x[1], max_n)
 if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
 if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
 if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
 if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
 ans <- double(max_n)
 sol <- try(.C("c_dsnig", x = as.double(x), mu = as.double(mu), 
              sigma = as.double(sigma), skew = as.double(skew), 
              shape = as.double(shape), ans = ans, n = as.integer(max_n), 
              logr = as.integer(log), PACKAGE = "tsdistributions"), silent = TRUE)
 if (inherits(sol, 'try-error')) {
  return(sol)
 } else {
  return(sol$ans)
 }
}

#' @rdname snig
#' @export
psnig <- function(q, mu = 0, sigma = 1, skew = 0, shape = 1, lower_tail = TRUE, log = FALSE) 
{
 return( psgh(q, mu, sigma, skew, shape, lambda = -0.5, lower_tail = lower_tail, log = log) )
}

#' @rdname snig
#' @export
qsnig = function(p, mu = 0, sigma = 1, skew = 0, shape = 1, lower_tail = TRUE, log = FALSE) 
{
 value_len <- c(length(p), length(mu), length(sigma), length(skew), length(shape))
 max_n <- max(value_len)
 if (value_len[1] != max_n) p <- rep(p[1], max_n)
 if (value_len[2] != max_n) mu <- rep(mu[1], max_n)
 if (value_len[3] != max_n) sigma <- rep(sigma[1], max_n)
 if (value_len[4] != max_n) skew <- rep(skew[1], max_n)
 if (value_len[5] != max_n) shape <- rep(shape[1], max_n)
 ans <- double(max_n)
 for (i in 1:max_n) {
  ans[i] = mu[i] + sigma[i]*.qsnigC(p[i], rho = skew[i], zeta = shape[i], lower_tail = lower_tail, log = log)
 }
 return(ans)
}

#' @rdname snig
#' @export
rsnig <- function(n, mu = 0, sigma = 1, skew = 0, shape = 1)
{
 value_len <- c(length(mu), length(sigma), length(skew), length(shape))
 if (value_len[1] != n) mu <- rep(mu[1], n)
 if (value_len[2] != n) sigma <- rep(sigma[1], n)
 if (value_len[3] != n) skew <- rep(skew[1], n)
 if (value_len[4] != n) shape <- rep(shape[1], n)
 ans <- double(n)
 sol <- try(.C("c_rsnig", n = as.integer(n), mu = as.double(mu), 
              sigma = as.double(sigma), skew = as.double(skew), 
              shape = as.double(shape),
              ans = ans, PACKAGE = "tsdistributions"), silent = TRUE)
 if (inherits(sol, 'try-error')) {
  return(sol)
 } else {
  return(sol$ans)
 }
}
#' Johnson's SU Distribution
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of Johnson's SU 
#' distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
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
qjsu = function(p, mu = 0, sigma = 1, skew = 1, shape = 0.5, lower_tail = TRUE, log = FALSE)
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
# ------------------------------------------------------------------------------
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

#---------------------------------------------------------------------------------
# functions for export:
#---------------------------------------------------------------------------------
#' Parameter Transformation
#' @description Transforms parameters from standardized representation to distribution
#' specific representation for the nig and ghyp distributions.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape shape parameter.
#' @param lambda additional shape parameter of the ghyp distribution.
#' @return The (alpha, beta, delta, mu) representation.
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
ghyptransform = function(mu = 0, sigma = 1,  skew = 0, shape = 3, lambda = -0.5)
{
 return(.ghypscale(mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda))
}

#' Distributions pqdr wrapper
#'
#' @description Density, distribution, quantile function and random number 
#' generation for the location scale invariant parameterization of all 
#' the distributions in the package.
#' @param distribution a valid distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean.
#' @param sigma standard deviation.
#' @param skew skew parameter.
#' @param shape  shape parameter.
#' @param lambda additional shape parameter determining subfamilies of the ghyp 
#' distribution.
#' @param log (logical) if TRUE, probabilities p are given as log(p).
#' @param lower_tail if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @return d gives the density, p gives the distribution function, q gives the quantile function 
#' and r generates random deviates. Output depends on x or q length, or n for the random number 
#' generator
#' @rdname ddist
#' @export
#' 
ddist <- function(distribution = "norm", x, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, log = FALSE)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 ans <- switch(distribution,
              norm = dnorm(x, mean = mu, sd = sigma, log = log),
              snorm = dsnorm(x, mu = mu, sigma = sigma, skew = skew, log = log),
              std = dstd(x, mu = mu, sigma = sigma, shape = shape, log = log),
              sstd = dsstd(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log),
              ged = dged(x, mu = mu, sigma = sigma, shape = shape, log = log),
              sged = dsged(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log),
              nig = dsnig((x - mu)/sigma, skew = skew, shape = shape, log = log)/sigma,
              ghyp = dsgh((x - mu)/sigma, skew = skew, shape = shape, lambda = lambda, log = log)/sigma,
              jsu = djsu(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log),
              ghst = dsghst(x, mu = mu, sigma = sigma, skew = skew, shape = shape, log = log)
 )
 return(ans)
}

#' @rdname ddist
#' @export
pdist <- function(distribution = "norm", q, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, lower_tail = TRUE, log = FALSE)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 ans <- switch(distribution,
              norm = pnorm(q, mean = mu, sd = sigma, lower.tail = lower_tail, log.p = log),
              snorm = psnorm(q, mu = mu, sigma = sigma, skew = skew, lower_tail = lower_tail, log = log),
              std = pstd(q, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log),
              sstd = psstd(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              ged = pged(q, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log),
              sged = psged(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              nig = psnig((q - mu)/sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              ghyp = psgh((q - mu)/sigma, skew = skew, shape = shape, lambda = lambda, lower_tail = lower_tail, log = log),
              jsu = pjsu(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              ghst = psghst(q, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log)
 )
 return(ans)
}

#' @rdname ddist
#' @export
qdist <- function(distribution = "norm", p, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5, lower_tail = TRUE, log = FALSE)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 ans <- switch(distribution,
              norm = qnorm(p, mean = mu, sd = sigma, lower.tail = lower_tail, log.p = log),
              snorm = qsnorm(p, mu = mu, sigma = sigma, skew = skew, lower_tail = lower_tail, log = log),
              std = qstd(p, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log),
              sstd = qsstd(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              ged = qged(p, mu = mu, sigma = sigma, shape = shape, lower_tail = lower_tail, log = log),
              sged = qsged(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              nig = qsnig(p, skew = skew, shape = shape, lower_tail = lower_tail, log = log)*sigma + mu,
              ghyp = qsgh(p, skew = skew, shape = shape, lambda = lambda, lower_tail = lower_tail, log = log)*sigma + mu,
              jsu = qjsu(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
              ghst = qsghst(p, mu = mu, sigma = sigma, skew = skew, shape = shape, lower_tail = lower_tail, log = log),
 )
 return(ans)
}

#' @rdname ddist
#' @export
rdist <- function(distribution = "norm", n, mu = 0, sigma = 1, skew = 1, shape = 5, lambda = -0.5)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 
 ans <- switch(distribution,
              norm = rnorm(n, mean = mu, sd = sigma),
              snorm = rsnorm(n, mu = mu, sigma = sigma, skew = skew),
              std = rstd(n, mu = mu, sigma = sigma, shape = shape),
              sstd = rsstd(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
              ged = rged(n, mu = mu, sigma = sigma, shape = shape),
              sged = rsged(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
              nig =  mu + sigma*rsnig(n, skew = skew, shape = shape),
              ghyp = mu + sigma*rsgh(n, skew = skew, shape = shape, lambda = lambda),
              jsu = rjsu(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
              ghst = rsghst(n, mu = mu, sigma = sigma, skew = skew, shape = shape)
 )
 return(ans)
}

#' Distribution skewness and kurtosis
#'
#' @description Calculates the skewness and kurtosis of the distribution given
#' a set of parameters.
#' @param distribution a valid distribution.
#' @param skew skew parameter.
#' @param shape  shape parameter.
#' @param lambda additional shape parameter determining subfamilies of the ghyp 
#' distribution.
#' @return A numeric value.
#' @rdname dskewness
#' @export
#'
dskewness <- function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 f <- Vectorize(.dskewness)
 ans <- f(distribution, skew, shape, lambda)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.dskewness <- function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
 ans <- switch(distribution,
              norm 	= 0,
              snorm 	= .snormskew(skew = skew),
              std 	= 0,
              sstd 	= .sstdskew(skew = skew, shape = shape),
              ged 	= 0,
              sged 	= .sgedskew(skew = skew, shape = shape),
              nig 	= .snigskew(skew = skew, shape = shape),
              ghyp 	= .sghypskew(skew = skew, shape = shape, lambda = lambda),
              jsu 	= .jsuskew(skew = skew, shape = shape),
              ghst	= .ghstskew(skew, shape)
 )
 return(as.numeric(ans))
}

#' @rdname dskewness
#' @export
dkurtosis <- function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
 distribution <- match.arg(distribution[1], valid_distributions())
 f <- Vectorize(.dkurtosis)
 ans <- f(distribution, skew, shape, lambda)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.dkurtosis <- function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
 ans <- switch(distribution,
              norm 	= 0,
              snorm 	= 0,
              std 	= .stdexkurt(shape = shape),
              sstd 	= .sstdexkurt(skew = skew, shape = shape),
              ged 	= .gedexkurt(shape = shape),
              sged 	= .sgedexkurt(skew = skew, shape = shape),
              nig 	= .snigexkurt(skew = skew, shape = shape),
              ghyp 	= .sghypexkurt(skew = skew, shape = shape, lambda = lambda),
              jsu 	= .jsuexkurt(skew = skew, shape = shape),
              ghst	= .ghstexkurt(skew, shape)
 )
 return(as.numeric(ans))
}


# NIG Moments
.nigmu <- function(alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 ans <- mu + (delta*beta)/gm
 return(ans)
}

.nigsigma <- function(alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 ans <- sqrt((delta*alpha^2)/(gm^3))
 return(ans)
}

.snigskew <- function(skew, shape){
 fun <- function(skew, shape){
  pars <- .paramGH(rho = skew, zeta = shape, lambda = -0.5)
  return(.nigskew(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
 }
 f <- Vectorize(fun)
 ans <- f(skew, shape)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.nigskew <- function(alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 ans <- 3*beta/(alpha*sqrt(delta*gm))
 return(ans)
}

.snigexkurt <- function(skew, shape){
 fun <- function(skew, shape){
  pars <- .paramGH(rho = skew, zeta = shape, lambda = -0.5)
  return(.nigexkurt(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
 }
 f <- Vectorize(fun)
 ans <- f(skew, shape)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.nigexkurt <- function(alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 ans <- 3*(1 + 4*(beta^2)/(alpha^2))/(delta*gm)
 return(ans)
}


.ghypmu <- function(lambda, alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 ans <- mu + (delta * beta * besselK( delta * gm, lambda + 1) )/( gm * besselK( delta * gm, lambda) ) 
 return(ans)
}

.ghypsigma <- function(lambda, alpha, beta, delta, mu){
 gm <- sqrt(alpha^2 - beta^2)
 x1 <- delta * besselK( delta * gm, lambda + 1) / ( gm * besselK( delta * gm, lambda) )
 x2 <- ( (beta^2 * delta^2) / (gm^2) ) * ( ( besselK( delta * gm, lambda + 2) / besselK( delta * gm, lambda) ) 
                                          - ( besselK( delta * gm, lambda + 1)^2 / besselK( delta * gm, lambda)^2 ) )
 ans <- sqrt(x1 + x2)
 return(ans)
}

.sghypskew <- function(skew, shape, lambda){
 fun <- function(skew, shape, lambda){
  pars <- .paramGH(rho = skew, zeta = shape, lambda = lambda)
  return(.ghypskew(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
 }
 f <- Vectorize(fun)
 ans <- f(skew, shape, lambda)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.ghypskew <- function(lambda, alpha, beta, delta, mu){
 skew <- ghypMom(3, lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu, momType = "central")/(.ghypsigma(lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu)^3)
 return(skew)
}

.sghypexkurt <- function(skew, shape, lambda){
 fun <- function(skew, shape, lambda){
  pars <- .paramGH(rho = skew, zeta = shape, lambda = lambda)
  return(.ghypexkurt(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
 }
 f <- Vectorize(fun)
 ans <- f(skew, shape, lambda)
 if (NCOL(ans) == 1) ans <- as.numeric(ans)
 return(ans)
}

.ghypexkurt <- function(lambda, alpha, beta, delta, mu){
 kurt <- ghypMom(4, lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu, momType = "central")/(.ghypsigma(lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu)^4) - 3
 return(kurt)
}

.norm2snorm1 <- function(mu, sigma, skew)
{
 m1 <- 2/sqrt(2 * pi)
 ans <- mu + m1 * (skew - 1/skew)*sigma
 return(ans)
}

.norm2snorm2 <- function(mu, sigma, skew)
{
 m1 <- 2/sqrt(2 * pi)
 m2 <- 1
 sigx <- sqrt((1 - m1^2)*(skew^2 + 1/skew^2) + 2*m1^2 - 1)
 ans <- sigx * sigma
 return(ans)
}

.snormskew <- function(skew)
{
 m1 <- 2/sqrt(2 * pi)
 m2 <- 1
 m3 <- 4/sqrt(2 * pi)
 ans <- (skew - 1/skew) * (( m3 + 2 * m1^3 - 3 * m1 * m2) * (skew^2 + (1/skew^2)) + 3 * m1 * m2 - 4 * m1^3 )/(((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2)^(3/2))
 return(ans)
}

.snormexkurt <- function(skew)
{
  return(0)
}

.stdskew <- function(shape)
{
 return(0)
}

.stdexkurt <- function(shape)
{
 ans <- ifelse(shape > 4, 6/(shape - 4), NA)
 return(ans)
}

.sstdskew <- function(skew, shape){
 # Theoretical moments based on bijection betweeen Fernandez and Steel versions
 # and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
 if (shape > 2) {
  eta <- shape
  k2 <- skew^2
  lda <- (k2 - 1)/(k2 + 1)
  ep1 <- (eta + 1)/2
  lnc <- lgamma(ep1) - lgamma(eta/2) - 0.5*log(pi*(eta - 2))
  c <- exp(lnc)
  a <- 4*lda*c*(eta - 2)/(eta - 1)
  b <- sqrt(1 + 3*lda^2 - a^2)
  my2 <- 1 + 3*lda^2
  my3 <- 16*c*lda*(1 + lda^2)*((eta - 2)^2)/((eta - 1)*(eta - 3))
  my4 <- 3*(eta - 2)*(1 + 10*lda^2 + 5*lda^4)/(eta - 4)
  m3 <- (my3 - 3*a*my2 + 2*a^3)/(b^3)
 } else {
  m3 <- NA
 }
 return(m3)
}

.sstdexkurt <- function(skew, shape)
{
 # Theoretical moments based on bijection betweeen Fernandez and Steel versions
 # and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
 if (shape > 4) {
  eta <- shape
  k2 <- skew^2
  lda <- (k2 - 1)/(k2 + 1)
  ep1 <- (eta + 1)/2
  lnc <- lgamma(ep1) - lgamma(eta/2) - 0.5*log(pi*(eta - 2))
  c <- exp(lnc)
  a <- 4*lda*c*(eta - 2)/(eta - 1)
  b <- sqrt(1 + 3*lda^2 - a^2)
  my2 <- 1 + 3*lda^2
  my3 <- 16*c*lda*(1 + lda^2)*((eta - 2)^2)/((eta - 1)*(eta - 3))
  my4 <- 3*(eta - 2)*(1 + 10*lda^2 + 5*lda^4)/(eta - 4)
  m3 <- (my3 - 3*a*my2 + 2*a^3)/(b^3)
  m4 <- -3 + (my4 - 4*a*my3 + 6*(a^2)*my2 - 3*a^4)/(b^4)
 } else {
  m4 <- NA
 }
 return(m4)
}

.gedskew <- function(shape)
{
 return(0)
}

.gedexkurt <- function(shape)
{
  ans <- (((gamma(1/shape)/gamma(3/shape))^2) * (gamma(5/shape)/gamma(1/shape))) - 3
  return(ans)
}

.sgedskew <- function(skew, shape)
{
 lambda <- sqrt(2^(-2/shape)*gamma(1/shape)/gamma(3/shape))
 m1 <- ((2^(1/shape)*lambda)^1*gamma(2/shape)/gamma(1/shape))
 m2 <- 1
 m3 <- ((2^(1/shape)*lambda)^3*gamma(4/shape)/gamma(1/shape))
 ans <- (skew - 1/skew) * ((m3 + 2 * m1^3 - 3 * m1 * m2) * (skew^2 + (1/skew^2)) + 3 * m1 * m2 - 4 * m1^3)/(((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2))
 return(ans) 
}

.sgedexkurt <- function(skew, shape)
{
 lambda <- sqrt(2^(-2/shape) * gamma(1/shape)/gamma(3/shape))
 m1 <- ((2^(1/shape)*lambda)^1 * gamma(2/shape)/gamma(1/shape))
 m2 <- 1
 m3 <- ((2^(1/shape)*lambda)^3 * gamma(4/shape)/gamma(1/shape))
 m4 <- ((2^(1/shape)*lambda)^4 * gamma(5/shape)/gamma(1/shape))
 cm4 <- (-3 * m1^4 * (skew - 1/skew)^4) + 
  (6 * m1^2 * (skew - 1/skew)^2 * m2*(skew^3 + 1/skew^3) )/(skew + 1/skew) - 
  (4 * m1*(skew - 1/skew) * m3 * (skew^4 - 1/skew^4))/(skew + 1/skew) + 
  (m4 * (skew^5 + 1/skew^5))/(skew + 1/skew)
 ans <- (cm4/(((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ 2)) - 3
 return(ans)
}

.jsuskew <- function(mu = 0, sigma = 1, skew, shape)
{	
 Omega <- -skew/shape
 w <- exp(shape^-2)
 s3 <- -0.25*sqrt(w)*((w - 1)^2)*(w*(w + 2)*sinh(3*Omega) + 3*sinh(Omega))
 ans <- s3/(0.5*(w - 1)*(w*cosh(2*Omega) + 1))^(3/2)
 return(ans)
}

.jsuexkurt <- function(mu = 0, sigma = 1, skew, shape)
{
 Omega <- -skew/shape
 w <- exp(shape^-2)
 s4 <- 0.125 * (w - 1)^2*(w^2*(w^4 + 2*w^3 + 3*w^2 - 3)*cosh(4*Omega) + 4*w^2*(w + 2)*cosh(2*Omega) + 3*(2*w + 1))
 ans <- s4/(0.5*(w - 1)*(w*cosh(2*Omega) + 1))^2
 return(ans - 3)
}

.ghstskew <- function(skew, shape){
 if (shape < 6) {
  ans <- NA
 } else {
  params <- .paramGHST(nu = shape, betabar = skew)
  delta <- params[2]
  beta <- params[3]
  nu <- params[4]
  beta2 <- beta*beta
  delta2 <- delta*delta
  ans <- ((2 * sqrt(nu - 4)*beta*delta)/((2*beta2*delta2 + (nu - 2)*(nu - 4))^(3/2))) * (3*(nu - 2) + ((8*beta2*delta2)/(nu - 6)))
 }
 return(ans)
}

.ghstexkurt <- function(skew, shape){
 if (shape < 8) {
  ans <- NA
 } else {
  params <- .paramGHST(nu = shape, betabar = skew)
  delta <- params[2]
  beta <- params[3]
  nu <- params[4]
  beta2 <- beta*beta
  delta2 <- delta*delta
  k1 <- 6/((2*beta2*delta2 + (nu - 2)*(nu - 4))^2)
  k21 <- (nu - 2)*(nu - 2)*(nu - 4)
  k22 <- (16*beta2*delta2*(nu - 2)*(nu - 4))/(nu - 6)
  k23 <- (8*(beta2^2)*(delta2^2)*(5*nu - 22))/((nu - 6)*(nu - 8))
  ans <- k1*(k21 + k22 + k23)
 }
 return( ans )
}

validate_bounds <- function(distribution, sigma = 1, skew = 0.2, shape = 5, lambda = 1, return_table = FALSE)
{
  distribution <- match.arg(distribution[1], valid_distributions())
  value <- lower <- upper <- lower_check <- upper_check <- NULL
  p <- distribution_bounds(distribution)
  p <- p[-which(p$parameter == "mu")]
  newd <- data.table(parameter = c("sigma","skew","shape","lambda"), value = c(sigma, skew, shape, lambda))
  out <- merge(p, newd, by = "parameter", all.x = T)
  out[,lower_check := value >= lower]
  out[,upper_check := value <= upper]
  if (return_table) {
    return(out)
  } else {
    if (any(!out$lower_check) | any(!out$upper_check)) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

#' Distribution Authorized Domain
#'
#' @description Calculated the region of Skewness-Kurtosis for which a density exists.
#' @param distribution a valid distribution with skew and shape parameters.
#' @param max_kurt the maximum kurtosis for which to determine the bounds for 
#' the skewness-kurtosis domain.
#' @param n  the number of points between the lower and upper bounds of the 
#' skew and shape parameters for which to evaluate the skewness and excess kurtosis. 
#' This determines the kurtosis interval (3 - max_kurt) for which to 
#' calculate (solver based) the maximum skewness.
#' @param lambda additional shape parameter determining subfamilies of the ghyp 
#' distribution.
#' @return A list with the lower half of the skewness and kurtosis values.
#' @rdname authorized_domain
#' @export
#'
authorized_domain <- function(distribution, max_kurt = 30, n = 25, lambda = 1) {
  parameter <- NULL
  valid_d <- c("nig","ghyp","jsu","sstd","ghst")
  distribution <- match.arg(distribution[1], valid_d)
  di <- skdomain_bounds(distribution)
  k <- seq(5, max_kurt, length = n)
  dpars <- c(di[parameter == "skew"]$value, di[parameter == "shape"]$value)
  # sstd skew = 1 is equivalent to zero skewness
  skew_min <- switch(distribution, "jsu" = 0, "nig" = 0, "sstd" = 1, "ghyp" = 0, "ghst" = 0)
  if (distribution == "ghyp") {
    maxkurt <- dkurtosis(distribution, skew = skew_min, shape = di[parameter == "shape"]$upper, lambda = lambda)
    f1 <- function(x, kurt, xlambda){
      -dskewness(distribution, skew = x[1], shape = x[2], lambda = xlambda)
    }
    fin1 <- function(x, kurt, xlambda){
      dkurtosis(distribution, skew = x[1], shape = x[2], lambda = xlambda) + 3 - maxkurt - kurt
    }
    pars <- matrix(NA, ncol = 4, nrow = n)
    for (i in 1:length(k)) {
      sol <- try(solnp(pars = dpars, fun = f1, eqfun = fin1, eqB = 0, LB = c(di[parameter == "skew"]$lower, di[parameter == "shape"]$lower),
                   UB = c(di[parameter == "skew"]$upper, di[parameter == "shape"]$upper), 
                   control = list(trace = 0, outer.iter = 25), kurt = k[i], xlambda = lambda), silent = TRUE)
      if (inherits(sol, 'try-error')) {
        pars[i,1:2] <- rep(NA, 2)
        pars[i,3] <- NA
        pars[i,4] <- NA
      } else {
        if (any(is.na(sol$pars))) {
          pars[i,1:2] <- rep(NA, 2)
          pars[i,3] <- NA
          pars[i,4] <- NA
        } else {
          pars[i,1:2] <- sol$pars
          pars[i,3] <- tail(sol$value,1)
          pars[i,4] <- dkurtosis(distribution, skew =  sol$pars[1], shape =  sol$pars[2], lambda = lambda) + 3
        }
      }
    }
    pars <- rbind(matrix(c(0, di[parameter == "shape"]$upper, dskewness(distribution, 0,  di[parameter == "shape"]$upper, lambda = lambda), 
                           3 + dkurtosis(distribution, 0,  di[parameter == "shape"]$upper, lambda = lambda)), ncol = 4), pars)
  } else {
    maxkurt <- dkurtosis(distribution, skew = skew_min, shape = di[parameter == "shape"]$upper)
    f2 <- function(x, kurt){
      -dskewness(distribution, skew = x[1], shape = x[2])
    }
    fin2 <- function(x, kurt){
      dkurtosis(distribution, skew = x[1], shape = x[2]) + 3 - maxkurt - kurt
    }
    pars <- matrix(NA, ncol = 4, nrow = n)
    for (i in 1:length(k)) {
      sol <- try(solnp(pars = dpars, fun = f2, eqfun = fin2, eqB = 0, LB = c(di[parameter == "skew"]$lower, di[parameter == "shape"]$lower),
                   UB = c(di[parameter == "skew"]$upper, di[parameter == "shape"]$upper), 
                   control = list(trace = 0, outer.iter = 25), kurt = k[i]), silent = TRUE)
      if (inherits(sol, 'try-error')) {
        pars[i,1:2] <- rep(NA, 2)
        pars[i,3] <- NA
        pars[i,4] <- NA
      } else {
        if (any(is.na(sol$pars))) {
          pars[i,1:2] <- rep(NA, 2)
          pars[i,3] <- NA
          pars[i,4] <- NA
        } else {
          pars[i,1:2] <- sol$pars
          pars[i,3] <- tail(sol$value,1)
          pars[i,4] <- dkurtosis(distribution, skew =  sol$pars[1], shape =  sol$pars[2]) + 3
        }
      }
    }
    pars <- rbind(matrix(c(0, di[parameter == "shape"]$upper, dskewness(distribution, skew_min,  di[parameter == "shape"]$upper), 3 + dkurtosis(distribution, skew_min,  di[parameter == "shape"]$upper)), ncol = 4), pars)
  }
  ans <- spline(pars[,4], pars[,3], method = "fmm")
  return(list(Skewness = ans$y, Kurtosis = ans$x, pars = pars[,1:2]))
}

# valid bounds for the existence of a fourth moment
skdomain_bounds <- function(distribution)
{
  parameter <- NULL
  d <- distribution_bounds(distribution)
  if (distribution == "sstd") {
    out <- rbind(data.table(parameter = "skew", lower = 1, upper = 60, value = 1.1),
          data.table(parameter = "shape", lower = 4.01, upper = 300, value = 4.05))
  } else if (distribution == "nig") {
    out <- rbind(data.table(parameter = "skew", lower = 0.05, upper = d[parameter == "skew"]$upper, value = 0.1),
                 data.table(parameter = "shape", lower = d[parameter == "shape"]$lower, upper = d[parameter == "shape"]$upper, value = 0.5))
  } else if (distribution == "ghst") {
    out <- rbind(data.table(parameter = "skew", lower = 0.1, upper = d[parameter == "skew"]$upper, value = 0.2),
                 data.table(parameter = "shape", lower = 8.01, upper = d[parameter == "shape"]$upper, value = 8.1))
  } else if (distribution == "ghyp") {
    out <- rbind(data.table(parameter = "skew", lower = 0.05, upper = d[parameter == "skew"]$upper, value = 0.1),
                 data.table(parameter = "shape", lower = d[parameter == "shape"]$lower, upper = d[parameter == "shape"]$upper, value = 0.5))
  } else if (distribution == "jsu") {
    out <- rbind(data.table(parameter = "skew", lower = 0.05, upper = d[parameter == "skew"]$upper, value = 0.1),
                 data.table(parameter = "shape", lower = 0.1, upper = d[parameter == "shape"]$upper, value = 0.5))
  }
  return(out)
}