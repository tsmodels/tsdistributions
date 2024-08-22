
#' Distribution skewness and kurtosis
#'
#' @description Calculates the skewness and excess kurtosis of the distribution given
#' a set of parameters.
#' @param distribution a valid distribution.
#' @param skew skew parameter.
#' @param shape  shape parameter.
#' @param lambda additional shape parameter for the Generalized Hyperbolic
#' distribution.
#' @returns A numeric value for the skewness and excess kurtosis.
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
                  "norm" = 0, 
                  "snorm" = .snormskew(skew = skew), 
                  "std" = 0,
                  "sstd" = .sstdskew(skew = skew, shape = shape), 
                  "ged" = 0, 
                  "sged" = .sgedskew(skew = skew, shape = shape), 
                  "nig" = .snigskew(skew = skew, shape = shape),
                  "gh" = .sghypskew(skew = skew, shape = shape, lambda = lambda), 
                  "jsu" = .jsuskew(skew = skew, shape = shape), 
                  "ghst" = .ghstskew(skew, shape)
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
                  "norm" = 0, 
                  "snorm" = 0, 
                  "std" = .stdexkurt(shape = shape), 
                  "sstd" = .sstdexkurt(skew = skew, shape = shape), 
                  "ged" = .gedexkurt(shape = shape), 
                  "sged" = .sgedexkurt(skew = skew, shape = shape), 
                  "nig" = .snigexkurt(skew = skew, shape = shape), 
                  "gh" = .sghypexkurt(skew = skew, shape = shape, lambda = lambda), 
                  "jsu" = .jsuexkurt(skew = skew, shape = shape), 
                  "ghst" = .ghstexkurt(skew, shape)
    )
    return(as.numeric(ans))
}


# NIG Moments
.nigmu <- function(alpha, beta, delta, mu){
    gm <- sqrt(alpha^2 - beta^2)
    ans <- mu + (delta * beta)/gm
    return(ans)
}

.nigsigma <- function(alpha, beta, delta, mu){
    gm <- sqrt(alpha^2 - beta^2)
    ans <- sqrt((delta * alpha^2)/(gm^3))
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
    ans <- 3*beta / (alpha * sqrt(delta*gm))
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
    ans <- 3 * (1 + 4 * (beta^2) / (alpha^2)) / (delta * gm)
    return(ans)
}


.ghypmu <- function(lambda, alpha, beta, delta, mu){
    gm <- sqrt(alpha^2 - beta^2)
    ans <- mu + (delta * beta * besselK(delta * gm, lambda + 1)) / (gm * besselK(delta * gm, lambda))
    return(ans)
}

.ghypsigma <- function(lambda, alpha, beta, delta, mu){
    gm <- sqrt(alpha^2 - beta^2)
    x1 <- delta * besselK(delta * gm, lambda + 1) / (gm * besselK(delta * gm, lambda))
    x2 <- ((beta^2 * delta^2) / (gm^2)) * ((besselK(delta * gm, lambda + 2) / besselK(delta * gm, lambda)) 
                                           - (besselK(delta * gm, lambda + 1)^2 / besselK(delta * gm, lambda)^2))
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
    skew <- ghypMom(3, lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu, momType = "central") / 
        (.ghypsigma(lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu)^3)
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
    kurt <- ghypMom(4, lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu, momType = "central") / 
        (.ghypsigma(lambda = lambda, alpha = alpha, beta = beta, delta = delta, mu = mu)^4) - 3
    return(kurt)
}

.norm2snorm1 <- function(mu, sigma, skew)
{
    m1 <- 2 / sqrt(2 * pi)
    ans <- mu + m1 * (skew - 1 / skew) * sigma
    return(ans)
}

.norm2snorm2 <- function(mu, sigma, skew)
{
    m1 <- 2 / sqrt(2 * pi)
    m2 <- 1
    sigx <- sqrt((1 - m1^2) * (skew^2 + 1 / skew^2) + 2 * m1^2 - 1)
    ans <- sigx * sigma
    return(ans)
}

.snormskew <- function(skew)
{
    m1 <- 2 / sqrt(2 * pi)
    m2 <- 1
    m3 <- 4 / sqrt(2 * pi)
    ans <- (skew - 1/skew) * ((m3 + 2 * m1^3 - 3 * m1 * m2) * (skew^2 + (1/skew^2)) + 3 * m1 * m2 - 4 * m1^3) / 
        (((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2)^(3/2))
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
    ans <- ifelse(shape > 4, 6 / (shape - 4), NA)
    return(ans)
}

.sstdskew <- function(skew, shape){
    # Theoretical moments based on bijection betweeen Fernandez and Steel versions
    # and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
    if (shape > 2) {
        eta <- shape
        k2 <- skew^2
        lda <- (k2 - 1) / (k2 + 1)
        ep1 <- (eta + 1) / 2
        lnc <- lgamma(ep1) - lgamma(eta/2) - 0.5*log(pi * (eta - 2))
        cx <- exp(lnc)
        a <- 4 * lda * cx * (eta - 2) / (eta - 1)
        b <- sqrt(1 + 3 * lda^2 - a^2)
        my2 <- 1 + 3 * lda^2
        my3 <- 16 * cx * lda * (1 + lda^2) * ((eta - 2)^2) / ((eta - 1) * (eta - 3))
        my4 <- 3 * (eta - 2) * (1 + 10 * lda^2 + 5 * lda^4) / (eta - 4)
        m3 <- (my3 - 3 * a * my2 + 2 * a^3) / (b^3)
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
        lda <- (k2 - 1) / (k2 + 1)
        ep1 <- (eta + 1) / 2
        lnc <- lgamma(ep1) - lgamma(eta/2) - 0.5 * log(pi * (eta - 2))
        cx <- exp(lnc)
        a <- 4 * lda * cx * (eta - 2) / (eta - 1)
        b <- sqrt(1 + 3 * lda^2 - a^2)
        my2 <- 1 + 3 * lda^2
        my3 <- 16 * cx * lda * (1 + lda^2) * ((eta - 2)^2) / ((eta - 1) * (eta - 3))
        my4 <- 3 * (eta - 2) * (1 + 10 * lda^2 + 5 * lda^4) / (eta - 4)
        m4 <- -3 + (my4 - 4 * a * my3 + 6 * (a^2) * my2 - 3 * a^4) / (b^4)
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
    ans <- (((gamma(1/shape) / gamma(3/shape))^2) * (gamma(5/shape) / gamma(1/shape))) - 3
    return(ans)
}

.sgedskew <- function(skew, shape)
{
    lambda <- sqrt(2^(-2/shape) * gamma(1/shape) / gamma(3/shape))
    m1 <- ((2^(1/shape) * lambda)^1 * gamma(2/shape) / gamma(1/shape))
    m2 <- 1
    m3 <- ((2^(1/shape) * lambda)^3 * gamma(4/shape) / gamma(1/shape))
    ans <- (skew - 1/skew) * ((m3 + 2 * m1^3 - 3 * m1 * m2) * (skew^2 + (1/skew^2)) + 3 * m1 * m2 - 4 * m1^3) / 
        (((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2)^(3/2))
    return(ans)
}

.sgedexkurt <- function(skew, shape)
{
    lambda <- sqrt(2^(-2/shape) * gamma(1/shape) / gamma(3/shape))
    m1 <- ((2^(1/shape) * lambda)^1 * gamma(2/shape) / gamma(1/shape))
    m2 <- 1
    m3 <- ((2^(1/shape) * lambda)^3 * gamma(4/shape) / gamma(1/shape))
    m4 <- ((2^(1/shape) * lambda)^4 * gamma(5/shape) / gamma(1/shape))
    cm4 <- (-3 * m1^4 * (skew - 1/skew)^4) +
        (6 * m1^2 * (skew - 1/skew)^2 * m2 * (skew^3 + 1/skew^3)) / (skew + 1/skew) -
        (4 * m1 * (skew - 1/skew) * m3 * (skew^4 - 1/skew^4)) / (skew + 1/skew) +
        (m4 * (skew^5 + 1/skew^5)) / (skew + 1/skew)
    ans <- (cm4 / (((m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ 2)) - 3
    return(ans)
}

.jsuskew <- function(mu = 0, sigma = 1, skew, shape)
{
    omega <- -skew / shape
    w <- exp(shape^-2)
    s3 <- -0.25 * sqrt(w) * ((w - 1)^2) * (w * (w + 2) * sinh(3 * omega) + 3 * sinh(omega))
    ans <- s3 / (0.5 * (w - 1) * (w * cosh(2*omega) + 1))^(3/2)
    return(ans)
}

.jsuexkurt <- function(mu = 0, sigma = 1, skew, shape)
{
    omega <- -skew / shape
    w <- exp(shape^-2)
    s4 <- 0.125 * (w - 1)^2 * (w^2 * (w^4 + 2 * w^3 + 3 * w^2 - 3) * cosh(4 * omega) + 
                                   4 * w^2 * (w + 2) * cosh(2 * omega) + 3 * (2 * w + 1))
    ans <- s4 / (0.5 * (w - 1) * (w * cosh(2 * omega) + 1))^2
    return(ans - 3)
}

.ghstskew <- function(skew, shape)
{
    if (shape < 6) {
        ans <- NA
    } else {
        params <- .paramGHST(nu = shape, betabar = skew)
        delta <- params[2]
        beta <- params[3]
        nu <- params[4]
        beta2 <- beta * beta
        delta2 <- delta * delta
        ans <- ((2 * sqrt(nu - 4) * beta * delta) / ((2 * beta2 * delta2 + (nu - 2) * (nu - 4))^(3/2))) * 
            (3 * (nu - 2) + ((8 * beta2 * delta2)/(nu - 6)))
    }
    return(ans)
}

.ghstexkurt <- function(skew, shape)
{
    if (shape < 8) {
        ans <- NA
    } else {
        params <- .paramGHST(nu = shape, betabar = skew)
        delta <- params[2]
        beta <- params[3]
        nu <- params[4]
        beta2 <- beta * beta
        delta2 <- delta * delta
        k1 <- 6 / ((2 * beta2 * delta2 + (nu - 2) * (nu - 4))^2)
        k21 <- (nu - 2) * (nu - 2) * (nu - 4)
        k22 <- (16 * beta2 * delta2 * (nu - 2) * (nu - 4)) / (nu - 6)
        k23 <- (8 * (beta2^2) * (delta2^2) * (5 * nu - 22)) / ((nu - 6) * (nu - 8))
        ans <- k1 * (k21 + k22 + k23)
    }
    return(ans)
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
#' @param lambda additional shape parameter for the Generalized Hyperbolic
#' distribution.
#' @returns A list with the lower half of the skewness and kurtosis values.
#' @rdname authorized_domain
#' @export
#'
authorized_domain <- function(distribution, max_kurt = 30, n = 25, lambda = 1) 
{
    parameter <- NULL
    valid_d <- c("nig","gh","jsu","sstd","ghst")
    distribution <- match.arg(distribution[1], valid_d)
    di <- skdomain_bounds(distribution)
    k <- seq(5, max_kurt, length = n)
    dpars <- c(di[parameter == "skew"]$value, di[parameter == "shape"]$value)
    # sstd skew = 1 is equivalent to zero skewness
    skew_min <- switch(distribution, "jsu" = 0, "nig" = 0, "sstd" = 1, "gh" = 0, "ghst" = 0)
    if (distribution == "gh") {
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
        pars <- rbind(matrix(c(0, di[parameter == "shape"]$upper, 
                               dskewness(distribution, skew_min,  di[parameter == "shape"]$upper), 
                               3 + dkurtosis(distribution, skew_min,  di[parameter == "shape"]$upper)), ncol = 4), 
                      pars)
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
    } else if (distribution == "gh") {
        out <- rbind(data.table(parameter = "skew", lower = 0.05, upper = d[parameter == "skew"]$upper, value = 0.1), 
                     data.table(parameter = "shape", lower = d[parameter == "shape"]$lower, upper = d[parameter == "shape"]$upper, value = 0.5))
    } else if (distribution == "jsu") {
        out <- rbind(data.table(parameter = "skew", lower = 0.05, upper = d[parameter == "skew"]$upper, value = 0.1), 
                     data.table(parameter = "shape", lower = 0.1, upper = d[parameter == "shape"]$upper, value = 0.5))
    }
    return(out)
}
