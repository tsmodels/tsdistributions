
bdiag <- function(...) 
{
    if (nargs() == 1) x <- as.list(...)
    else x <- list(...)
    n <- length(x)
    if (n == 0) return(NULL)
    x <- lapply(x, function(y) if (length(y)) 
        as.matrix(y)
        else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2], 
                                                 (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}


.interp_pdf <- function(mid_values, kernel_x, kernel_z)
{
    # m = middle values of dspd(x,...)
    # kernel_x = kernel interpolated density of original data
    # kernel_z = empirical cdf probability [0,1] (cumsum of kernel_y/sum of kernel_y)
    m <- length(mid_values)
    bins <- vector(mode = "numeric", length = m)
    nbin <- length(kernel_x)
    nx <- vector(mode = "numeric", length = nbin)
    z <- as.data.frame(1:nbin)
    kk <- apply(z, 1, FUN = function(x) which(mid_values >= kernel_x[x]))
    for (i in 1:nbin) {
        bins[kk[[i]]] <- i
    }
    nx <- apply(z, 1, FUN = function(x) length(kk[[x]]))
    kk <- which(mid_values > kernel_x[nbin])
    bins[kk] <- 0
    nx[nbin + 1] <- length(kk)
    counts <- -diff(nx)
    bin <- bins
    bin[mid_values == kernel_x[1]] <- 1
    bin[mid_values == kernel_x[length(kernel_x)]] <- length(counts) - 1
    fx <- vector(mode = "numeric", length = length(mid_values))
    tx <- bin > 0
    bin <- bin[tx]
    # find the PDF probability
    tmp <- (kernel_z[bin + 1] - kernel_z[bin]) / (kernel_x[bin + 1] - kernel_x[bin])
    tna <- which(!is.na(tmp))
    tx <- tx[tna]
    fx[tna] <- tmp[!is.na(tmp)]
    return(fx)
}

.find_threshold <- function(x, exceed)
{
    x <- rev(sort(x))
    uniq <- unique(x)
    idx <- match(x[exceed], uniq)
    idx <- pmin(idx + 1, length(uniq))
    return(uniq[idx])
}
