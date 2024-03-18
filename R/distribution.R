#' Distribution Bounds
#'
#' @details Returns the upper a lower bounds for the parameters of a distribution.
#' @param distribution A valid distribution
#' @returns A data.table of the parameters and their default bounds.
#' @export
#'
#'
#'
distribution_bounds <- function(distribution = "norm")
{
    tmp <- rbind(
        data.table(parameter = "mu", lower = -Inf, upper = Inf),
        data.table(parameter = "sigma", lower = 1e-12, upper = Inf))
    if (distribution == "norm") {
        tmp[,distribution := "norm"]
        return(tmp)
    }
    if (distribution == "ged") {
        tmp <- rbind(tmp, data.table(parameter = "shape", lower = 0.1, upper = 100))
        tmp[,distribution := "ged"]
        return(tmp)
    }
    if (distribution == "std") {
        tmp <- rbind(tmp,
                     data.table(parameter = "shape", lower = 2.01, upper = 100))
        tmp[,distribution := "std"]
        return(tmp)
    }
    if (distribution == "snorm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = 0.1, upper = 10))
        tmp[,distribution := "snorm"]
        return(tmp)
    }
    if (distribution == "sged") {
        tmp <- rbind(
            data.table(parameter = "skew", lower = 0.01, upper = 30),
            data.table(parameter = "shape", lower = 0.1, upper = 100))
        tmp[,distribution := "sged"]
        return(tmp)
    }
    if (distribution == "sstd") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = 0.01, upper = 30),
                     data.table(parameter = "shape", lower = 2.01, upper = 100))
        tmp[,distribution := "sstd"]
        return(tmp)
    }
    if (distribution == "nig") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = -0.99, upper = 0.99),
                     data.table(parameter = "shape", lower = 0.01, upper = 100))
        tmp[,distribution := "nig"]
        return(tmp)
    }
    if (distribution == "gh") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = -0.99, upper = 0.99),
                     data.table(parameter = "shape", lower = 0.25, upper = 100),
                     data.table(parameter = "lambda", lower = -30, upper = 30))
        tmp[,distribution := "gh"]
        return(tmp)
    }
    if (distribution == "jsu") {
        # johnson has 2 shape parameters. The second one we model with the "skew"
        # representation in rugarch
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = -20, upper = 20),
                     data.table(parameter = "shape", lower = 0.1, upper = 100))
        tmp[,distribution := "jsu"]
        return(tmp)
    }
    if (distribution == "ghst") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", lower = -80, upper = 80),
                     data.table(parameter = "shape", lower = 4.01, upper = 100))
        tmp[,distribution := "ghst"]
        return(tmp)
    }
}



distribution_class <- function(distribution)
{
    switch(distribution,
           "norm" = 1,
           "std" = 2,
           "snorm" = 3,
           "sstd" = 4,
           "ged" = 5,
           "sged" = 6,
           "nig" = 7,
           "gh" = 8,
           "jsu" = 9,
           "ghst" = 10
           )
}

valid_distributions <- function()
{
    c("norm", "std", "snorm", "sstd", "ged", "sged", "nig", "gh", "jsu", "ghst")
}

