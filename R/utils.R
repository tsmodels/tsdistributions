.make_standard_errors <- function(pmatrix, H)
{
    pars <- pmatrix[estimate == 1]$value
    se <- sqrt(diag(solve(H)))
    tvalues <- pars/se
    pvalues <- 2*(1 - pnorm(abs(tvalues)))
    return(data.frame("Std. Error" = se,"t value" = tvalues, "Pr(>|t|)" = pvalues, check.names = FALSE))
}
