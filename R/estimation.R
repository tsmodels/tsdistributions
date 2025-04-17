#' Estimates the parameters of a distribution using autodiff.
#'
#' @param object an object of class \dQuote{tsdistribution.spec}.
#' @param solver only \dQuote{nlminb} currently supported.
#' @param control solver control parameters.
#' @param use_hessian whether to use the hessian in the calculation.
#' @param ... additional parameters passed to the estimation function
#' @returns An object of class \dQuote{tsdistribution.estimate} with slots for 
#' the estimated coefficients, gradients, scores etc.
#' @details The estimation makes use of the TMB package for minimizing
#' the negative of the log-likelihood using automatic differentiation.
#' @method estimate tsdistribution.spec
#' @aliases estimate
#' @rdname estimate.tsdistribution.spec
#' @export
#'
#'
estimate.tsdistribution.spec <- function(object, solver = "nlminb", control = list(trace = 0, eval.max = 300, iter.max = 500), use_hessian = TRUE, ...)
{
    
    spec_list <- tsdistribution_tmb(object, use_hessian = use_hessian)
    other_opts <- list(...)
    if (!is.null(other_opts$silent)) {
        silent <- other_opts$silent
    } else {
        silent <- TRUE
    }
    fun <- try(MakeADFun(data = spec_list$data, parameters = spec_list$par_list, DLL = "tsdistributions_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent), silent = FALSE)
    if (inherits(fun, 'try-error')) {
        stop("\nestimate.tsdistribuition.modelspec found an error. Please use non ad version of estimator and contact developer with reproducible example.")
    }
    # currently we only support nlminb and optim (L-BFGS-B)
    solver <- match.arg(solver[1],c("nlminb","optim"))
    tsdistenv <- new.env()
    tsdistenv$llh <- 1
    tsdistenv$grad <- NULL
    tsdistenv$parameter_names <- names(fun$par)
    tsdistenv$estimation_names <- object$parmatrix[estimate == 1]$parameter
    tsdistenv$parmatrix <- object$parmatrix
    if (solver == "nlminb") {
        sol <- nlminb(start = fun$par, objective = spec_list$llh_fun, 
                      gradient = spec_list$grad_fun, hessian = spec_list$hess_fun, 
                      lower = spec_list$lower, upper = spec_list$upper,  control = control, 
                      fun = fun, tsdistenv = tsdistenv)
        pars <- sol$par
    } else if (solver == "optim") {
        sol <- optim(par = fun$par, fn = spec_list$llh_fun, gr = spec_list$grad_fun, 
                     lower = spec_list$lower, upper = spec_list$upper, 
                     method = "L-BFGS-B", fun = fun, tsdistenv = tsdistenv)
        pars <- sol$par
    }
    parmatrix <- object$parmatrix
    parmatrix[estimate == 1]$value <- pars
    par_list <- as.list(parmatrix$value)
    names(par_list) <- parmatrix$parameter
    object$parmatrix <- parmatrix
    loglik <- fun$fn(pars)
    gradient <- fun$gr(pars)
    hessian <- fun$he(pars, atomic = fun$env$usingAtomics())
    names(pars) <- tsdistenv$estimation_names
    colnames(gradient) <- tsdistenv$estimation_names
    colnames(hessian) <- rownames(hessian) <- tsdistenv$estimation_names
    scores <- score_function(object, pars, use_hessian)
    nobs <- NROW(object$target$y)
    df <- length(pars)
    out <- list(pars = pars, parmatrix = parmatrix, loglik = loglik, gradient = gradient, hessian = hessian, scores = scores, nobs = nobs, df = df, solver_out = sol, spec = object)
    class(out) <- "tsdistribution.estimate"
    return(out)
}

# convert specification to object ready for estimation using TMB
tsdistribution_tmb <- function(object, use_hessian = FALSE, silent = TRUE)
{
    include <- value <- parameter <- NULL
    setup <- copy(object$parmatrix)
    data_list <- list(y = object$target$y, dclass = distribution_class(object$distribution), model = "distribution")
    setup[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(setup[,list(P = 1:.N * include), by = "parameter"], by = "parameter", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(setup[,list(value, parameter)], by = "parameter", keep.by = FALSE), function(x) as.numeric(x$value))
    par_list <- parameters
    names(par_list) <- setup$parameter
    # create the functions
    parameter_names <- setup[estimate == 1]$parameter
    llh_fun <- function(pars, fun, tsdistenv)
    {
        names(pars) <- tsdistenv$parameter_names
        llh <- tsdistenv$nll
        llh <- try(fun$fn(pars), silent = TRUE)
        if (inherits(llh, 'try-error')) {
            llh <- llh + 0.2 * abs(llh)
        } else if (is.na(llh) | !is.finite(llh)) {
            llh <- abs(tsdistenv$nll) + 0.2 * abs(tsdistenv$nll)
        }
        tsdistenv$nll <- llh
        return(llh)
    }
 
    grad_fun <- function(pars, fun, tsdistenv)
    {
        names(pars) <- tsdistenv$parameter_names
        grad <- try(fun$gr(pars), silent = TRUE)
        if (inherits(grad, 'try-error')) {
            grad <- matrix(1e6, ncol = length(pars), nrow = 1)
        } else if (any(is.na(grad)) | !any(is.finite(grad))) {
            grad <- matrix(1e6, ncol = length(pars), nrow = 1)
        }
        return(grad)
    }
    if (use_hessian) {
        hess_fun <- function(pars, fun, tsdistenv)
        {
            names(pars) <- tsdistenv$parameter_names
            return(fun$he(pars, atomic = fun$env$usingAtomics()))
        }
    } else {
        hess_fun <- NULL
    }
    out <- list(data = data_list, par_list = par_list, map = map, llh_fun = llh_fun,
             grad_fun = grad_fun, hess_fun = hess_fun, parameter_names = parameter_names, 
             lower = setup[estimate == 1]$lower, upper = setup[estimate == 1]$upper)
    return(out)
}

score_function <- function(object, pars, use_hessian)
{
    y <- object$target$y
    object$target$y <- y[1]
    spec_list <- tsdistribution_tmb(object, use_hessian = use_hessian)
    silent <- TRUE
    fun <- try(MakeADFun(data = spec_list$data, parameters = spec_list$par_list, DLL = "tsdistributions_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent), silent = FALSE)
    m <- length(fun$par)
    n <- length(y)
    jac <- matrix(0, ncol = m, nrow = n)
    for (i in 1:n) {
        fun$env$data$y <- y[i]
        jac[i,] <- fun$gr(pars)
    }
    return(jac)    
}
