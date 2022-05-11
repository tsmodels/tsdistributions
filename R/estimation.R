#' Estimates the parameters of a distribution using autodiff.
#'
#' @param object an object of class tsdistribution.spec.
#' @param solver only \dQuote{nlminb} currently supported.
#' @param control solver control parameters.
#' @param use_hessian whether to use the hessian in the calculation.
#' @param ... additional parameters passed to the estimation function
#' @return A list of coefficients and other information.
#' @method estimate tsdistribution.spec
#' @details The estimation makes use of the TMB package for minimizing
#' the negative of the log-likelihood using automatic differentiation.
#' @aliases estimate
#' @rdname estimate
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
 tsgarchenv <- new.env()
 tsgarchenv$llh <- 1
 tsgarchenv$grad <- NULL
 tsgarchenv$parameter_names <- names(fun$par)
 tsgarchenv$estimation_names <- object$parmatrix[estimate == 1]$parameter
 tsgarchenv$parmatrix <- object$parmatrix
 if (solver == "nlminb") {
  sol <- nlminb(start = fun$par, objective = spec_list$llh_fun,
                gradient = spec_list$grad_fun, hessian = spec_list$hess_fun,
                lower = spec_list$lower,
                upper = spec_list$upper,  control = control,
                fun = fun, tsgarchenv = tsgarchenv)
  pars <- sol$par
 } else if (solver == "optim") {
  sol <- optim(par = fun$par, fn = spec_list$llh_fun, gr = spec_list$grad_fun,
               lower = spec_list$lower, upper = spec_list$upper, method = "L-BFGS-B", fun = fun, tsgarchenv = tsgarchenv)
  pars <- sol$par
 }
 parmatrix <- object$parmatrix
 parmatrix[estimate == 1]$value <- pars
 par_list <- as.list(parmatrix$value)
 names(par_list) <- parmatrix$parameter
 object$parmatrix <- parmatrix
 llh <- spec_list$llh_fun(pars, fun, tsgarchenv)
 gradient <- spec_list$grad_fun(pars, fun, tsgarchenv)
 hessian <- spec_list$hess_fun(pars, fun, tsgarchenv)
 names(pars) <- tsgarchenv$estimation_names
 colnames(gradient) <- tsgarchenv$estimation_names
 colnames(hessian) <- rownames(hessian) <- tsgarchenv$estimation_names
 out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, solver_out = sol, spec = object)
 class(out) <- "tsdistribution.estimate"
 return(out)
}


# convert specification to object ready for estimation using TMB
tsdistribution_tmb <- function(object, use_hessian = FALSE, silent = TRUE)
{
 setup <- object$parmatrix
 data_list <- list(y = object$target$y, dclass = distribution_class(object$distribution), model = "distribution")
 if (any(setup$include == 0)) {
  fixed_values <- setup[estimate == 0]$parameter
  fixed <- vector(mode = "list", length = length(fixed_values))
  names(fixed) <- fixed_values
  for (i in 1:length(fixed)) fixed[[i]] <- factor(NA)
 } else {
  fixed <- NULL
 }
 par_list <- as.list(as.numeric(setup$value))
 names(par_list) <- setup$parameter
 # create the functions
 parameter_names <- setup[estimate == 1]$parameter
 llh_fun <- function(pars, fun, tsgarchenv)
 {
  names(pars) <- tsgarchenv$parameter_names
  llh <- tsgarchenv$llh
  llh <- try(-fun$fn(pars), silent = TRUE)
  if (inherits(llh, 'try-error')) {
   llh <- llh + 0.2 * abs(llh)
  } else if (is.na(llh) | !is.finite(llh)) {
   llh <- llh + 0.2 * abs(llh)
  }
  tsgarchenv$llh <- llh
  return(llh)
 }
 
 grad_fun <- function(pars, fun, tsgarchenv)
 {
  names(pars) <- tsgarchenv$parameter_names
  grad <- try(-fun$gr(pars), silent = TRUE)
  if (inherits(grad, 'try-error')) {
    grad <- matrix(1e6, ncol = length(pars), nrow = 1)
  } else if (any(is.na(grad)) | !any(is.finite(grad))) {
    grad <- matrix(1e6, ncol = length(pars), nrow = 1)
  }
  return(grad)
 }
 if (use_hessian) {
  hess_fun <- function(pars, fun, tsgarchenv)
  {
   names(pars) <- tsgarchenv$parameter_names
   return(-fun$he(pars))
  }
 } else {
  hess_fun <- NULL
 }
 
 out <- list(data = data_list, par_list = par_list,
             map = fixed, llh_fun = llh_fun,
             grad_fun = grad_fun, hess_fun = hess_fun,
             parameter_names = parameter_names, lower = setup[estimate == 1]$lower,
             upper = setup[estimate == 1]$upper)
 return(out)
}

.make_standard_errors <- function(pmatrix, H)
{
 pars <- pmatrix[estimate == 1]$value
 se <- sqrt(diag(solve(H)))
 tvalues <- pars/se
 pvalues <- 2*(1 - pnorm(abs(tvalues)))
 return(data.frame("Std. Error" = se,"t value" = tvalues, "Pr(>|t|)" = pvalues, check.names = FALSE))
}
