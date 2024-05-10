#' @rawNamespace useDynLib(tsdistributions, .registration=TRUE); useDynLib(tsdistributions_TMBExports)
#' @keywords internal
#' @importFrom TMB MakeADFun
#' @import data.table
#' @import methods
#' @importFrom zoo coredata
#' @importFrom SkewHyperbolic pskewhyp qskewhyp
#' @importFrom GeneralizedHyperbolic ghypMom rghyp
#' @importFrom KernSmooth bkde
#' @importFrom Rsolnp solnp
#' @importFrom sandwich estfun bwNeweyWest vcovHAC vcovOPG bread
#' @importFrom stats dnorm integrate optim pnorm qnorm rnorm runif sd uniroot spline nlminb na.omit coef logLik printCoefmat ar lm logLik vcov AIC BIC ppoints approx rexp var
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom mev fit.gpd dgp pgp qgp rgp
#' @importFrom utils tail
#' @importFrom Rdpack reprompt
#' @import tsmethods
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
