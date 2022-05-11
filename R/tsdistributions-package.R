#' @rawNamespace useDynLib(tsdistributions, .registration=TRUE); useDynLib(tsdistributions_TMBExports)
#' @keywords internal
#' @importFrom TMB MakeADFun
#' @importFrom stats nlminb na.omit
#' @import data.table
#' @import methods
#' @importFrom zoo coredata
#' @importFrom knitr kable
#' @importFrom SkewHyperbolic pskewhyp qskewhyp
#' @importFrom GeneralizedHyperbolic ghypMom
#' @importFrom Rsolnp solnp
#' @importFrom stats dnorm integrate optim pnorm qnorm rnorm runif sd uniroot
#' @import tsmethods
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
