// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_dghst
NumericVector c_dghst(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dghst(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dghst(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_rghst
NumericVector c_rghst(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rghst(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rghst(n, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dghyp
NumericVector c_dghyp(NumericVector x, NumericVector alpha, NumericVector beta, NumericVector delta, NumericVector mu, NumericVector lambda, int logr);
RcppExport SEXP _tsdistributions_c_dghyp(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP deltaSEXP, SEXP muSEXP, SEXP lambdaSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dghyp(x, alpha, beta, delta, mu, lambda, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_dgh
NumericVector c_dgh(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, NumericVector lambda, int logr);
RcppExport SEXP _tsdistributions_c_dgh(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP lambdaSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dgh(x, mu, sigma, skew, shape, lambda, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_dnig
NumericVector c_dnig(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dnig(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dnig(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_rstd
NumericVector c_rstd(int n, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rstd(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rstd(n, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dstd
NumericVector c_dstd(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dstd(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dstd(x, mu, sigma, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_pstd
NumericVector c_pstd(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_pstd(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pstd(q, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_qstd
NumericVector c_qstd(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_qstd(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qstd(p, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_rsstd
NumericVector c_rsstd(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rsstd(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rsstd(n, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dsstd
NumericVector c_dsstd(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dsstd(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dsstd(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_psstd
NumericVector c_psstd(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_psstd(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_psstd(q, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_qsstd
NumericVector c_qsstd(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_qsstd(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qsstd(p, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_djsu
NumericVector c_djsu(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_djsu(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_djsu(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_qjsu
NumericVector c_qjsu(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_qjsu(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qjsu(p, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_pjsu
NumericVector c_pjsu(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_pjsu(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pjsu(q, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_rjsu
NumericVector c_rjsu(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rjsu(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rjsu(n, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_rsnorm
NumericVector c_rsnorm(int n, NumericVector mu, NumericVector sigma, NumericVector skew);
RcppExport SEXP _tsdistributions_c_rsnorm(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rsnorm(n, mu, sigma, skew));
    return rcpp_result_gen;
END_RCPP
}
// c_dsnorm
NumericVector c_dsnorm(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, int logr);
RcppExport SEXP _tsdistributions_c_dsnorm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dsnorm(x, mu, sigma, skew, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_psnorm
NumericVector c_psnorm(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew);
RcppExport SEXP _tsdistributions_c_psnorm(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    rcpp_result_gen = Rcpp::wrap(c_psnorm(q, mu, sigma, skew));
    return rcpp_result_gen;
END_RCPP
}
// c_qsnorm
NumericVector c_qsnorm(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew);
RcppExport SEXP _tsdistributions_c_qsnorm(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qsnorm(p, mu, sigma, skew));
    return rcpp_result_gen;
END_RCPP
}
// c_rged
NumericVector c_rged(int n, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rged(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rged(n, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dged
NumericVector c_dged(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dged(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dged(x, mu, sigma, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_pged
NumericVector c_pged(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_pged(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pged(q, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_qged
NumericVector c_qged(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector shape);
RcppExport SEXP _tsdistributions_c_qged(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qged(p, mu, sigma, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_rsged
NumericVector c_rsged(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_rsged(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rsged(n, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dsged
NumericVector c_dsged(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dsged(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dsged(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_psged
NumericVector c_psged(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_psged(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_psged(q, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_qsged
NumericVector c_qsged(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape);
RcppExport SEXP _tsdistributions_c_qsged(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qsged(p, mu, sigma, skew, shape));
    return rcpp_result_gen;
END_RCPP
}
// c_dhyp
NumericVector c_dhyp(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr);
RcppExport SEXP _tsdistributions_c_dhyp(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP skewSEXP, SEXP shapeSEXP, SEXP logrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type skew(skewSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type logr(logrSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dhyp(x, mu, sigma, skew, shape, logr));
    return rcpp_result_gen;
END_RCPP
}
// c_rghyp
NumericVector c_rghyp(int n, double mu, double delta, double alpha, double beta, double lambda);
RcppExport SEXP _tsdistributions_c_rghyp(SEXP nSEXP, SEXP muSEXP, SEXP deltaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rghyp(n, mu, delta, alpha, beta, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tsdistributions_c_dghst", (DL_FUNC) &_tsdistributions_c_dghst, 6},
    {"_tsdistributions_c_rghst", (DL_FUNC) &_tsdistributions_c_rghst, 5},
    {"_tsdistributions_c_dghyp", (DL_FUNC) &_tsdistributions_c_dghyp, 7},
    {"_tsdistributions_c_dgh", (DL_FUNC) &_tsdistributions_c_dgh, 7},
    {"_tsdistributions_c_dnig", (DL_FUNC) &_tsdistributions_c_dnig, 6},
    {"_tsdistributions_c_rstd", (DL_FUNC) &_tsdistributions_c_rstd, 4},
    {"_tsdistributions_c_dstd", (DL_FUNC) &_tsdistributions_c_dstd, 5},
    {"_tsdistributions_c_pstd", (DL_FUNC) &_tsdistributions_c_pstd, 4},
    {"_tsdistributions_c_qstd", (DL_FUNC) &_tsdistributions_c_qstd, 4},
    {"_tsdistributions_c_rsstd", (DL_FUNC) &_tsdistributions_c_rsstd, 5},
    {"_tsdistributions_c_dsstd", (DL_FUNC) &_tsdistributions_c_dsstd, 6},
    {"_tsdistributions_c_psstd", (DL_FUNC) &_tsdistributions_c_psstd, 5},
    {"_tsdistributions_c_qsstd", (DL_FUNC) &_tsdistributions_c_qsstd, 5},
    {"_tsdistributions_c_djsu", (DL_FUNC) &_tsdistributions_c_djsu, 6},
    {"_tsdistributions_c_qjsu", (DL_FUNC) &_tsdistributions_c_qjsu, 5},
    {"_tsdistributions_c_pjsu", (DL_FUNC) &_tsdistributions_c_pjsu, 5},
    {"_tsdistributions_c_rjsu", (DL_FUNC) &_tsdistributions_c_rjsu, 5},
    {"_tsdistributions_c_rsnorm", (DL_FUNC) &_tsdistributions_c_rsnorm, 4},
    {"_tsdistributions_c_dsnorm", (DL_FUNC) &_tsdistributions_c_dsnorm, 5},
    {"_tsdistributions_c_psnorm", (DL_FUNC) &_tsdistributions_c_psnorm, 4},
    {"_tsdistributions_c_qsnorm", (DL_FUNC) &_tsdistributions_c_qsnorm, 4},
    {"_tsdistributions_c_rged", (DL_FUNC) &_tsdistributions_c_rged, 4},
    {"_tsdistributions_c_dged", (DL_FUNC) &_tsdistributions_c_dged, 5},
    {"_tsdistributions_c_pged", (DL_FUNC) &_tsdistributions_c_pged, 4},
    {"_tsdistributions_c_qged", (DL_FUNC) &_tsdistributions_c_qged, 4},
    {"_tsdistributions_c_rsged", (DL_FUNC) &_tsdistributions_c_rsged, 5},
    {"_tsdistributions_c_dsged", (DL_FUNC) &_tsdistributions_c_dsged, 6},
    {"_tsdistributions_c_psged", (DL_FUNC) &_tsdistributions_c_psged, 5},
    {"_tsdistributions_c_qsged", (DL_FUNC) &_tsdistributions_c_qsged, 5},
    {"_tsdistributions_c_dhyp", (DL_FUNC) &_tsdistributions_c_dhyp, 6},
    {"_tsdistributions_c_rghyp", (DL_FUNC) &_tsdistributions_c_rghyp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_tsdistributions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}