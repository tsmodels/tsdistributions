#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

// Declare all the functions
double dnorm_std(const double x);
double signum(const double x);
double dhuge(void);
double heaviside(const double x, const double a);
double depsilon(void);
double kappagh(const double x, const double lambda);
double deltakappagh(const double x, const double lambda);
Rcpp::NumericVector paramgh(const double rho, const double zeta, const double lambda);
Rcpp::NumericVector paramghst(const double betabar, const double nu);
double dghst_std(const double x, const double betabar, const double nu);
Rcpp::NumericVector c_dghst(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double rsghst_std(const double betabar, const double nu);
Rcpp::NumericVector c_rghst(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double dgh(const double x, const double alpha, const double beta, const double delta, const double mu, const double lambda);
Rcpp::NumericVector c_dghyp(Rcpp::NumericVector x, Rcpp::NumericVector alpha, Rcpp::NumericVector beta, Rcpp::NumericVector delta, Rcpp::NumericVector mu, Rcpp::NumericVector lambda, int logr);
double dgh_std(const double x, const double rho, const double zeta, const double lambda);
Rcpp::NumericVector c_dgh(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, Rcpp::NumericVector lambda, int logr);
double dnig_std(const double x, const double rho, const double zeta);
Rcpp::NumericVector c_dnig(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double rstd_std(const double nu);
Rcpp::NumericVector c_rstd(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double xdt(const double x, const double nu);
double dstd_std(const double x, const double nu);
Rcpp::NumericVector c_dstd(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape, int logr);
double pstd_std(const double q, const double mu, const double sigma, const double nu);
Rcpp::NumericVector c_pstd(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double qstd_std(const double p, const double mu, const double sigma, const double nu);
Rcpp::NumericVector c_qstd(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double rsstd_std(const double xi, const double nu);
Rcpp::NumericVector c_rsstd(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double dsstd_std(const double x, const double xi, const double nu);
Rcpp::NumericVector c_dsstd(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double psstd_std(const double q, const double mu, const double sigma, const double xi, const double nu);
Rcpp::NumericVector c_psstd(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double qsstd_std(const double p, const double xi, const double nu);
Rcpp::NumericVector c_qsstd(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double djsu_std(const double x, const double nu, const double tau);
Rcpp::NumericVector c_djsu(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double qjsu_std(const double p, const double nu, const double tau);
Rcpp::NumericVector c_qjsu(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double pjsu_std(const double q, const double mu, const double sigma, const double nu, const double tau);
Rcpp::NumericVector c_pjsu(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double rsnorm_std(const double xi);
Rcpp::NumericVector c_rsnorm(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew);
double dsnorm_std(const double x, const double xi);
Rcpp::NumericVector c_dsnorm(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, int logr);
double psnorm_std(const double q, const double mu, const double sigma, const double xi);
Rcpp::NumericVector c_psnorm(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew);
double qsnorm_std(const double p, const double xi);
Rcpp::NumericVector c_qsnorm(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew);
double rged_std(const double nu);
Rcpp::NumericVector c_rged(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double dged_std(const double x, const double nu);
Rcpp::NumericVector c_dged(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape, int logr);
double pged_std(const double q, const double mu, const double sigma, const double nu);
Rcpp::NumericVector c_pged(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double qged_std(const double p, const double shape);
Rcpp::NumericVector c_qged(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector shape);
double rsged_std(const double xi, const double nu);
Rcpp::NumericVector c_rsged(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double dsged_std(const double x, const double xi, const double nu);
Rcpp::NumericVector c_dsged(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double psged_std(const double q, const double mu, const double sigma, const double xi, const double nu);
Rcpp::NumericVector c_psged(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double qsged_std(const double p, const double xi, const double nu);
Rcpp::NumericVector c_qsged(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape);
double dhyp(const double x, const double alpha, const double beta, const double delta, const double mu);
double dhyp_std(const double x, const double rho, const double zeta);
Rcpp::NumericVector c_dhyp(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector skew, Rcpp::NumericVector shape, int logr);
double g_function(double y, double beta, double m, double lambda);
double f_function(double y, double beta, double m, double lambda);
Rcpp::NumericVector rgig1(int n, Rcpp::NumericVector param);
Rcpp::NumericVector rgig(int n, Rcpp::NumericVector param);
Rcpp::NumericVector c_rghyp(int n, double mu, double delta, double alpha, double beta, double lambda);

#endif
