#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
/*
 * -----------------------------------------
 * Key Functions
 * -----------------------------------------
 */
double dnorm_std(const double );
double signum(const double );
double dhuge(void);
double Heaviside(const double , const double );
double depsilon(void);
double kappagh(const double , const double );
double deltakappagh(const double , const double );
double* paramgh(const double , const double , const double );
double* paramghskt(const double , const double );
/*
 * -----------------------------------------
 * GH Skew Student Distribution
 * -----------------------------------------
 */
double dghst_std(const double , const double , const double );
void c_dghst(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int* logr);
double rsghst_std(const double , const double );
void c_rghst(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans);
/*
 * -----------------------------------------
 * GH Distribution
 * -----------------------------------------
 */
// alpha-beta-delta-mu parameterization
double dgh(const double , const double , const double , const double , const double , const double);
void c_dghyp(double *, double *, double *, double *, double *, double *, double *, int *, int *);
// rho-zeta parameterization
double dgh_std(const double , const double , const double , const double);
void c_dgh(double *x, double *mu, double *sigma, double *skew, double *shape, double *lambda, double *ans, int *n, int *logr);
/*
 * -----------------------------------------
 * NIG Distribution
 * -----------------------------------------
 */
double dnig_std(const double , const double , const double);
void c_dnig(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr);
/*
 * -----------------------------------------
 * Student Distribution
 * -----------------------------------------
 */
double rstd_std(const double );
void c_rstd(int *n, double *mu, double *sigma, double *shape, double *ans);
double xdt(const double , const double );
double dstd_std(const double , const double );
void c_dstd(double *x, double *mu, double *sigma, double *shape, double *ans, int *n, int *logr);
double pstd_std(const double , const double , const double , const double );
void c_pstd(double *q, double *mu, double *sigma, double *shape, double *ans, int *n);
double qstd_std(const double , const double , const double , const double );
void c_qstd(double *p, double *mu, double *sigma, double *shape, double *ans, int *n);
/*
 * -----------------------------------------
 * Skew Student Distribution (Fernandez & Steel)
 * -----------------------------------------
 */
double rsstd_std(const double , const double );
void c_rsstd(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans);
double dsstd_std(const double , const double , const double );
void c_dsstd(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr);
double psstd_std(const double , const double , const double , const double , const double );
void c_psstd(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
double qsstd_std(const double , const double , const double );
void c_qsstd(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
/*
 * -----------------------------------------
 * Johnson's SU Distribution
 * -----------------------------------------
 */
double djsu_std(const double , const double , const double );
void c_djsu(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr);
double qjsu_std(const double , const double , const double );
void c_qjsu(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
double rjsu_std(const double n, const double );
void c_rjsu(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans);
double pjsu_std(const double , const double , const double , const double , const double );
void c_pjsu(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
/*
 * -----------------------------------------
 * Skew Normal Distribution
 * -----------------------------------------
 */
double rsnorm_std(const double );
void c_rsnorm(int *n, double *mu, double *sigma, double *skew, double *ans);
double dsnorm_std(const double , const double );
void c_dsnorm(double *x, double *mu, double *sigma, double *skew, double *ans, int *n, int *logr);
double psnorm_std(const double , const double , const double , const double );
void c_psnorm(double *q, double *mu, double *sigma, double *skew, double *ans, int *n);
double qsnorm_std(const double , const double );
void c_qsnorm(double *p, double *mu, double *sigma, double *skew, double *ans, int *n);
/*
 * -----------------------------------------
 * Generalized Error Distribution
 * -----------------------------------------
 */
double rged_std(const double );
void c_rged(int *n, double *mu, double *sigma, double *shape, double *ans);
double dged_std(const double , const double );
void c_dged(double *x, double *mu, double *sigma, double *shape, double *ans, int *n, int *logr);
double pged_std(const double , const double , const double , const double );
void c_pged(double *q, double *mu, double *sigma, double *shape, double *ans, int *n);
double qged_std(const double , const double );
void c_qged(double *p, double *mu, double *sigma, double *shape, double *ans, int *n);
/*
 * -----------------------------------------
 * Skew Generalized Error Distribution (Fernandez & Steel)
 * -----------------------------------------
 */
double rsged_std(const double , const double );
void c_rsged(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans);
double dsged_std(const double , const double , const double );
void c_dsged(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr);
double psged_std(const double , const double , const double , const double , const double );
void c_psged(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
double qsged_std(const double , const double, const double );
void c_qsged(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n);
/*
 * -----------------------------------------
 * Hypebolic Distribution
 * -----------------------------------------
 */
double dhyp(const double , const double , const double , const double , const double);
double dhyp_std(const double ,  const double , const double);
void c_dhyp(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr);
#endif /* DISTRIBUTIONS_H */
