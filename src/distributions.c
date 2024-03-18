# include "distributions.h"
# include <R.h>
# include <limits.h>
# include <math.h>
# include <Rmath.h>

/*
 * Key Functions
 */
double dnorm_std(const double x)
{
  double pdf;
  pdf = exp (-0.5 * x * x) / sqrt (2.0 * PI);
  if(pdf == 0.0) pdf = 0.0 + 2.22507e-24;
  return pdf;
}

double signum(const double x)
{
	double res = -(x < 0.0) + (x > 0.0);
	return res;
}

double dhuge(void)
{
  return HUGE_VAL;
}

double heaviside(const double x, const double a){
    double res = (signum(x - a) + 1.0) / 2.0;
	return res;
}

double depsilon(void)
{
  double r;
  r = 1.0;
  while( 1.0 < (double)(1.0 + r))
  {
      r = r / 2.0;
  }
  double res = 2.0 * r;
  return res;
}

double kappagh(const double x, const double lambda)
{
	double kappa = 0.0;
	if(lambda == -0.5) {
	    kappa = 1.0 / x;
	} else {
	    kappa = (bessel_k(x, lambda + 1.0, 2.0) / bessel_k(x, lambda, 2.0)) / x;
	}
	return kappa;
}

double deltakappagh(const double x, const double lambda)
{
	double deltakappa = 0.0;
    deltakappa = kappagh(x, lambda + 1.0) - kappagh(x, lambda);
	return deltakappa;
}

double* paramgh(const double rho, const double zeta, const double lambda)
{
	double *param = malloc(4 * sizeof(double));
	double rho2 = 1.0 - pow(rho, 2.0);
	double alpha = pow(zeta, 2.0) * kappagh(zeta, lambda) / rho2;
	alpha = alpha * (1.0 + pow(rho, 2.0) * pow(zeta, 2.0) * deltakappagh(zeta, lambda) / rho2);
	alpha = sqrt(alpha);
	double beta = alpha * rho;
	double delta = zeta / (alpha * sqrt(rho2));
	double mu = -beta * pow(delta, 2.0) * kappagh(zeta, lambda);
	param[0] = (double) alpha;
	param[1] = (double) beta;
	param[2] = (double) delta;
	param[3] = (double) mu;
	return param;
}

double* paramghst(const double betabar, const double nu)
{
	double *param = malloc(4 * sizeof(double));
	double delta = sqrt(1.0 / (((2.0 * betabar * betabar) / ((nu - 2.0) * (nu - 2.0) * (nu - 4.0))) + (1.0/(nu - 2.0))));
	double beta = betabar / delta;
	double mu = -1.0 * ((beta * (delta * delta)) / (nu - 2.0));
	param[0] = (double) nu;
	param[1] = (double) beta;
	param[2] = (double) delta;
	param[3] = (double) mu;
	return param;
}

/*
 * GH Skew Student Distribution
 */
double dghst_std(const double x, const double betabar, const double nu)
{
	double *param;
	param = paramghst(betabar, nu);
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	double betasqr = beta * beta;
	double deltasqr = delta * delta;
	double res = x - mu;
	double pdf = ((1.0 - nu) / 2.0) * log(2) + nu * log(delta) + ((nu + 1.0) / 2.0) * log(fabs(beta)) + 
	    log(bessel_k(sqrt(betasqr * (deltasqr + res * res)), (nu + 1.0)/2.0, 2.0)) - sqrt(betasqr * 
	    (deltasqr + res * res)) + beta * res - lgammafn(nu/2.0) - log(PI)/2.0 - ((nu + 1.0) / 2.0) * 
	    log(deltasqr + res * res) / 2.0;
	free(param);
	pdf = exp(pdf);
	return pdf;
}

void c_dghst(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dghst_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double rsghst_std(const double betabar, const double nu)
{
	// Note: for values of nu<5 it is likely that sd(r) will have a very large variance
	// Existence of moment conditions (vs lower parameter bounds) are defined in the paper.
	double *param;
	param = paramghst(betabar, nu);
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	double y = 1.0 / rgamma(nu / 2.0, 2.0 / (delta * delta));
	double sigma = sqrt(y);
	double z = rnorm(0,1);
	double ans =  mu + beta * sigma * sigma + sigma * z;
	free(param);
	return ans;
}

void c_rghst(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rsghst_std(skew[i], shape[i]) * sigma[i];
	}
	PutRNGstate();
}

/*
 * GH Distribution
 */
double dgh(const double x, const double alpha, const double beta, const double delta, const double mu, const double lambda)
{
	double pdf = 0.0;
	if(alpha <= 0.0) {
		return pdf = 0.0;
	}
	if(delta <= 0.0) {
		return pdf = 0.0;
	}
	if(fabs(beta) >= alpha) {
		return pdf = 0.0;
	}
	double arg = delta * sqrt(pow(alpha, 2.0) - pow(beta, 2.0));
	double a = (lambda / 2.0) * log(pow(alpha, 2.0) - pow(beta, 2.0)) - (log(sqrt(2.0 * PI)) + (lambda - 0.5) * log(alpha) + 
             lambda * log(delta) + log(bessel_k(arg, lambda, 2.0)) - arg);
	double f = ((lambda - 0.5) / 2.0) * log(pow(delta, 2.0) + pow((x - mu), 2.0));
	arg = alpha * sqrt(pow(delta, 2.0) + pow((x - mu), 2.0));
	double k = log(bessel_k(arg, lambda - 0.5, 2.0)) - arg;
	double e = beta * (x - mu);
	pdf = exp(a + f + k + e);
	return pdf;
}

void c_dghyp(double *x, double *alpha, double *beta, double *delta, double *mu, double *lambda, double *ans, int *n, int *logr)
{
    int i;
    for(i = 0; i<*n; i++)
    {
        ans[i] = dgh(x[i], alpha[i], beta[i], delta[i], mu[i], lambda[i]);
        if(*logr == 1) ans[i] = log(ans[i]);
    }
}


double dgh_std(const double x, const double rho, const double zeta, const double lambda)
{
	double pdf;
	double *param;
	param = paramgh(rho, zeta, lambda);
	double alpha = param[0];
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	pdf = dgh(x, alpha, beta, delta, mu, lambda);
	free(param);
	return pdf;
}

void c_dgh(double *x, double *mu, double *sigma, double *skew, double *shape, double *lambda, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dgh_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i], lambda[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}


/*
 * NIG Distribution
 */

double dnig_std(const double x, const double rho, const double zeta)
{
	double pdf = 0.0;
	double lambda = -0.5;
	double *param;
	param = paramgh(rho, zeta, lambda);
	double alpha = param[0];
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	double deltasq = delta * delta;
	double res = x - mu;
	pdf = -log(PI) + log(alpha) + log(delta) + log(bessel_k(alpha * sqrt(deltasq + res * res), 1, 1)) + 
	    delta * sqrt(alpha * alpha - beta * beta) + beta * res -0.5 * log(deltasq + res * res);
	pdf = exp(pdf);
	free(param);
	return pdf;
}

void c_dnig(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dnig_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}


/*
 * Student Distribution
 */

double rstd_std(const double nu)
{
	double ans = 0.0;
	if(nu > 2.0) {
	    double s = sqrt(nu / (nu - 2.0));
	    ans = rt(nu) * 1.0 / s;
	}
	return ans;
}

void c_rstd(int *n, double *mu, double *sigma, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rstd_std(shape[i]) * sigma[i];
	}
	PutRNGstate();
}

double xdt(const double x, const double nu)
{
	double a, b, pdf;
	a = gammafn((nu + 1.0)/2.0) / sqrt(PI * nu);
    b = gammafn(nu / 2.0) * pow((1.0 + (x * x) / nu), ((nu + 1.0) / 2.0));
    pdf = a / b;
	return pdf;
}

double dstd_std(const double x, const double nu)
{
	double pdf, s;
	if(nu <= 2) {
		pdf = 999;
	} else {
		s = sqrt(nu / (nu - 2.0));
		pdf = s * xdt(x * s, nu);
	}
	return pdf;
}

void c_dstd(double *x, double *mu, double *sigma, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dstd_std((x[i] - mu[i]) / sigma[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double pstd_std(const double q, const double mu, const double sigma, const double nu)
{
	double s = sqrt(nu / (nu - 2.0));
	double z = (q - mu) / sigma;
	double p = pt(z * s, nu, 1, 0);
	return p;
}

void c_pstd(double *q, double *mu, double *sigma, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = pstd_std(q[i], mu[i], sigma[i], shape[i]);
	}
}

double qstd_std(const double p, const double mu, const double sigma, const double nu)
{
	double s = sqrt(nu / (nu - 2.0));
	double q = qt(p, nu, 1, 0) * sigma / s + mu;
	return q;
}

void c_qstd(double *p, double *mu, double *sigma, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = qstd_std(p[i], mu[i], sigma[i], shape[i]);
	}
}

/*
 * Skew Student Distribution (Fernandez & Steel)
 */
double rsstd_std(const double xi, const double nu)
{
	double weight, z, rr, m1, mu, sigma, xx, ans;
	ans = 0.0;
	weight = xi / (xi + 1.0 / xi);
	z = runif(-1.0 * weight, 1.0 - weight);
	xx = (z < 0) ? 1.0 / xi : xi;
	rr = -1.0 * fabs(rstd_std(nu)) / xx * sign(z);
	m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta(0.5, 0.5 * nu);
	mu = m1 * (xi - 1.0 / xi);
	sigma =  sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0 / (xi * xi)) + 2 * (m1 * m1) - 1.0);
	ans =  (rr - mu) / sigma;
	return ans;
}

void c_rsstd(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rsstd_std(skew[i], shape[i]) * sigma[i];
	}
	PutRNGstate();
}

double dsstd_std(const double x, const double xi, const double nu)
{
	double mu, m1,beta, sigma, z, g,pdf,a,b, xxi;
	xxi = xi;
	a = 0.5;
	b = nu / 2.0;
	beta = (gammafn(a) / gammafn(a + b)) * gammafn(b);
	m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta;
	mu = m1 * (xi - 1.0 / xi);
	sigma = sqrt((1.0 - pow(m1, 2.0)) * (pow(xi, 2.0) + 1.0 / (pow(xi, 2.0))) + 2.0 * pow(m1, 2.0) - 1.0);
	z = x * sigma + mu;
	if(z == 0) {
		xxi = 1.0;
	}
	if(z < 0) {
		xxi = 1.0 / xi;
	}
	g = 2.0 / (xi + 1.0 / xi);
	pdf = g * dstd_std(z / xxi, nu) * sigma;
	return pdf;
}

void c_dsstd(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dsstd_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double psstd_std(const double q, const double mu, const double sigma, const double xi, const double nu)
{
	double qx = (q - mu) / sigma;
	double m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta(0.5, nu / 2.0);
	double mux = m1 * (xi - 1.0 / xi);
	double sig =  sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double z = qx * sig + mux;
	double Xi = (z < 0) ? 1.0 / xi : xi;
	double g = 2.0 / (xi + 1.0 / xi);
	double p = heaviside(z, 0) - signum(z) * g * Xi * pstd_std(-fabs(z) / Xi, 0.0, 1.0, nu);
	return p;
}

void c_psstd(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = psstd_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
	}
}

double qsstd_std(const double p, const double xi, const double nu)
{
	double m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta(0.5, nu / 2.0);
	double mu = m1 * (xi - 1.0 / xi);
	double sigma =  sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double g = 2.0 / (xi + 1.0 / xi);
	double z = p - (1.0/(1.0 + xi * xi));
	double Xi = pow(xi, signum(z));
	double tmp = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
	double q = (-signum(z) * qstd_std(tmp, 0.0, 1.0, nu) * Xi - mu) / sigma;
	return q;
}

void c_qsstd(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = qsstd_std(p[i], skew[i], shape[i]) * sigma[i] + mu[i];
	}
}

/*
 * Johnson's SU Distribution
 */
//nu = skew, tau = shape (!)
double djsu_std(const double x, const double nu, const double tau)
{
	double w, z, r, omega, c, pdf = 0.0;
	double rtau = 1.0 / tau;
	if(rtau < 0.0000001) {
		w = 1.0;
	} else {
		w = exp(rtau * rtau);
	}
	omega= -nu * rtau;
	c = sqrt(1.0 / (0.5 * (w - 1.0) * (w * cosh(2.0 * omega) + 1.0)));
	z = (x - (c * sqrt(w) * sinh(omega))) / c;
	r = -nu + asinh(z) / rtau;
	pdf = -log(c) - log(rtau) - 0.5 * log(z * z + 1.0) - 0.5 * log(2.0 * PI) - 0.5 * r * r;
	pdf = exp(pdf);
	return pdf;
}

void c_djsu(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = djsu_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double qjsu_std(const double p, const double nu, const double tau)
{
	double rtau, rr, z, w, omega, cc, ans;
	ans = 0.0;
	rtau = 1.0 / tau;
	rr = qnorm(p, 0.0, 1.0, 1, 0);
	z = sinh(rtau * (rr + nu));
	w = (rtau < 0.0000001) ? 1 : exp(rtau * rtau);
	omega = -1.0 * nu * rtau;
	cc = sqrt(1.0 / (0.5 * (w - 1.0) * (w * cosh(2.0 * omega) + 1.0)));
	ans = (cc * sqrt(w) * sinh(omega)) + cc * z;
	return ans;
}

void c_qjsu(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + qjsu_std(p[i], skew[i], shape[i]) * sigma[i];
	}
}

double rjsu_std(const double nu, const double tau)
{
	double x, ans;
	ans = 0.0;
	x = runif(0, 1);
	ans = qjsu_std(x, nu, tau);
	return ans;
}

void c_rjsu(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rjsu_std(skew[i], shape[i]) * sigma[i];
	}
	PutRNGstate();
}

double pjsu_std(const double q, const double mu, const double sigma, const double nu, const double tau)
{
	double rtau = 1.0 / tau;
	double w = (rtau < 0.0000001) ? 1.0 : exp(rtau * rtau);
	double omega = -1.0 * nu * rtau;
	double c = 1.0 / sqrt(0.5 * (w - 1.0) * (w * cosh(2.0 * omega) + 1.0));
	double z = (q - (mu + c * sigma * sqrt(w) * sinh(omega))) / (c * sigma);
	double r = -1.0 * nu + asinh(z) / rtau;
	double p = pnorm(r, 0, 1, 1, 0);
	return p;
}

void c_pjsu(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = pjsu_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
	}
}

/*
 * Skew Normal Distribution
 */
double rsnorm_std(const double xi)
{
	double weight, z, rr, m1, mu, sigma, xx, ans;
	weight = xi / (xi + 1.0/xi);
	z = runif(-weight, 1.0 - weight);
	xx = (z < 0) ? 1.0 / xi : xi;
	rr = -1.0 * fabs(rnorm(0, 1)) / xx * sign(z);
	m1 = 2.0 / sqrt(2.0 * PI);
	mu = m1 * (xi - 1.0 / xi);
	sigma = sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0 / (xi * xi) ) + 2 * (m1 * m1) - 1.0);
	ans = (rr - mu ) / sigma;
	return ans;
}

void c_rsnorm(int *n, double *mu, double *sigma, double *skew, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rsnorm_std(skew[i]) * sigma[i];
	}
	PutRNGstate();
}

double dsnorm_std(const double x, const double xi)
{
	double pdf;
	double mu, sigma,z, xxi, g;
	double m1 = 2.0 / sqrt(2.0 * PI);
	double m12 = m1 * m1;
	double xi2 = xi * xi;
	mu = m1 * (xi - 1.0 / xi);
	sigma = sqrt((1.0 - m12) * (xi2 + 1.0 / xi2) + 2.0 * m12 - 1.0);
	z = x * sigma + mu;
	xxi = xi;
	if(z == 0.0) {
	 xxi = 1.0;
	}
	if(z < 0.0) {
	 xxi = 1.0 / xi;
	}
	g = 2.0 / (xi + 1.0 / xi);
	pdf = g * dnorm_std(z / xxi) * sigma;
	return pdf;
}

void c_dsnorm(double *x, double *mu, double *sigma, double *skew, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dsnorm_std((x[i] - mu[i]) / sigma[i], skew[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double psnorm_std(const double q, const double mu, const double sigma, const double xi)
{
	double qx = (q - mu) / sigma;
	double m1 = 2.0 / sqrt(2.0 * PI);
	double mux = m1 * (xi - 1.0 / xi);
	double sig = sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double z = qx * sig + mux;
	double Xi = (z < 0) ? 1.0 / xi : xi;
	double g = 2.0 / (xi + 1.0 / xi);
	double p = heaviside(z, 0) - signum(z) * g * Xi * pnorm(-fabs(z)/Xi, 0, 1, 1, 0);
	return p;
}

void c_psnorm(double *q, double *mu, double *sigma, double *skew, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = psnorm_std(q[i], mu[i], sigma[i], skew[i]);
	}
}

double qsnorm_std(const double p, const double xi)
{
	double m1 = 2.0 / sqrt(2.0 * PI);
	double mu = m1 * (xi - 1.0 / xi);
	double sigma = sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double g = 2.0 / (xi + 1.0 / xi);
	double z = p - (1.0 / (1.0 + xi * xi));
	double Xi = pow(xi, signum(z));
	double tmp = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
	double q = (-1.0 * signum(z) * qnorm(tmp, 0, Xi, 1, 0) - mu) / sigma;
	return q;
}

void c_qsnorm(double *p, double *mu, double *sigma, double *skew, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + qsnorm_std(p[i], skew[i]) * sigma[i];
	}
}

/*
 * Generalized Error Distribution
 */
double rged_std(const double nu)
{
	double lambda, rr, ans;
	ans = 0.0;
	lambda = sqrt(pow(0.5, 2.0 / nu) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
	rr = rgamma(1.0 / nu, 1.0);
	ans =  lambda * pow(2.0 * rr, 1.0 / nu) * sign(runif(0, 1) - 0.5);
	return ans;
}

void c_rged(int *n, double *mu, double *sigma, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rged_std(shape[i]) * sigma[i];
	}
	PutRNGstate();
}

double dged_std(const double x, const double nu)
{
	double lambda, g, pdf;
	lambda = sqrt(pow(1.0 / 2.0, 2.0 / nu) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
    g = nu / (lambda * (pow(2.0, 1.0 + (1.0 / nu))) * gammafn(1.0 / nu));
	pdf = g * exp(-0.5 * pow(fabs(x / lambda), nu));
	return pdf;
}

void c_dged(double *x, double *mu, double *sigma, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dged_std((x[i] - mu[i]) / sigma[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double pged_std(const double q, const double mu, const double sigma, const double nu)
{
	double qx = (q - mu) / sigma;
	double lambda = sqrt(1.0 / pow(2.0, (2.0 / nu)) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
	double g  = nu / (lambda * (pow(2.0, (1.0 + 1.0 / nu))) * gammafn(1.0 / nu));
	double h = pow(2.0, (1.0 / nu)) * lambda * g * gammafn(1.0 / nu) / nu;
	double s = 0.5 * pow(fabs(qx) / lambda , nu);
	double p = 0.5 + signum(qx) * h * pgamma(s, 1.0 / nu, 1, 1, 0);
	return p;
}

void c_pged(double *q, double *mu, double *sigma, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = pged_std(q[i], mu[i], sigma[i], shape[i]);
	}
}

double qged_std(const double p, const double shape)
{
	double y = 2.0 * p - 1.0;
	double lambda = sqrt(1.0 / pow(2.0, (2.0 / shape)) * gammafn(1.0 / shape) / gammafn(3.0 / shape));
	double q = lambda * pow(2.0 * qgamma(fabs(y), 1.0 / shape, 1, 1, 0), 1.0 / shape);
	q = q * signum(y);
	return q;
}

void c_qged(double *p, double *mu, double *sigma, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = qged_std(p[i], shape[i]) * sigma[i] + mu[i];
	}
}

/*
 * Skew Generalized Error Distribution (Fernandez & Steel)
 */
double rsged_std(const double xi, const double nu)
{
	double weight, lambda, z, rr, m1, mu, sigma, xx, ans;
	weight = xi / (xi + 1.0/xi);
	z = runif(-1.0 * weight, 1.0 - weight);
	xx = (z < 0) ? 1.0 / xi : xi;
	rr = -1.0 * fabs(rged_std(nu)) / xx * sign(z);
	lambda = sqrt (pow(0.5, 2.0/nu) * gammafn(1.0/nu) / gammafn(3.0/nu));
	m1 = pow(2.0, 1.0 / nu) * lambda * gammafn(2.0 / nu) / gammafn(1.0 / nu);
	mu = m1 * (xi - 1.0 / xi);
	sigma = sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0 / (xi* xi)) + 2.0 * (m1 * m1) - 1.0);
	ans = (rr - mu) / sigma;
	return ans;
}

void c_rsged(int *n, double *mu, double *sigma, double *skew, double *shape, double *ans)
{
	GetRNGstate();
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = mu[i] + rsged_std(skew[i], shape[i]) * sigma[i];
	}
	PutRNGstate();
}

double dsged_std(const double x, const double xi, const double nu)
{
	double lambda, m1, mu, sigma, z, g, pdf, xxi;
	xxi = xi;
	lambda = sqrt(pow(1.0 / 2.0, 2.0 / nu) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
	m1 = pow(2.0, 1.0 / nu) * lambda * gammafn(2.0 / nu) / gammafn(1.0 / nu);
	mu = m1 * (xi - 1.0 / xi);
	sigma = (1.0 - pow(m1, 2.0)) * (pow(xi, 2.0) + 1.0 / (pow(xi, 2.0))) + 2.0 * (pow(m1, 2.0)) - 1.0;
	sigma = sqrt(sigma);
	z = x * sigma + mu;
	if(z == 0) {
		xxi = 1.0;
	}
	if(z < 0) {
		xxi = 1 / xi;
	}
	g = 2.0 / (xi + 1.0 / xi);
	pdf = g * dged_std(z / xxi, nu) * sigma;
	return pdf;
}

void c_dsged(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dsged_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}

double psged_std(const double q, const double mu, const double sigma, const double xi, const double nu)
{
	double qx = (q - mu) / sigma;
	double lambda = sqrt (1.0/pow(2.0, 2.0 / nu) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
	double m1 = pow(2.0, 1.0 / nu) * lambda * gammafn(2.0 / nu) / gammafn(1.0 / nu);
	double mux = m1 * (xi - 1.0 / xi);
	double sig =  sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double z = qx * sig + mux;
	double Xi = (z < 0) ? 1.0 / xi : xi;
	double g = 2.0 / (xi + 1.0 / xi);
	double p = heaviside(z, 0) - signum(z) * g * Xi * pged_std(-fabs(z) / Xi, 0.0, 1.0, nu);
	return p;
}

void c_psged(double *q, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = psged_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
	}
}

double qsged_std(const double p, const double xi, const double nu)
{
	double lambda = sqrt (1.0 / pow(2.0, 2.0 / nu) * gammafn(1.0 / nu) / gammafn(3.0 / nu));
	double m1 = pow(2.0, 1.0 / nu) * lambda * gammafn(2.0 / nu) / gammafn(1.0 / nu);
	double mu = m1 * (xi - 1.0 / xi);
	double sigma = sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
	double g = 2.0 / (xi + 1.0 / xi);
	double z = p - (1.0 / (1.0 + xi * xi));
	double Xi = pow(xi, signum(z));
	double q = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
	q = (-signum(z) * qged_std(q, nu) * Xi - mu) / sigma;
	return q;
}

void c_qsged(double *p, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = qsged_std(p[i], skew[i], shape[i]) * sigma[i] + mu[i];
	}
}

/*
 * -----------------------------------------
 * Hypebolic Distribution
 * -----------------------------------------
 */
double dhyp(const double x, const double alpha, const double beta, const double delta, const double mu)
{
	double pdf = 0.0;
	if(alpha <= 0.0) {
		return pdf;
	}
	if(delta <= 0) {
		return pdf;
	}
	if(fabs(beta) >= alpha) {
		return pdf;
	}
	double g = alpha * alpha - beta * beta;
	double e = x - mu;
	pdf = 0.5 * log(g) - log(2.0 * alpha * delta * bessel_k(delta * sqrt(g), 1, 2)) - alpha * sqrt(delta * delta + e * e)  + beta * e;
	pdf = exp(pdf);
	return pdf;
}

double dhyp_std(const double x,  const double rho, const double zeta)
{
	double pdf;
	double *param;
	param = paramgh(rho, zeta, 1.0);
	double alpha = param[0];
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	pdf = dhyp(x, alpha, beta, delta, mu);
	free(param);
	return pdf;
}

void c_dhyp(double *x, double *mu, double *sigma, double *skew, double *shape, double *ans, int *n, int *logr)
{
	int i;
	for(i = 0; i<*n; i++) {
		ans[i] = dhyp_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
		if(*logr == 1) ans[i] = log(ans[i]);
	}
}
