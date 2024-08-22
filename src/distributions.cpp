# include "distributions.h"
using namespace Rcpp;

// Standard normal density function
double dnorm_std(const double x) {
    double pdf = std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
    if (pdf == 0.0) pdf = std::numeric_limits<double>::min(); // Smallest positive double
    return pdf;
}

// Sign function
double signum(const double x) {
    return -(x < 0.0) + (x > 0.0);
}

// Return a huge value (infinity)
double dhuge(void) {
    return R_PosInf;
}

// Heaviside step function
double heaviside(const double x, const double a) {
    return (signum(x - a) + 1.0) / 2.0;
}

// Machine epsilon
double depsilon(void) {
    double r = 1.0;
    while (1.0 < (1.0 + r)) {
        r /= 2.0;
    }
    return 2.0 * r;
}

// Kappa function for Generalized Hyperbolic Distribution
double kappagh(const double x, const double lambda) {
    double kappa = 0.0;
    if (lambda == -0.5) {
        kappa = 1.0 / x;
    } else {
        kappa = (R::bessel_k(x, lambda + 1.0, 2.0) / R::bessel_k(x, lambda, 2.0)) / x;
    }
    return kappa;
}

// Delta Kappa function
double deltakappagh(const double x, const double lambda) {
    return kappagh(x, lambda + 1.0) - kappagh(x, lambda);
}

// Parametrization for Generalized Hyperbolic Distribution
NumericVector paramgh(const double rho, const double zeta, const double lambda) {
    double rho2 = 1.0 - std::pow(rho, 2.0);
    double alpha = std::pow(zeta, 2.0) * kappagh(zeta, lambda) / rho2;
    alpha *= (1.0 + std::pow(rho, 2.0) * std::pow(zeta, 2.0) * deltakappagh(zeta, lambda) / rho2);
    alpha = std::sqrt(alpha);
    double beta = alpha * rho;
    double delta = zeta / (alpha * std::sqrt(rho2));
    double mu = -beta * std::pow(delta, 2.0) * kappagh(zeta, lambda);
    
    return NumericVector::create(alpha, beta, delta, mu);
}

NumericVector paramghst(const double betabar, const double nu) {
    double delta = std::sqrt(1.0 / (((2.0 * betabar * betabar) / ((nu - 2.0) * (nu - 2.0) * (nu - 4.0))) + (1.0 / (nu - 2.0))));
    double beta = betabar / delta;
    double mu = -beta * (delta * delta) / (nu - 2.0);
    
    return NumericVector::create(nu, beta, delta, mu);
}


double dghst_std(const double x, const double betabar, const double nu) {
    NumericVector param = paramghst(betabar, nu);
    double beta = param[1];
    double delta = param[2];
    double mu = param[3];
    
    double betasqr = beta * beta;
    double deltasqr = delta * delta;
    double res = x - mu;
    
    double pdf = ((1.0 - nu) / 2.0) * std::log(2) + nu * std::log(delta) +
        ((nu + 1.0) / 2.0) * std::log(std::abs(beta)) +
        std::log(R::bessel_k(std::sqrt(betasqr * (deltasqr + res * res)), (nu + 1.0) / 2.0, 2.0)) -
        std::sqrt(betasqr * (deltasqr + res * res)) + beta * res -
        R::lgammafn(nu / 2.0) - std::log(M_PI) / 2.0 - ((nu + 1.0) / 2.0) * 
        std::log(deltasqr + res * res) / 2.0;
    
    pdf = std::exp(pdf);
    
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dghst(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dghst_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if(logr == 1) ans[i] = std::log(ans[i]);
    }
    
    return ans;
}

double rsghst_std(const double betabar, const double nu) {
    NumericVector param = paramghst(betabar, nu);
    double beta = param[1];
    double delta = param[2];
    double mu = param[3];
    
    // Sample y from the inverse gamma distribution
    double y = 1.0 / R::rgamma(nu / 2.0, 2.0 / (delta * delta));
    double sigma = std::sqrt(y);
    
    // Sample z from standard normal distribution
    double z = R::rnorm(0, 1);
    
    // Calculate the final result
    double ans = mu + beta * sigma * sigma + sigma * z;
    
    return ans;
}

// [[Rcpp::export]]
NumericVector c_rghst(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = mu[i] + rsghst_std(skew[i], shape[i]) * sigma[i];
    }
    
    return ans;
}

double dgh(const double x, const double alpha, const double beta, const double delta, const double mu, const double lambda) {
    if (alpha <= 0.0 || delta <= 0.0 || std::abs(beta) >= alpha) {
        return 0.0;
    }
    
    double arg = delta * std::sqrt(std::pow(alpha, 2.0) - std::pow(beta, 2.0));
    double a = (lambda / 2.0) * std::log(std::pow(alpha, 2.0) - std::pow(beta, 2.0)) -
        (std::log(std::sqrt(2.0 * M_PI)) + (lambda - 0.5) * std::log(alpha) +
        lambda * std::log(delta) + std::log(R::bessel_k(arg, lambda, 2.0)) - arg);
    
    double f = ((lambda - 0.5) / 2.0) * std::log(std::pow(delta, 2.0) + std::pow((x - mu), 2.0));
    arg = alpha * std::sqrt(std::pow(delta, 2.0) + std::pow((x - mu), 2.0));
    double k = std::log(R::bessel_k(arg, lambda - 0.5, 2.0)) - arg;
    double e = beta * (x - mu);
    
    return std::exp(a + f + k + e);
}

// [[Rcpp::export]]
NumericVector c_dghyp(NumericVector x, NumericVector alpha, NumericVector beta, NumericVector delta, NumericVector mu, NumericVector lambda, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dgh(x[i], alpha[i], beta[i], delta[i], mu[i], lambda[i]);
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double dgh_std(const double x, const double rho, const double zeta, const double lambda) {
    NumericVector param = paramgh(rho, zeta, lambda);
    double alpha = param[0];
    double beta = param[1];
    double delta = param[2];
    double mu = param[3];
    
    return dgh(x, alpha, beta, delta, mu, lambda);
}

// [[Rcpp::export]]
NumericVector c_dgh(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, NumericVector lambda, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dgh_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i], lambda[i]) / sigma[i];
        if(logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double dnig_std(const double x, const double rho, const double zeta) {
    double lambda = -0.5;
    NumericVector param = paramgh(rho, zeta, lambda);
    double alpha = param[0];
    double beta = param[1];
    double delta = param[2];
    double mu = param[3];
    
    double deltasq = delta * delta;
    double res = x - mu;
    
    // Calculate the PDF using the Normal Inverse Gaussian formula
    double pdf = -std::log(M_PI) + std::log(alpha) + std::log(delta) + 
        std::log(R::bessel_k(alpha * std::sqrt(deltasq + res * res), 1.0, 1)) +
        delta * std::sqrt(alpha * alpha - beta * beta) + 
        beta * res - 0.5 * std::log(deltasq + res * res);
    
    pdf = std::exp(pdf);
    
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dnig(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dnig_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double rstd_std(const double nu) {
    if (nu > 2.0) {
        double s = std::sqrt(nu / (nu - 2.0));
        return R::rt(nu) * 1.0 / s;
    }
    return 0.0; // Default to 0 if nu <= 2.0 (undefined)
}

// [[Rcpp::export]]
NumericVector c_rstd(int n, NumericVector mu, NumericVector sigma, NumericVector shape) {
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = mu[i] + rstd_std(shape[i]) * sigma[i];
    }
    
    return ans;
}

double xdt(const double x, const double nu) {
    double a = R::gammafn((nu + 1.0) / 2.0) / std::sqrt(M_PI * nu);
    double b = R::gammafn(nu / 2.0) * std::pow((1.0 + (x * x) / nu), (nu + 1.0) / 2.0);
    return a / b;
}

double dstd_std(const double x, const double nu) {
    if (nu <= 2.0) {
        return 999.0; // Undefined behavior for nu <= 2.0
    } else {
        double s = std::sqrt(nu / (nu - 2.0));
        return s * xdt(x * s, nu);
    }
}

// [[Rcpp::export]]
NumericVector c_dstd(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dstd_std((x[i] - mu[i]) / sigma[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double pstd_std(const double q, const double mu, const double sigma, const double nu) {
    double s = std::sqrt(nu / (nu - 2.0));
    double z = (q - mu) / sigma;
    double p = R::pt(z * s, nu, 1, 0);  // 1 for lower.tail = TRUE, 0 for log.p = FALSE
    return p;
}

// [[Rcpp::export]]
NumericVector c_pstd(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector shape) {
    int n = q.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = pstd_std(q[i], mu[i], sigma[i], shape[i]);
    }
    
    return ans;
}

double qstd_std(const double p, const double mu, const double sigma, const double nu) {
    double s = std::sqrt(nu / (nu - 2.0));
    double q = R::qt(p, nu, 1, 0) * sigma / s + mu;  // 1 for lower.tail = TRUE, 0 for log.p = FALSE
    return q;
}

// [[Rcpp::export]]
NumericVector c_qstd(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector shape) {
    int n = p.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = qstd_std(p[i], mu[i], sigma[i], shape[i]);
    }
    
    return ans;
}

double rsstd_std(const double xi, const double nu) {
    double weight = xi / (xi + 1.0 / xi);
    double z = R::runif(-1.0 * weight, 1.0 - weight);
    double xx = (z < 0) ? 1.0 / xi : xi;
    double rr = -1.0 * std::fabs(rstd_std(nu)) / xx * ((z < 0) ? -1.0 : 1.0);
    
    double a = 0.5, b = nu / 2.0;
    double beta_val = (R::gammafn(a) / R::gammafn(a + b)) * R::gammafn(b);
    
    double m1 = 2.0 * std::sqrt(nu - 2.0) / (nu - 1.0) / beta_val;
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    return (rr - mu) / sigma;
}

// [[Rcpp::export]]
NumericVector c_rsstd(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    NumericVector ans(n);
    for(int i = 0; i < n; i++) {
        ans[i] = mu[i] + rsstd_std(skew[i], shape[i]) * sigma[i];
    }
    return ans;
}

double dsstd_std(const double x, const double xi, const double nu) {
    double a = 0.5, b = nu / 2.0;
    double beta_val = (R::gammafn(a) / R::gammafn(a + b)) * R::gammafn(b);
    
    double m1 = 2.0 * std::sqrt(nu - 2.0) / (nu - 1.0) / beta_val;
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    double z = x * sigma + mu;
    double xxi = (z < 0) ? 1.0 / xi : xi;
    
    double g = 2.0 / (xi + 1.0 / xi);
    double pdf = g * dstd_std(z / xxi, nu) * sigma;
    
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dsstd(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = dsstd_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double psstd_std(const double q, const double mu, const double sigma, const double xi, const double nu) {
    double qx = (q - mu) / sigma;
    
    // Calculate m1, mu, and sigma based on xi and nu
    double m1 = 2.0 * std::sqrt(nu - 2.0) / (nu - 1.0) / (R::beta(0.5, nu / 2.0));
    double mux = m1 * (xi - 1.0 / xi);
    double sig = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    double z = qx * sig + mux;
    double Xi = (z < 0) ? 1.0 / xi : xi;
    double g = 2.0 / (xi + 1.0 / xi);
    
    // Compute the CDF using the heaviside and signum functions
    double p = heaviside(z, 0) - signum(z) * g * Xi * pstd_std(-std::fabs(z) / Xi, 0.0, 1.0, nu);
    return p;
}

// [[Rcpp::export]]
NumericVector c_psstd(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = q.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = psstd_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
    }
    
    return ans;
}

double qsstd_std(const double p, const double xi, const double nu) {
    double m1 = 2.0 * std::sqrt(nu - 2.0) / (nu - 1.0) / (R::beta(0.5, nu / 2.0));
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    double g = 2.0 / (xi + 1.0 / xi);
    
    double z = p - (1.0 / (1.0 + xi * xi));
    double Xi = std::pow(xi, signum(z));
    
    double tmp = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
    double q = (-signum(z) * qstd_std(tmp, 0.0, 1.0, nu) * Xi - mu) / sigma;
    
    return q;
}

// [[Rcpp::export]]
NumericVector c_qsstd(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = p.size();
    NumericVector ans(n);
    
    for(int i = 0; i < n; i++) {
        ans[i] = qsstd_std(p[i], skew[i], shape[i]) * sigma[i] + mu[i];
    }
    
    return ans;
}

double djsu_std(const double x, const double nu, const double tau) {
    double w, z, r, omega, c, pdf;
    double rtau = 1.0 / tau;
    
    if (rtau < 1e-7) {
        w = 1.0;
    } else {
        w = std::exp(rtau * rtau);
    }
    
    omega = -nu * rtau;
    c = std::sqrt(1.0 / (0.5 * (w - 1.0) * (w * std::cosh(2.0 * omega) + 1.0)));
    z = (x - (c * std::sqrt(w) * std::sinh(omega))) / c;
    r = -nu + std::asinh(z) / rtau;
    
    pdf = -std::log(c) - std::log(rtau) - 0.5 * std::log(z * z + 1.0) - 0.5 * std::log(2.0 * M_PI) - 0.5 * r * r;
    return std::exp(pdf);
}

// [[Rcpp::export]]
NumericVector c_djsu(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = djsu_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double qjsu_std(const double p, const double nu, const double tau) {
    double rtau = 1.0 / tau;
    double rr = R::qnorm(p, 0.0, 1.0, 1, 0);
    double z = std::sinh(rtau * (rr + nu));
    
    double w = (rtau < 1e-7) ? 1.0 : std::exp(rtau * rtau);
    double omega = -nu * rtau;
    double cc = std::sqrt(1.0 / (0.5 * (w - 1.0) * (w * std::cosh(2.0 * omega) + 1.0)));
    
    return (cc * std::sqrt(w) * std::sinh(omega)) + cc * z;
}

// [[Rcpp::export]]
NumericVector c_qjsu(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = p.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + qjsu_std(p[i], skew[i], shape[i]) * sigma[i];
    }
    
    return ans;
}

double pjsu_std(const double q, const double mu, const double sigma, const double nu, const double tau) {
    double rtau = 1.0 / tau;
    double w = (rtau < 1e-7) ? 1.0 : std::exp(rtau * rtau);
    double omega = -nu * rtau;
    double c = 1.0 / std::sqrt(0.5 * (w - 1.0) * (w * std::cosh(2.0 * omega) + 1.0));
    
    double z = (q - (mu + c * sigma * std::sqrt(w) * std::sinh(omega))) / (c * sigma);
    double r = -nu + std::asinh(z) / rtau;
    
    return R::pnorm(r, 0, 1, 1, 0);  // pnorm(r, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
}

// [[Rcpp::export]]
NumericVector c_pjsu(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = q.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = pjsu_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
    }
    
    return ans;
}


double rjsu_std(const double nu, const double tau) {
    double x = R::runif(0.0, 1.0);
    return qjsu_std(x, nu, tau);
}

// [[Rcpp::export]]
NumericVector c_rjsu(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    NumericVector ans(n);
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + rjsu_std(skew[i], shape[i]) * sigma[i];
    }
    return ans;
}

double rsnorm_std(const double xi) {
    double weight = xi / (xi + 1.0 / xi);
    double z = R::runif(-weight, 1.0 - weight);
    double xx = (z < 0) ? 1.0 / xi : xi;
    double rr = -1.0 * std::fabs(R::rnorm(0.0, 1.0)) / xx * ((z < 0) ? -1.0 : 1.0);
    
    double m1 = 2.0 / std::sqrt(2.0 * M_PI);
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * ((xi * xi) + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    return (rr - mu) / sigma;
}

// [[Rcpp::export]]
NumericVector c_rsnorm(int n, NumericVector mu, NumericVector sigma, NumericVector skew) {
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + rsnorm_std(skew[i]) * sigma[i];
    }
    
    return ans;
}

double dsnorm_std(const double x, const double xi) {
    double m1 = 2.0 / std::sqrt(2.0 * M_PI);
    double m12 = m1 * m1;
    double xi2 = xi * xi;
    
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m12) * (xi2 + 1.0 / xi2) + 2.0 * m12 - 1.0);
    
    double z = x * sigma + mu;
    double xxi = (z == 0.0) ? 1.0 : ((z < 0.0) ? 1.0 / xi : xi);
    
    double g = 2.0 / (xi + 1.0 / xi);
    double pdf = g * dnorm_std(z / xxi) * sigma;
    
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dsnorm(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = dsnorm_std((x[i] - mu[i]) / sigma[i], skew[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double psnorm_std(const double q, const double mu, const double sigma, const double xi) {
    double qx = (q - mu) / sigma;
    double m1 = 2.0 / std::sqrt(2.0 * M_PI);
    double mux = m1 * (xi - 1.0 / xi);
    double sig = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    double z = qx * sig + mux;
    double Xi = (z < 0) ? 1.0 / xi : xi;
    double g = 2.0 / (xi + 1.0 / xi);
    
    double p = heaviside(z, 0) - signum(z) * g * Xi * R::pnorm(-std::fabs(z) / Xi, 0.0, 1.0, 1, 0);
    return p;
}

// [[Rcpp::export]]
NumericVector c_psnorm(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew) {
    int n = q.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = psnorm_std(q[i], mu[i], sigma[i], skew[i]);
    }
    
    return ans;
}

double qsnorm_std(const double p, const double xi) {
    double m1 = 2.0 / std::sqrt(2.0 * M_PI);
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    double g = 2.0 / (xi + 1.0 / xi);
    double z = p - (1.0 / (1.0 + xi * xi));
    double Xi = std::pow(xi, signum(z));
    
    double tmp = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
    double q = (-signum(z) * R::qnorm(tmp, 0.0, Xi, 1, 0) - mu) / sigma;
    
    return q;
}

// [[Rcpp::export]]
NumericVector c_qsnorm(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew) {
    int n = p.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + qsnorm_std(p[i], skew[i]) * sigma[i];
    }
    
    return ans;
}

double rged_std(const double nu) {
    double lambda = std::sqrt(std::pow(0.5, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double rr = R::rgamma(1.0 / nu, 1.0);
    double ans = lambda * std::pow(2.0 * rr, 1.0 / nu) * ((R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0);
    return ans;
}

// [[Rcpp::export]]
NumericVector c_rged(int n, NumericVector mu, NumericVector sigma, NumericVector shape) {
    NumericVector ans(n);
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + rged_std(shape[i]) * sigma[i];
    }
    return ans;
}

double dged_std(const double x, const double nu) {
    double lambda = std::sqrt(std::pow(1.0 / 2.0, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double g = nu / (lambda * (std::pow(2.0, 1.0 + (1.0 / nu))) * R::gammafn(1.0 / nu));
    double pdf = g * std::exp(-0.5 * std::pow(std::fabs(x / lambda), nu));
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dged(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = dged_std((x[i] - mu[i]) / sigma[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double pged_std(const double q, const double mu, const double sigma, const double nu) {
    double qx = (q - mu) / sigma;
    double lambda = std::sqrt(1.0 / std::pow(2.0, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double g = nu / (lambda * (std::pow(2.0, 1.0 + 1.0 / nu)) * R::gammafn(1.0 / nu));
    double h = std::pow(2.0, 1.0 / nu) * lambda * g * R::gammafn(1.0 / nu) / nu;
    double s = 0.5 * std::pow(std::fabs(qx) / lambda, nu);
    double p = 0.5 + ((qx >= 0.0) ? 1.0 : -1.0) * h * R::pgamma(s, 1.0 / nu, 1.0, 1, 0);
    return p;
}

// [[Rcpp::export]]
NumericVector c_pged(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector shape) {
    int n = q.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = pged_std(q[i], mu[i], sigma[i], shape[i]);
    }
    
    return ans;
}

double qged_std(const double p, const double shape) {
    double y = 2.0 * p - 1.0;
    double lambda = std::sqrt(std::pow(0.5, 2.0 / shape) * R::gammafn(1.0 / shape) / R::gammafn(3.0 / shape));
    double q = lambda * std::pow(2.0 * R::qgamma(std::fabs(y), 1.0 / shape, 1.0, 1, 0), 1.0 / shape);
    q = q * ((y >= 0) ? 1.0 : -1.0);
    return q;
}

// [[Rcpp::export]]
NumericVector c_qged(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector shape) {
    int n = p.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = qged_std(p[i], shape[i]) * sigma[i] + mu[i];
    }
    
    return ans;
}

double rsged_std(const double xi, const double nu) {
    double weight = xi / (xi + 1.0 / xi);
    double z = R::runif(-weight, 1.0 - weight);
    double xx = (z < 0) ? 1.0 / xi : xi;
    double rr = -1.0 * std::fabs(rged_std(nu)) / xx * ((z < 0) ? -1.0 : 1.0);
    
    double lambda = std::sqrt(std::pow(0.5, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double m1 = std::pow(2.0, 1.0 / nu) * lambda * R::gammafn(2.0 / nu) / R::gammafn(1.0 / nu);
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    return (rr - mu) / sigma;
}

// [[Rcpp::export]]
NumericVector c_rsged(int n, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = mu[i] + rsged_std(skew[i], shape[i]) * sigma[i];
    }
    
    return ans;
}

double dsged_std(const double x, const double xi, const double nu) {
    double lambda = std::sqrt(std::pow(1.0 / 2.0, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double m1 = std::pow(2.0, 1.0 / nu) * lambda * R::gammafn(2.0 / nu) / R::gammafn(1.0 / nu);
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    
    double z = x * sigma + mu;
    double xxi = (z == 0) ? 1.0 : ((z < 0) ? 1.0 / xi : xi);
    
    double g = 2.0 / (xi + 1.0 / xi);
    double pdf = g * dged_std(z / xxi, nu) * sigma;
    
    return pdf;
}

// [[Rcpp::export]]
NumericVector c_dsged(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = dsged_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}


double psged_std(const double q, const double mu, const double sigma, const double xi, const double nu) {
    double qx = (q - mu) / sigma;
    double lambda = std::sqrt(1.0 / std::pow(2.0, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double m1 = std::pow(2.0, 1.0 / nu) * lambda * R::gammafn(2.0 / nu) / R::gammafn(1.0 / nu);
    double mux = m1 * (xi - 1.0 / xi);
    double sig = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    double z = qx * sig + mux;
    double Xi = (z < 0) ? 1.0 / xi : xi;
    double g = 2.0 / (xi + 1.0 / xi);
    double p = heaviside(z, 0) - signum(z) * g * Xi * pged_std(-std::fabs(z) / Xi, 0.0, 1.0, nu);
    return p;
}

// [[Rcpp::export]]
NumericVector c_psged(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = q.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = psged_std(q[i], mu[i], sigma[i], skew[i], shape[i]);
    }
    
    return ans;
}

double qsged_std(const double p, const double xi, const double nu) {
    double lambda = std::sqrt(1.0 / std::pow(2.0, 2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
    double m1 = std::pow(2.0, 1.0 / nu) * lambda * R::gammafn(2.0 / nu) / R::gammafn(1.0 / nu);
    double mu = m1 * (xi - 1.0 / xi);
    double sigma = std::sqrt((1.0 - m1 * m1) * (xi * xi + 1.0 / (xi * xi)) + 2.0 * m1 * m1 - 1.0);
    double g = 2.0 / (xi + 1.0 / xi);
    double z = p - (1.0 / (1.0 + xi * xi));
    double Xi = std::pow(xi, signum(z));
    double q = (heaviside(z, 0) - signum(z) * p) / (g * Xi);
    q = (-signum(z) * qged_std(q, nu) * Xi - mu) / sigma;
    
    return q;
}

// [[Rcpp::export]]
NumericVector c_qsged(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape) {
    int n = p.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = qsged_std(p[i], skew[i], shape[i]) * sigma[i] + mu[i];
    }
    
    return ans;
}


double dhyp(const double x, const double alpha, const double beta, const double delta, const double mu) {
    double pdf = 0.0;
    if (alpha <= 0.0 || delta <= 0 || std::fabs(beta) >= alpha) {
        return pdf;
    }
    
    double g = alpha * alpha - beta * beta;
    double e = x - mu;
    double besselK = R::bessel_k(delta * std::sqrt(g), 1.0, 2);
    
    if (besselK == 0.0) return 0.0;
    
    pdf = 0.5 * std::log(g) - std::log(2.0 * alpha * delta * besselK) - alpha * std::sqrt(delta * delta + e * e) + beta * e;
    return std::exp(pdf);
}

double dhyp_std(const double x, const double rho, const double zeta) {
    NumericVector param = paramgh(rho, zeta, 1.0);
    double alpha = param[0];
    double beta = param[1];
    double delta = param[2];
    double mu = param[3];
    
    return dhyp(x, alpha, beta, delta, mu);
}

// [[Rcpp::export]]
NumericVector c_dhyp(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector skew, NumericVector shape, int logr) {
    int n = x.size();
    NumericVector ans(n);
    
    for (int i = 0; i < n; i++) {
        ans[i] = dhyp_std((x[i] - mu[i]) / sigma[i], skew[i], shape[i]) / sigma[i];
        if (logr == 1) {
            ans[i] = std::log(ans[i]);
        }
    }
    
    return ans;
}

double g_function(double y, double beta, double m, double lambda) {
    return 0.5 * beta * std::pow(y, 3) - std::pow(y, 2) * (0.5 * beta * m + lambda + 1) + y * (-0.5 * beta) + 0.5 * beta * m;
}

double f_function(double y, double beta, double m, double lambda) {
    return 0.5 * beta * std::pow(y, 3) - std::pow(y, 2) * (0.5 * beta * m + lambda + 1) +
        y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m;
}

NumericVector rgig1(int n, NumericVector param) {
    double chi_val = param[0];
    double psi_val = param[1];
    double lambda = 1.0;
    
    double alpha = std::sqrt(psi_val / chi_val);
    double beta = std::sqrt(psi_val * chi_val);
    double m = std::abs(beta) / beta;
    
    Rcpp::Function uniroot("uniroot");
    
    List yM_result = uniroot(Named("f") = Rcpp::InternalFunction(&g_function),
                             Named("interval") = NumericVector::create(0.0, m),
                             Named("beta") = beta,
                             Named("m") = m,
                             Named("lambda") = lambda);
    
    double lower2 = m;
    double upper2 = 2 * m;
    
    while (g_function(lower2, beta, m, lambda) * g_function(upper2, beta, m, lambda) > 0) {
        upper2 *= 2;
        if (upper2 > 1000) {
            stop("No root found for yP within a reasonable interval.");
        }
    }
    
    List yP_result = uniroot(Named("f") = Rcpp::InternalFunction(&g_function),
                             Named("interval") = NumericVector::create(lower2, upper2),
                             Named("beta") = beta,
                             Named("m") = m,
                             Named("lambda") = lambda);
    
    double yM = as<double>(yM_result["root"]);
    double yP = as<double>(yP_result["root"]);
    
    double a = (yP - m) * std::exp(-0.25 * beta * (yP + 1 / yP - m - 1 / m));
    double b = (yM - m) * std::exp(-0.25 * beta * (yM + 1 / yM - m - 1 / m));
    double c = -0.25 * beta * (m + 1 / m);
    
    NumericVector output(n);
    
    for (int i = 0; i < n; i++) {
        bool needValue = true;
        
        while (needValue) {
            double R1 = R::runif(0.0, 1.0);
            double R2 = R::runif(0.0, 1.0);
            double Y = m + a * R2 / R1 + b * (1 - R2) / R1;
            
            if (Y > 0) {
                if (-std::log(R1) >= 0.25 * beta * (Y + 1 / Y) + c) {
                    needValue = false;
                    output[i] = Y;
                }
            }
        }
    }
    
    return output / alpha;
}


NumericVector rgig(int n, NumericVector param) {
    if (param.size() != 3) {
        stop("Parameter vector must have exactly 3 elements.");
    }
    double chi = param[0];
    double psi = param[1];
    double lambda = param[2];
    double alpha = std::sqrt(psi / chi);
    double beta = std::sqrt(psi * chi);
    double m = (lambda - 1 + std::sqrt((lambda - 1) * (lambda - 1) + beta * beta)) / beta;
    Rcpp::Function uniroot("uniroot");
    double upper = m;
    while (f_function(upper, beta, m, lambda) <= 0) {
        upper *= 2;
    }
    List yM_result = uniroot(Named("f") = Rcpp::InternalFunction(&f_function),
                             Named("interval") = NumericVector::create(0.0, m),
                             Named("beta") = beta,
                             Named("m") = m,
                             Named("lambda") = lambda);
    
    List yP_result = uniroot(Named("f") = Rcpp::InternalFunction(&f_function),
                             Named("interval") = NumericVector::create(m, upper),
                             Named("beta") = beta,
                             Named("m") = m,
                             Named("lambda") = lambda);
    
    double yM = as<double>(yM_result["root"]);
    double yP = as<double>(yP_result["root"]);
    
    double a = (yP - m) * std::pow(yP / m, 0.5 * (lambda - 1)) * std::exp(-0.25 * beta * (yP + 1 / yP - m - 1 / m));
    double b = (yM - m) * std::pow(yM / m, 0.5 * (lambda - 1)) * std::exp(-0.25 * beta * (yM + 1 / yM - m - 1 / m));
    double c = -0.25 * beta * (m + 1 / m) + 0.5 * (lambda - 1) * std::log(m);
    
    NumericVector output(n);
    
    for (int i = 0; i < n; i++) {
        bool needValue = true;
        
        while (needValue) {
            double R1 = R::runif(0.0, 1.0);
            double R2 = R::runif(0.0, 1.0);
            double Y = m + a * R2 / R1 + b * (1 - R2) / R1;
            
            if (Y > 0) {
                if (-std::log(R1) >= -0.5 * (lambda - 1) * std::log(Y) + 0.25 * beta * (Y + 1 / Y) + c) {
                    needValue = false;
                    output[i] = Y;
                }
            }
        }
    }
    
    return output / alpha;
}

// [[Rcpp::export]]
NumericVector c_rghyp(int n, double mu = 0, double delta = 1, double alpha = 1, double beta = 0, double lambda = 1) {
    double chi = delta * delta;
    double psi = alpha * alpha - beta * beta;
    NumericVector X;
    if (lambda == 1) {
        X = rgig1(n, NumericVector::create(chi, psi));
    } else {
        X = rgig(n, NumericVector::create(chi, psi, lambda));
    }
    NumericVector sigma = sqrt(X);
    NumericVector Z = rnorm(n);
    NumericVector Y = mu + beta * pow(sigma, 2) + sigma * Z;
    return Y;
}
