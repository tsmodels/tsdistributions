# tsdistributions 1.0.0

* Initial CRAN submission.
* Cleaned up c-code naming of distributions and removed redundant code.
* Removed gig.c and nig.c which were used for generating random numbers. Now
using imported functionality (`rghyp`) from GeneralizedHyperbolic package with 
transformed parameters.
* Cleaned up documentation and naming convention for file names with methods.
* Renamed Generalized Hyperbolic acronym to 'gh' (from 'ghyp').
* Added a `tsprofile` method.
* Added an html estimation example vignette.
* Added an html profile example vignette.
* Added a pdf vignette for location scale distributions outlining the
formulation of each distribution.

# tsdistributions 0.3.3

* Fix to the Generalized Hyperbolic Distribution estimation code as well as the 
scaled besselK function in the presence of fixed lambda. The fix somewhat negates 
the speedup from the previous version when lambda is fixed, but ensures correctness. 

# tsdistributions 0.3.2

* The Jacobian, which is used in the calculation of the QMLE and OPG covariance 
matrix is now calculated using TMB rather than numDeriv since that had problems 
for the ghst and gh distribution giving hugely erroneous results. The method of 
calculating it loops through each data point, and is rather hacky (but works). 
Using autodiff::jacobian directly in the TMB C++ code is proving challenging.
* Removed numDeriv dependency

# tsdistributions 0.3.1

* Fix to Fernandez and Steel skewed distributions quantile function (qsnorm, qsged and qsstd). 
* Added additional unit tests.

# tsdistributions 0.2.0

* Fix to snorm distribution (dsnorm) for cases when variable z == 0. Other skewed 
distributions used signum function but this distribution seems to have missed the 
0 case for some reason.
* Added vcov methods and other sandwich estimators (QMLE, OP, NW).
* Added BIC and AIC and created a print method for the summary.
* Fixed some code fomatting issues
