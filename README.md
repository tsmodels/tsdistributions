
[![R-CMD-check](https://github.com/tsmodels/tsdistributions/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tsmodels/tsdistributions/actions/workflows/R-CMD-check.yaml)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--05--17-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-0.2.0-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsdistributions)](https://cran.r-project.org/package=tsdistributions)

# tsdistributions

A number of location-scale invariant distributions with full d,p,q,r
methods and autodiff (TMB) backed estimation. These distributions are
parameterized in terms of the mean, standard deviation, skew and shape,
with the GHYP distribution having 2 shape parameters.

The currently implemented distributions are shown below.

| Distribution                        | Function | Parameters | Notes             |
|:------------------------------------|:--------:|-----------:|------------------:|
| Normal                              |   norm   |          2 |                   |
| Student                             |   std    |          3 |                   |
| GED                                 |   ged    |          3 |                   |
| Skew-Normal                         |  snorm   |          3 | Fernandez & Steel |
| Skew-Student                        |   sstd   |          4 | Fernandez & Steel |
| Skew-GED                            |   sged   |          4 | Fernandez & Steel |
| Normal Inverse Gaussian             |   nig    |          4 | Barndorff-Nielsen |
| Generalized Hyperbolic              |   ghyp   |          5 | Barndorff-Nielsen |
| Johnson’s SU                        |   jsu    |          4 | Johnson           |
| Generalized Hyperbolic Skew Student |   ghst   |          4 | Aas & Haff        |
