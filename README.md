
[![R-CMD-check](https://github.com/tsmodels/tsdistributions/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tsmodels/tsdistributions/actions/workflows/R-CMD-check.yaml)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--08--22-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.2-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsdistributions)](https://cran.r-project.org/package=tsdistributions)

# tsdistributions

A number of location-scale invariant distributions with full d,p,q,r
methods and autodiff (TMB) backed estimation. These distributions are
parameterized in terms of the mean, standard deviation, skew and shape,
with the GH distribution having 2 shape parameters.

Currently implemented location-scale family distributions:

| Distribution                        | Function | Parameters | Notes             |
|:------------------------------------|:--------:|-----------:|:------------------|
| Normal                              |   norm   |          2 |                   |
| Student                             |   std    |          3 |                   |
| GED                                 |   ged    |          3 |                   |
| Skew-Normal                         |  snorm   |          3 | Fernandez & Steel |
| Skew-Student                        |   sstd   |          4 | Fernandez & Steel |
| Skew-GED                            |   sged   |          4 | Fernandez & Steel |
| Normal Inverse Gaussian             |   nig    |          4 | Barndorff-Nielsen |
| Generalized Hyperbolic              |    gh    |          5 | Barndorff-Nielsen |
| Johnson’s SU                        |   jsu    |          4 | Johnson           |
| Generalized Hyperbolic Skew Student |   ghst   |          4 | Aas & Haff        |

A number of other distributions are also implemented or planned to be,
which do not technically belong to the location-scale family but are
useful in time-series and financial risk analysis.

The semi-parametric distribution (`spd`) models the upper and lower
tails using the Generalized Pareto distribution. This leads to 1 scale
and 1 shape parameters for each of the tails (4 parameters in total).

Other distributions:

| Distribution                 | Function | Parameters | Notes   |
|:-----------------------------|:--------:|-----------:|:--------|
| Semi Parametric Distribution |   spd    |          4 | Carmona |
