## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following a remaining test failure with the MKL
alternative BLAS/LAPACK implementation.

### Changes

* Changed equilibrium root calculations from minimizing squared residuals
  with `optimize()` to solving the corresponding equations directly with
  `uniroot()`, including both the package implementation and test
  calculations.
* Adjusted the tolerance for one test comparing independently calculated
  equilibrium roots from `sqrt(.Machine$double.eps)` to `1e-5` to account
  for numerical differences across BLAS/LAPACK implementations.
