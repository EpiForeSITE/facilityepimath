## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following a CRAN check failure with the MKL
alternative BLAS/LAPACK implementation.

### Changes

* Changed numerical root-finding in `facilityeq()` from minimizing a
  squared residual with `optimize()` to solving the equilibrium equation
  directly with `uniroot()`. This should resolve numerical differences in
  `facilityeq()` under alternative BLAS/LAPACK implementations.
* Corrected the use of `eigM$value` to `eigM$values` in `facilityR0()`.
