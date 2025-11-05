

<!-- README.md is generated from README.Rmd. Please edit that file -->

# facilityepimath

[![CRAN
status](https://www.r-pkg.org/badges/version/facilityepimath)](https://CRAN.R-project.org/package=facilityepimath)
[![CRANlogs
downloads](https://cranlogs.r-pkg.org/badges/facilityepimath)](https://cran.r-project.org/package=facilityepimath)
[![ForeSITE Group](https://github.com/EpiForeSITE/software/raw/e82ed88f75e0fe5c0a1a3b38c2b94509f122019c/docs/assets/foresite-software-badge.svg)](https://github.com/EpiForeSITE)

The goal of facilityepimath is to provide functions to calculate useful
quantities for a user-defined differential equation model of infectious
disease transmission among individuals in a healthcare facility,
including the basic facility reproduction number and model equilibrium.
A full description and derivation of the mathematical results
implemented in these functions can be found in the following manuscript:

Toth D, Khader K, Mitchell M, Samore M (2025). Transmission thresholds
for the spread of infections in healthcare facilities. PLoS
Computational Biology 21(10): e1013577.
https://doi.org/10.1371/journal.pcbi.1013577.

This work was supported by the Centers for Disease Control and
Prevention, Modeling Infectious Diseases in Healthcare Network award
U01CK000585 and Insight Net award number CDC-RFA-FT-23-0069.

## Installation

You can install the development version of facilityepimath from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("EpiForeSITE/facilityepimath")
```

You can install the facilitymath from CRAN with:

``` r
install.packages("facilityepimath")
```

## Example

The functions in this package analyze compartmental differential
equation models, where the compartments are states of individuals who
are co-located in a healthcare facility, for example, patients in a
hospital or patients/residents of long-term care facility.

Our example here is a hospital model with four compartments: two
compartments $C_1$ and $C_2$ representing patients who are colonized
with an infectious organism, and two compartments $S_1$ and $S_2$
representing patients who are susceptible to acquiring colonization with
that same organism. The two colonized and susceptible states are
distinguished by having potentially different infectivity to other
patients and vulnerability to acquisition from other patients,
respectively.

A system of differential equations with those 4 compartments may take
the following general form:

\[ = -(s\_{21}+(a\_{11}+a\_{21})+*1+h(t))S_1 + s*{12}S_2 + r\_{11}C_1 +
r\_{12}C_2\]

\[ = s\_{21}S_1 - (s\_{12}+(a\_{12}+a\_{22})+*2+h(t))S_2 + r*{21}C_1 +
r\_{22}C_2 \]

\[ = a\_{11}S_1 + a\_{12}S_2 - (c\_{21}+r\_{11}+r\_{21}+*3+h(t))C_1 +
c*{12}C_2 \]

\[ = a\_{21}S_1 + a\_{22}S_2 + c\_{21}C_1 -
(c\_{12}+r\_{12}+r\_{22}+\_4+h(t))C_2 \]

The acquisition rate $\alpha$ appearing in each equation, and governing
the transition rates between the S compartments and the C compartments,
is assumed to depend on the number of colonized patients in the
facility, as follows:

\[ = \_1 C_1 + \_2 C_2 \]

We will demonstrate how to calculate the basic reproduction number $R_0$
of this system using the `facilityR0` function. The following components
of the system are required as inputs to the function call below.

A matrix `S` governing the transitions between, and out of, the states
$S_1$ and $S_2$ in the absence of any colonized patients:

\[S = (
) \]

A matrix `C` governing the transitions between, and out of, the states
$C_1$ and $C_2$:

\[C = (
) \]

A matrix `A` describing the S-to-C state transitions when an acquisition
occurs:

\[A = (
) \]

A vector `transm` containing the $\beta$ coefficients (transmission
rates from each colonized compartment) appearing in the $\alpha$
equation: $(\beta_1,\beta_2)$

A vector `initS` containing the admission state probabilities for the
susceptible compartments only (i.e., the pre-invasion system before a
colonized patient is introduced): $(\theta_1,1-\theta_1)$

A function `mgf(x,deriv)` that is the moment-generating function (and
its derivatives) of the distribution for which the
time-of-stay-dependent removal rate `h(t)` is the hazard function. This
is the length of stay distribution when the state-dependent removal
rates $\omega$ are zero.

``` r
library(facilityepimath)
S <- rbind(c(-1,2),c(1,-2))
C <- rbind(c(-1.1,0),c(0.1,-0.9))
A <- rbind(c(1,0),c(0,2))
transm <- c(0.4,0.6)
initS <- c(0.9,0.1)

mgf <- function(x, deriv=0) MGFgamma(x, rate=0.01, shape=3.1, deriv)
facilityR0(S,C,A,transm,initS,mgf)
#> [1] 0.7244774
```
