
<!-- README.md is generated from README.Rmd. Please edit that file -->

# facilityepimath

<!-- badges: start -->

[![ForeSITE
Group](https://github.com/EpiForeSITE/software/blob/e82ed88f75e0fe5c0a1a3b38c2b94509f122019c/docs/assets/foresite-software-badge.svg)](https://github.com/EpiForeSITE)
<!-- badges: end -->

The goal of facilityepimath is to provide functions to calculate useful
quantities for a user-defined differential equation model of infectious
disease transmission among individuals in a healthcare facility,
including the basic facility reproduction number and model equilibrium.
A full description and derivation of the mathematical results
implemented in these functions can be found in the following manuscript:

Toth D, Khader K, Mitchell M, Samore M (2025). Transmission thresholds
for the spread of infections in healthcare facilities.
<https://doi.org/10.1101/2025.02.21.25322698>.

This work was supported by the Centers for Disease Control and
Prevention, Modeling Infectious Diseases in Healthcare Network award
U01CK000585.

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

$$ \frac{dS_1}{dt} = -(s_{21}+(a_{11}+a_{21})\alpha+\omega_1+h(t))S_1 + s_{12}S_2 + r_{11}C_1 + r_{12}C_2$$

$$\frac{dS_2}{dt} = s_{21}S_1 - (s_{12}+(a_{12}+a_{22})\alpha+\omega_2+h(t))S_2 + r_{21}C_1 + r_{22}C_2 $$

$$\frac{dC_1}{dt} = a_{11}\alpha S_1 + a_{12}\alpha S_2 - (c_{21}+r_{11}+r_{21}+\omega_3+h(t))C_1 + c_{12}C_2 $$

$$\frac{dC_2}{dt} = a_{21}\alpha S_1 + a_{22}\alpha S_2 + c_{21}C_1 - (c_{12}+r_{12}+r_{22}+\omega_4+h(t))C_2 $$

The acquisition rate $\alpha$ appearing in each equation, and governing
the transition rates between the S compartments and the C compartments,
is assumed to depend on the number of colonized patients in the
facility, as follows:

$$ \alpha = \beta_1 C_1 + \beta_2 C_2 $$

We will demonstrate how to calculate the basic reproduction number $R_0$
of this system using the `facilityR0` function. The following components
of the system are required as inputs to the function call below.

A matrix `S` governing the transitions between, and out of, the states
$S_1$ and $S_2$ in the absence of any colonized patients:

$$S = \left(
\begin{matrix}
    -s_{21}-\omega_1 & s_{12} \\
    s_{21} & -s_{12}-\omega_2
\end{matrix}\right)
$$

A matrix `C` governing the transitions between, and out of, the states
$C_1$ and $C_2$:

$$C = \left(
\begin{matrix}
    - (c_{21}+r_{11}+r_{21}+\omega_3) & c_{12} \\
    c_{21} & - (c_{12}+r_{12}+r_{22}+\omega_4)
\end{matrix}\right)
$$

A matrix `A` describing the S-to-C state transitions when an acquisition
occurs:

$$A = \left(
\begin{matrix}
    a_{11} & a_{12} \\
    a_{21} & a_{22}
\end{matrix}\right)
$$

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
