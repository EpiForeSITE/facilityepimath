% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facilityeq.R
\name{facilityeq}
\alias{facilityeq}
\title{Calculate the equilibrium of a facility transmission model}
\usage{
facilityeq(S, C, A, R, transm, init, mgf = NULL)
}
\arguments{
\item{S}{A matrix of state transition rates between and removal from the susceptible states in the absence of colonized individuals}

\item{C}{A matrix of state transition rates between and removal from the colonized states}

\item{A}{A matrix describing transitions from susceptible to colonized states at acquisition}

\item{R}{A matrix of recovery rates: state transition rates from colonized to susceptible states}

\item{transm}{A vector of transmission rates from each colonized state}

\item{init}{A vector of admission state probabilities to each state}

\item{mgf}{The moment generating function characterizing a time-of-stay-dependent removal hazard}
}
\value{
A vector with the proportion of patients in each state at equilibrium; the vector contains the equilibrium S states followed by C states
}
\description{
Calculate the equilibrium of a facility transmission model
}
\examples{
S <- 0
C <- rbind(c(-0.38,0),c(0.08,0))
A <- rbind(1,0)
R <- cbind(0.3,0)
transm <- c(0.1,0.05)
init <- c(0.99,0.01,0)
mgf <- function(x, deriv=0) MGFgamma(x, rate=0.2, shape=3, deriv)
facilityeq(S, C, A, R, transm, init, mgf)

}
