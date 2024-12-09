#' Calculate the equilibrium of a facility transmission model
#'
#' @param Sfun A function that produces a matrix of state transition rates between and removal from the susceptible states for a given value of the acquisition rate
#' @param C A matrix of state transition rates between and removal from the colonized states
#' @param Afun A function that produces a matrix of transition rates from susceptible to colonized states for a given value of the acquisition rate
#' @param R A matrix of recovery rates: state transition rates from colonized to susceptible states
#' @param transm A vector of transmission rates from each colonized state
#' @param init A vector of admission state probabilities to each state
#' @param mgf The moment generating function characterizing a time-of-stay-dependent removal hazard
#' @importFrom stats optimize
#' @return A vector with the proportion of patients in each state at equilibrium
#' @export
facilityeq <- function(Sfun,C,Afun,R,transm,init,mgf=NULL){

  mfun <- function(alpha) rbind(cbind(Sfun(alpha),R),cbind(Afun(alpha),C))
  colinds <- (nrow(as.matrix(R))+1):(nrow(as.matrix(R))+nrow(as.matrix(C)))

  getbeta <- function(alpha){
    eq <- equilib(mfun(alpha),init,mgf)
    alpha/sum(eq[colinds]*transm)
  }

  maxalpha <- 1
  while(getbeta(maxalpha) < 1) maxalpha <- maxalpha*10
  alpha <- optimize(f = function(x) (getbeta(x) - 1)^2, interval = c(0,maxalpha), tol=1e-10)$minimum

  equilib(mfun(alpha),init,mgf)
}
