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
  K <- NULL
  if(!is.null(mgf)){
    K <- function(x, deriv = 0)
      ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
    K <- Vectorize(K,'x')
  }

  mfun <- function(alpha) rbind(cbind(Sfun(alpha),R),cbind(Afun(alpha),C))
  colinds <- (nrow(as.matrix(R))+1):(nrow(as.matrix(R))+nrow(as.matrix(C)))

  getbeta <- function(alpha){
    eq <- equilib(mfun(alpha),init,K)
    alpha/sum(eq[colinds]*transm)
  }

  getalpha <- function(beta){
    err <- function(x) (getbeta(x) - beta)^2
    maxalpha <- 1
    while(getbeta(maxalpha) < beta) maxalpha <- maxalpha*10
    optimize(f = err, interval = c(0,maxalpha), tol=1e-10)$minimum
  }

  alpha <- getalpha(1)
  equilib(mfun(alpha),init,K)
}
