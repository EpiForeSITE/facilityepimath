#' Calculate the equilibrium of a facility transmission model
#'
#' @param S A matrix of state transition rates between and removal from the susceptible states in the absence of colonized individuals
#' @param C A matrix of state transition rates between and removal from the colonized states
#' @param A A matrix describing transitions from susceptible to colonized states at acquisition
#' @param R A matrix of recovery rates: state transition rates from colonized to susceptible states
#' @param transm A vector of transmission rates from each colonized state
#' @param init A vector of admission state probabilities to each state
#' @param mgf The moment generating function characterizing a time-of-stay-dependent removal hazard
#' @importFrom stats optimize
#' @return A vector with the proportion of patients in each state at equilibrium
#' @examples
#' S <- 0
#' C <- rbind(c(-0.38,0),c(0.08,0))
#' A <- rbind(1,0)
#' R <- cbind(0.3,0)
#' transm <- c(0.1,0.05)
#' init <- c(0.99,0.01,0)
#' mgf <- function(x, deriv=0) MGFgamma(x, rate=0.2, shape=3, deriv)
#' facilityeq(S, C, A, R, transm, init, mgf)
#'
#' @export
facilityeq <- function(S,C,A,R,transm,init,mgf=NULL){

  Sfun <- function(x) S-diag(colSums(as.matrix(A)))*x
  Afun <- function(x) A*x

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
