#' Calculate the equilibrium of a linear facility model
#'
#' @param M A matrix of state transition rates between facility patient states
#' @param init A vector of admission state probabilities to each state
#' @param mgf The moment generating function characterizing a time-of-stay-dependent removal hazard
#' @return A vector with the proportion of patients in each state at equilibrium
#' @examples
#' M <- rbind(c(-0.06,0.03,0),c(0.06,-0.08,0),c(0,0.05,0))
#' init <- c(0.95,0.05,0)
#' mgf <- function(x, deriv=0) MGFgamma(x, rate = 0.05, shape = 2.5, deriv)
#' equilib(M, init, mgf)
#' @export
equilib <- function(M, init, mgf=NULL){
  if(is.null(mgf)){
    MinvInit <- as.vector(solve(M,init))
    return(MinvInit/sum(MinvInit))
  }
  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
  K <- Vectorize(K,'x')
  eig <- eigen(M)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(sapply(eig$values,function(x) min(x,0)))
  numK <- as.vector(eig$vectors %*% (solveV * Klamb))
  numK/sum(numK)
}
