#' Calculate the equilibrium of a linear facility model
#'
#' @param M A matrix of state transition rates between facility patient states
#' @param init A vector of admission state probabilities to each state
#' @param mgf The moment generating function characterizing a time-of-stay-dependent removal hazard
#' @return A vector with the proportion of patients in each state at equilibrium
#' @export
equilib <- function(M, init, mgf){
  if(is.null(mgf)){
    MinvInit <- as.vector(solve(M,init))
    return(MinvInit/sum(MinvInit))
  }
  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
  K <- Vectorize(K,'x')
  eig <- eigen(M)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(eig$values)
  numK <- as.vector(eig$vectors %*% (solveV * Klamb))
  numK/sum(numK)
}
