#' Calculate the mean length of stay for a linear facility model
#'
#' @param M A matrix of state transition rates between facility patient states
#' @param init A vector of admission state probabilities to each state
#' @param mgf The moment generating function characterizing a time-of-stay-dependent removal hazard
#' @return The mean length of stay
#' @examples
#' M <- rbind(c(-1.1,2),c(1,-2.2))
#' init <- c(0.9,0.1)
#' mgf <- function(x, deriv=0) MGFgamma(x, rate=0.2, shape=3, deriv)
#' meanlengthofstay(M, init, mgf)
#' @export
meanlengthofstay <- function(M, init, mgf){
  if(is.null(mgf)){
    return(sum(as.vector(solve(M,init))))
  }
  if(all(abs(colSums(M))<.Machine$double.eps)){
    return(mgf(0,1))
  }
  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
  K <- Vectorize(K,'x')
  eig <- eigen(M)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(eig$values)
  sum(as.vector(eig$vectors %*% (solveV * Klamb)))
}
