equilib <- function(W, init, mgf){
  if(is.null(mgf)){
    WinvInit <- as.vector(solve(W,init))
    return(WinvInit/sum(WinvInit))
  }
  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
  K <- Vectorize(K,'x')
  eig <- eigen(W)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(eig$values)
  numK <- as.vector(eig$vectors %*% (solveV * Klamb))
  numK/sum(numK)
}
