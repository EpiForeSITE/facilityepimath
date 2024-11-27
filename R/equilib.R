equilib <- function(W, init, K){
  if(is.null(K)){
    WinvInit <- as.vector(solve(W,init))
    return(WinvInit/sum(WinvInit))
  }
  eig <- eigen(W)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(eig$values)
  numK <- as.vector(eig$vectors %*% (solveV * Klamb))
  numK/sum(numK)
}
