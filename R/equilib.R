equilib <- function(W, init, K){
  eig <- eigen(W)
  solveV <- solve(eig$vectors,init)
  Klamb <- K(eig$values)
  numK <- as.vector(eig$vectors %*% (solveV * Klamb))
  numK/sum(numK)
}
