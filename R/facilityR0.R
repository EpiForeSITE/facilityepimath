#' Calculate basic reproduction number R0
#'
#' @param S A matrix of state transition rates between and removal from the susceptible states in the absence of colonized individuals
#' @param C A matrix of state transition rates between and removal from the colonized states
#' @param A A matrix describing transitions from susceptible to colonized states at acquisition
#' @param transm A vector of transmission rates from each colonized state
#' @param initS A vector of admission state probabilities to each susceptible state
#' @param mgf The moment generating function characterizing the time-of-stay-dependent removal hazard
#'
#' @importFrom MASS ginv
#'
#' @return A number (R0)
#' @export
facilityR0 <- function(S,C,A,transm,initS,mgf=NULL){

  Smat <- as.matrix(S)
  impS <- (initS > 0)
  initSadj <- initS[impS]
  Sadj <- as.matrix(S)[impS,impS]
  Aadj <- as.matrix(as.matrix(A)[,impS])

  if(is.null(mgf)){
    Ssol <- solve(Sadj, initSadj)
    R0 <- -sum(transm * solve(C, Aadj %*% Ssol)) / sum(Ssol)
    return(R0)
  }

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))
  K <- Vectorize(K,'x')

  n <- nrow(as.matrix(Sadj))
  m <- nrow(as.matrix(C))
  M <- rbind(cbind(Sadj,matrix(0,n,m)),cbind(Aadj,C))
  eigM <- eigen(M, symmetric = FALSE)
  dup <- duplicated(eigM$values)
  PcolFilled <- !dup
  if(sum(dup)>0){
    P <- eigM$vectors
    KJ <- diag(K(eigM$values))
    for(lambrep in unique(eigM$values[dup])){
      Mlamb <- M-diag(lambrep,n+m)
      inds <- which(eigM$value == lambrep)
      uniqueEigVecs <- as.matrix(P[,inds[1]])
      uniqueEigVecInds <- inds[1]
      for(i in 2:length(inds)){
        if(all(apply(uniqueEigVecs, 2, function(v) !equaleigvec(v,P[,inds[i]])))){
          uniqueEigVecs <- cbind(uniqueEigVecs, P[,inds[i]])
          uniqueEigVecInds <- c(uniqueEigVecInds, inds[i])
          PcolFilled[inds[i]] <- TRUE
        }
      }
      if(!all(PcolFilled[inds])){
        gvecInds <- c(inds[1],inds[!PcolFilled[inds]])
        i <- 2; stillGoing <- TRUE
        while(i <= length(gvecInds) & stillGoing){
          newvec <- ginv(M-diag(lambrep,n+m)) %*% P[,gvecInds[i-1]]
          if(all(apply(P, 2, function(v) !equaleigvec(v,newvec)))){
            P[,gvecInds[i]] <- newvec
            PcolFilled[gvecInds[i]] <- TRUE
            for(j in 1:(i-1)) KJ[gvecInds[i-j],gvecInds[i]] <- K(lambrep,deriv=j)/factorial(j)
            i <- i + 1
          }else{
            stillGoing <- FALSE
          }
        }
      }
    }
    x <- P %*% KJ %*% solve(P, c(initSadj,rep(0,m)))
  }else{
    x <- eigM$vectors %*% (solve(eigM$vectors, c(initSadj,rep(0,m))) * K(eigM$values))
  }
  sum(transm * x[-(1:n)]) / sum(x[1:n])
}
