equaleigvec <- function(x, y , tol = sqrt(.Machine$double.eps)){
  xa <- x/sqrt(sum(x^2))
  ya <- y/sqrt(sum(y^2))
  max(abs(xa-ya)) < tol | max(abs(xa+ya)) < tol
}
