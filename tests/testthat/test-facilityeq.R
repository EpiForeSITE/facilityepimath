test_that("facilityeq() works for Model 1", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  pa <- 0.011

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  Ceq <- function(alpha) 1-(1-pa)*K(-alpha)/K(0)

  alpha <- optimize(f = function(x) (x-bet*Ceq(x))^2, interval = c(0,1), tol=1e-10)$minimum

  Ceqexact <- Ceq(alpha)
  eqexact <- c(1-Ceqexact, Ceqexact)

  facilityeqtest <- facilityeq(Sfun=function(x) -x,
                               C=0,
                               Afun=function(x) x,
                               R=0,
                               transm=bet,
                               init=c(1-pa,pa),
                               mgf=mgf)

  expect_equal(facilityeqtest[1], eqexact[1], tolerance = sqrt(.Machine$double.eps))
  expect_equal(facilityeqtest[2], eqexact[2], tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityeq() matrix version works for Model 1", {
  rate <- 0.0285
  bet <- 0.051
  pa <- 0.011

  mgf <- function(x, deriv=0) MGFexponential(x, rate, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  Ceq <- function(alpha) 1-(1-pa)*K(-alpha)/K(0)

  alpha <- optimize(f = function(x) (x-bet*Ceq(x))^2, interval = c(0,1), tol=1e-10)$minimum

  Ceqexact <- Ceq(alpha)
  eqexact <- c(1-Ceqexact, Ceqexact)

  facilityeqtest <- facilityeq(Sfun=function(x) -x-rate,
                               C=-rate,
                               Afun=function(x) x,
                               R=0,
                               transm=bet,
                               init=c(1-pa,pa))

  expect_equal(facilityeqtest[1], eqexact[1], tolerance = sqrt(.Machine$double.eps))
  expect_equal(facilityeqtest[2], eqexact[2], tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityeq() matrix/erlang version works for Model 1", {
  rate <- 0.0285
  shape <- 2
  bet <- 0.051
  pa <- 0.011

  mgf <- function(x, deriv=0) MGFgamma(x, rate, shape, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  Ceq <- function(alpha) 1-(1-pa)*K(-alpha)/K(0)

  alpha <- optimize(f = function(x) (x-bet*Ceq(x))^2, interval = c(0,1), tol=1e-10)$minimum

  Ceqexact <- Ceq(alpha)
  eqexact <- c(1-Ceqexact, Ceqexact)

  facilityeqtest <- facilityeq(Sfun = function(x) rbind(c(-x-rate,0),c(rate,-x-rate)),
                               C = rbind(c(-rate,0),c(rate,-rate)),
                               Afun = function(x) rbind(c(x,0),c(0,x)),
                               R = matrix(0,2,2),
                               transm = c(bet,bet),
                               init=c(1-pa,0,pa,0))

  expect_equal(sum(facilityeqtest[1:2]), eqexact[1], tolerance = sqrt(.Machine$double.eps))
  expect_equal(sum(facilityeqtest[3:4]), eqexact[2], tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityeq() works for Model 2", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026
  pa <- 0.011

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  Ceq <- function(alpha) alpha/(alpha+gam)-(alpha/(alpha+gam)-pa)*K(-alpha-gam)/K(0)

  alpha <- optimize(f = function(x) (x-bet*Ceq(x))^2, interval = c(0,1), tol=1e-10)$minimum

  Ceqexact <- Ceq(alpha)
  eqexact <- c(1-Ceqexact, Ceqexact)

  facilityeqtest <- facilityeq(Sfun = function(x) -x,
                               C = -gam,
                               Afun = function(x) x,
                               R = gam,
                               transm = bet,
                               init = c(1-pa, pa),
                               mgf = mgf)

  expect_equal(facilityeqtest[1], eqexact[1], tolerance = sqrt(.Machine$double.eps))
  expect_equal(facilityeqtest[2], eqexact[2], tolerance = sqrt(.Machine$double.eps))
})


test_that("facilityeq() works for Model 3", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026
  dc <- 0.00845
  eps <- 0.55
  pa <- 0.011

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  lamb <- function(alpha)
    (-(alpha + gam + dc) + c(-1,1) * sqrt((alpha + gam + dc)^2 - 4 * alpha * dc))/2

  eq <- function(alpha){
    l <- lamb(alpha)
    l1 <- l[1]; l2 <- l[2]
    Seq <- ((gam*pa-(alpha+l2)*(1-pa)))/(l1-l2)*K(l1)/K(0)+
      (((alpha+l1)*(1-pa)-gam*pa))/(l1-l2)*K(l2)/K(0)
    Ceq <- (alpha+l1)/(l1-l2)*(pa-(alpha+l2)*(1-pa)/gam)* K(l1)/K(0)+
      (alpha+l2)/(l1-l2)*((alpha+l1)*(1-pa)/gam-pa)*K(l2)/K(0)
    c(Seq,Ceq,1-Seq-Ceq)
  }

  alpha <- optimize(f = function(x) (x-sum(bet*c(1,1-eps)*eq(x)[2:3]))^2, interval = c(0,1), tol=1e-10)$minimum

  eqexact <- eq(alpha)

  facilityeqtest <- facilityeq(Sfun = function(x) -x,
                               C = rbind(c(-gam-dc,0),c(dc,0)),
                               Afun = function(x) rbind(x,0),
                               R = cbind(gam,0),
                               transm = bet*c(1,1-eps),
                               init = c(1-pa,pa,0),
                               mgf = mgf)

  expect_equal(facilityeqtest[1], eqexact[1], tolerance = sqrt(.Machine$double.eps))
  expect_equal(facilityeqtest[2], eqexact[2], tolerance = sqrt(.Machine$double.eps))
  expect_equal(facilityeqtest[3], eqexact[3], tolerance = sqrt(.Machine$double.eps))
})
