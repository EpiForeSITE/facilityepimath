test_that("facilityeq() works for Model 1", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  pa <- 0.011

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

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

test_that("facilityeq() works for Model 2", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026
  pa <- 0.011

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

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
