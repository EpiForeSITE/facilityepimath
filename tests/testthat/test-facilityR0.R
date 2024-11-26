test_that("facilityR0() works for Model 1", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  R0exact <- bet * mgf(0,2) / mgf(0,1) / 2

  expect_equal(facilityR0(S = 0, C = 0, A = 1, transm = bet, initS = 1, mgf = mgf), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() works for Model 2", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x) (mgf(x)-1)/x

  R0exact <- bet / gam * (1 - K(-gam) / mgf(0,1))

  expect_equal(facilityR0(S = 0, C = -gam, A = 1, transm = bet, initS = 1, mgf = mgf), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() works for Model 3", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026
  dc <- 0.00845
  eps <- 0.55

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x) (mgf(x)-1)/x

  R0exact <- bet/(dc+gam)*((1-dc*(1-eps)/(dc+gam))*(1-K(-dc-gam)/mgf(0,1))+dc*(1-eps)*mgf(0,2)/mgf(0,1)/2)

  expect_equal(facilityR0(S = 0, C = rbind(c(-dc-gam,0),c(dc,0)), A = rbind(1,0), transm = bet*c(1,1-eps), initS = 1, mgf = mgf), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() works for typical Model 4", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.0026
  dc <- 0.00845
  eps <- 0.55
  ds <- 0.07
  gamd <- 0.026
  gd <- gamd-gam

  expect_false(ds == gd)

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x) (mgf(x)-1)/x

  R0exact <- bet/(ds+dc+gam)*(1+(1-eps)/(ds-gd)*(dc*gd/(ds+dc+gam)-ds))*(1-K(-ds-dc-gam)/mgf(0,1)) +
    bet*(1-eps)/(ds-gd)*(ds*gamd/(dc+gamd)^2*(1-K(-dc-gamd)/mgf(0,1)) +
                           (dc*ds/(dc+gamd)-dc*gd/(ds+dc+gam))*mgf(0,2)/mgf(0,1)/2)

  expect_equal(facilityR0(S = rbind(c(0,0),c(0,0)),
                          C = rbind(c(-ds-dc-gam,0,0),c(ds,-dc-gamd,0),c(dc,dc,0)),
                          A = rbind(c(1,0),c(0,1-eps),c(0,0)),
                          transm = bet*c(1,1-eps,1-eps), initS = c(1,0), mgf = mgf), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() works for special case Model 4", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051
  gam <- 0.03
  dc <- 0.00845
  eps <- 0.55
  ds <- 0.07
  gamd <- 0.1
  gd <- gamd-gam

  expect_true(ds == gd)

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x, deriv = 0)
    ifelse(x == 0, mgf(0, deriv+1)/(deriv+1), ifelse(deriv == 0, (mgf(x)-1)/x, (mgf(x, deriv) - deriv * K(x, deriv-1))/x))

  R0exact <- bet/(ds+dc+gam)*((1-(1-eps)*((dc-ds)/(ds+dc+gam)+(2*ds*dc)/(ds+dc+gam)^2))*(1-K(-ds-dc-gam)/K(0))+
                                      (1-eps)*(((ds*dc)/(ds+dc+gam)-ds)*(K(-ds-dc-gam,1))/K(0)+(dc+(ds*dc)/(ds+dc+gam))*K(0,1)/K(0)))

  expect_equal(facilityR0(S = rbind(c(0,0),c(0,0)),
                          C = rbind(c(-ds-dc-gam,0,0),c(ds,-dc-gamd,0),c(dc,dc,0)),
                          A = rbind(c(1,0),c(0,1-eps),c(0,0)),
                          transm = bet*c(1,1-eps,1-eps), initS = c(1,0), mgf = mgf), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() matrix version works for Model 1", {
  prob <- 1
  rate <- 0.0285
  shape <- 1
  bet <- 0.051

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  R0exact <- bet * mgf(0,2) / mgf(0,1) / 2

  expect_equal(facilityR0(S = -rate, C = -rate, A = 1, transm = bet, initS = 1), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() matrix version works for Model 2", {
  prob <- 1
  rate <- 0.0285
  shape <- 1
  bet <- 0.051
  gam <- 0.0026

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x) (mgf(x)-1)/x

  R0exact <- bet / gam * (1 - K(-gam) / mgf(0,1))

  expect_equal(facilityR0(S = -rate, C = -gam-rate, A = 1, transm = bet, initS = 1), R0exact, tolerance = sqrt(.Machine$double.eps))
})

test_that("facilityR0() matrix version works for Model 3", {
  prob <- 1
  rate <- 0.0285
  shape <- 1
  bet <- 0.051
  gam <- 0.0026
  dc <- 0.00845
  eps <- 0.55

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  K <- function(x) (mgf(x)-1)/x

  R0exact <- bet/(dc+gam)*((1-dc*(1-eps)/(dc+gam))*(1-K(-dc-gam)/mgf(0,1))+dc*(1-eps)*mgf(0,2)/mgf(0,1)/2)

  expect_equal(facilityR0(S = -rate, C = rbind(c(-dc-gam-rate,0),c(dc,-rate)), A = rbind(1,0), transm = bet*c(1,1-eps), initS = 1), R0exact, tolerance = sqrt(.Machine$double.eps))
})
