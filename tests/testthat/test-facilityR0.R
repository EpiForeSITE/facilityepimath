test_that("Model 1 works", {
  prob <- c(0.58, 1-0.58)
  rate <- c(0.0285, 0.179)
  shape <- c(1, 5.74)
  bet <- 0.051

  MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
    sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob, rate, shape, deriv)

  R0exact <- bet * mgf(0,2) / mgf(0,1) / 2

  expect_equal(facilityR0(S = 0, C = 0, A = 1, transm = bet, initS = 1, mgf = mgf), R0exact)
})
