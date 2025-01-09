#' Evaluate the moment generating function (MGF) of the mixed gamma distribution or a derivative of the MGF
#'
#' @param x The value at which to evaluate the MGF
#' @param prob A vector of probabilities of following each gamma distribution in the mixture
#' @param rate A vector of rate parameter values for each gamma distribution in the mixture
#' @param shape A vector of shape parameter values for each gamma distribution in the mixture
#' @param deriv An integer, the number of derivatives of the MGF to apply
#' @return The number resulting from the function evaluation
#' @examples
#' # MGF of a 40/60 mixture of two gamma distributions, evaluated at -0.1:
#' MGFmixedgamma(-0.1, prob = c(0.4,0.6), rate = c(0.4,0.7), shape = c(0.5,3))
#' # Second moment of the distribution (second derivative evaluated at zero):
#' MGFmixedgamma(0, prob = c(0.4,0.6), rate = c(0.4,0.7), shape = c(0.5,3), deriv = 2)
#' @export
MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
  sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))
