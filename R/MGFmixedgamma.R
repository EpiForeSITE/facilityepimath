#' Evaluate the moment generating function (MGF) of the mixed gamma distribution or a derivative of the MGF
#'
#' @param x The value at which to evaluate the MGF
#' @param prob A vector of probabilities of following each gamma distribution in the mixture
#' @param rate A vector of rate parameter values for each gamma distribution in the mixture
#' @param shape A vector of shape parameter values for each gamma distribution in the mixture
#' @param deriv An integer, the number of derivatives of the MGF to apply
#' @return The number resulting from the function evaluation
#' @export
MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
  sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))
