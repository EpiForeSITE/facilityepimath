#' Evaluate the moment generating function (MGF) of the gamma distribution or a derivative of the MGF
#'
#' @param x The value at which to evaluate the MGF
#' @param rate The rate parameter value of the gamma distribution
#' @param shape The shape parameter values of the gamma distribution
#' @param deriv An integer, the number of derivatives of the MGF to apply
#' @return The number resulting from the function evaluation
#' @examples
#' # MGF of a gamma distributions, evaluated at -0.1:
#' MGFgamma(-0.1, rate = 0.7, shape = 3)
#' # Second moment of the distribution (second derivative evaluated at zero):
#' MGFgamma(0, rate = 0.7, shape = 3, deriv = 2)
#' @export
MGFgamma <- function(x, rate, shape, deriv=0)
  MGFmixedgamma(x, 1, rate, shape, deriv)
