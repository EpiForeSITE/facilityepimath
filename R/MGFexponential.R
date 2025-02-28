#' Evaluate the moment generating function (MGF) of the exponential distribution or a derivative of the MGF
#'
#' @param x The value at which to evaluate the MGF
#' @param rate The rate parameter value of the exponential distribution
#' @param deriv An integer, the number of derivatives of the MGF to apply
#' @return The number resulting from the function evaluation
#' @examples
#' # MGF of an exponential distribution, evaluated at -0.1:
#' MGFexponential(-0.1, rate = 0.05)
#' # Second moment of the distribution (second derivative evaluated at zero):
#' MGFexponential(0, rate = 0.05, deriv = 2)
#' @export
MGFexponential <- function(x, rate, deriv=0)
  MGFgamma(x, rate, 1, deriv)
