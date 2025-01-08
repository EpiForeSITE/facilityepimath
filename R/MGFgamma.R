#' Evaluate the moment generating function (MGF) of the gamma distribution or a derivative of the MGF
#'
#' @param x The value at which to evaluate the MGF
#' @param rate The rate parameter value of the gamma distribution
#' @param shape The shape parameter values of the gamma distribution
#' @param deriv An integer, the number of derivatives of the MGF to apply
#' @return The number resulting from the function evaluation
#' @export
MGFgamma <- function(x, rate, shape, deriv=0)
  MGFmixedgamma(x, 1, rate, shape, deriv)
