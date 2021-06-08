##' Theoretical circadian power calculation based on F-statistics
##'
##' Theoretically calculate power of circaidan data under indepednent Normal error assumption using F-statistics.
##' @title F-statistics-based Theoretical Circadian Power Calculation
##' @param n Sample size.
##' @param A Amplitude of the sine curve: \eqn{A * sin(2\pi/period * (phi + cts)) + C}.
##' @param sigma \sigma in the independent Normal error term \eqn{N(0,\sigma)}.
##' @param phi Phase of the sine curve. Default is 0.
##' @param period Period of the sine curve. Default is 24.
##' @param cts Circadian time design vector.
##' @param alpha Type I error control. Default is 0.05.
##' @return Theoretucal power based on F-statistics.
##' @author Wei Zong, Zhiguang Huo
##' @export
##' @examples
##' n = 100
##' cts = seq(0,24,length.out = n+1)[-1]
##' A = 1
##' sigma = 1
##' phi = 0
##' C = 10
##' F_powerCircadian(n, A, sigma, phi=0, period = 24, cts=cts, alpha = 0.05)

F_powerCircadian = function(n, A, sigma, phi=0, period = 24, cts=NULL, alpha = 0.05){

  if(is.null(cts)){
    stop("the circadian time design (cts) has to be provided for this version of R pacakge")
  }

  stopifnot(length(cts) == n)

  w = 2*pi / period

  r = 2
  df1 = r
  df2 = n - r - 1

  ncp = A^2/sigma^2 * sum(sin(w*(cts+phi))^2)

  Z = qf(1 - alpha, df1, df2, ncp = 0)
  beta = pf(Z, df1, df2, ncp = ncp)

  power = 1 - beta
  return(power)
}


