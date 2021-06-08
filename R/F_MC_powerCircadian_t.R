##' Monte-Carlo circadian power calculation based on F-statistics assuming independent t-distributed errors
##'
##' Calculate power of circaidan data from Monte-Carlo simulated data assuming independent t-distributed errors.
##' @title F-statistics-based Monte-Carlo Circadian Power Calculation assuming independent t-distributed errors
##' @param B_MC Numer of Monte-Carlo simulations. Default is 10000.
##' @param n Sample size.
##' @param A Amplitude of the sine curve: \eqn{A * sin(2\pi/period * (phi + cts)) + C}.
##' @param df Degree of freedom in the independent t-distributed error term \eqn{t(df)}.
##' @param phi Phase of the sine curve. Default is 0.
##' @param period Period of the sine curve. Default is 24.
##' @param C Offset of the sine curve. Default is 10.
##' @param cts Circadian time design vector.
##' @param alpha Type I error control. Default is 0.05.
##' @return A vector of empirical type I error and Monte-Carlo power based on F-statistics.
##' @author Wei Zong, Zhiguang Huo
##' @export
##' @examples
##' B_MC = 10000
##' n = 100
##' cts = seq(0,24,length.out = n+1)[-1]
##' A = 1
##' df = 3
##' phi = 0
##' C = 10
##' F_MC_powerCircadian_t(B_MC, n, A, df, phi=0, period = 24, cts=cts, alpha = 0.05)

F_MC_powerCircadian_t <- function(B_MC=10000, n, A, df, phi=0, period = 24, C=10, cts=NULL, alpha = 0.05){
  ##TypeI error
  MC.data0 = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rt(n, df = df)
    yy = C + error
    MC.data0 = rbind(MC.data0, yy)
  }
  pvalues.LR0 = apply(MC.data0, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_typeIerror = mean(pvalues.LR0<alpha)

  ##Power
  MC.data = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rt(n, df = df)
    yy = A*sin(2*pi/period*(cts + phi)) + C + error
    MC.data = rbind(MC.data, yy)
  }
  pvalues.LR = apply(MC.data, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_power = mean(pvalues.LR < alpha)

  return(c(F_MC_typeIerror=F_MC_typeIerror, F_MC_power=F_MC_power))
}
