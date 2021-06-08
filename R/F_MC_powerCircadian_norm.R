##' Monte-Carlo circadian power calculation based on F-statistics assuming independent Normal errors
##'
##' Calculate power of circaidan data from Monte-Carlo simulated data assuming independent Normal errors.
##' @title F-statistics-based Monte-Carlo Circadian Power Calculation assuming independent Normal errors
##' @param B_MC Numer of Monte-Carlo simulations. Default is 10000.
##' @param n Sample size.
##' @param A Amplitude of the sine curve: \eqn{A * sin(2\pi/period * (phi + cts)) + C}.
##' @param sigma \sigma in the independent Normal error term \eqn{N(0,\sigma)}.
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
##' sigma = 1
##' phi = 0
##' C = 10
##' F_MC_powerCircadian_norm(B_MC, n, A, sigma, phi=0, period = 24, cts=cts, alpha = 0.05)

F_MC_powerCircadian_norm <- function(B_MC=10000, n, A, sigma, phi=0, period = 24, C=10, cts=NULL, alpha = 0.05){
  ##TypeI error
  MC.data0 = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rnorm(n,0,sigma)
    yy = 0 + C + error
    MC.data0 = rbind(MC.data0, yy)
  }

  pvalues.LR0 = apply(MC.data0, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_typeIerror = mean(pvalues.LR0<alpha)


  ##Power
  MC.data = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rnorm(n,0,sigma)
    yy = A*sin(2*pi/period*(cts + phi)) + C + error
    MC.data = rbind(MC.data, yy)
  }

  pvalues.LR = apply(MC.data, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_power = mean(pvalues.LR < alpha)

  return(c(F_MC_typeIerror=F_MC_typeIerror, F_MC_power=F_MC_power))
}

#internal functions from 'diffCircadian'
LR_rhythmicity <- function(tt,yy,period=24,method="LR",FN=TRUE){
  if(method=="Wald"){
    WaldTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else if(method=="LR"){
    LRTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else(("Please check your input! Method only supports 'Wald','LR','F' or 'Permutation'."))
}

LRTest <- function(tt,yy, period = 24,type=TRUE){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss

  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset

  sigma02 <- 1/(n)*sum((yy-mean(yy))^2)
  sigmaA2 <- 1/(n)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  l0 <- -n/2*log(2*pi*sigma02)-1/(2*sigma02)*sum((yy-mean(yy))^2)
  l1 <- -n/2*log(2*pi*sigmaA2)-1/(2*sigmaA2)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  dfdiff <- (n-1)-(n-3)
  LR_stat <- -2*(l0-l1)

  if(type==FALSE){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  }
  else if(type==TRUE){
    r <- 2
    k <- 3
    LR_stat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(LR_stat,df1 = r, df2 = n-k, lower.tail = F)
  }
  R2 <- 1-rss/tss
  res <- list(
    amp = amp,
    phase = phase,
    peakTime = (6 - phase) %% period,
    offset = offset,
    sigma02=sigma02, sigmaA2=sigmaA2,
    l0=l0,
    l1=l1,
    stat=LR_stat,
    pvalue=pvalue,R2=R2)
  return(res)
}

fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){

  getPred <- function(parS, tt) {
    parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p,tt)

  nls.out <-minpack.lm::nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)

  apar <- nls.out$par

  amp0 <- apar$amp
  asign <- sign(amp0)
  ## restrict amp > 0
  amp <- amp0 * asign

  phase0 <- apar$phase
  #phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period
  phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period
  offset <- apar$offset

  peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
  if(peak > period/4*3) peak = peak - period

  A <- amp0 * cos(2*pi/period * phase0)
  B <- amp0 * sin(2*pi/period * phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss/tss

  if(F){
    amp <- apar$amp
    phase <- apar$phase
    offset <- apar$offset
  }

  res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
  res
}


