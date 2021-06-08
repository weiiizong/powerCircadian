##' Sample circadian time of death from the density of observed time of death in real data.
##'
##' Generate circadian time of death at various sample sizes from kernel density estimated
##' from observed time of death in real data.
##' @title Circadian Time of Death from Kernal Density of Real Data
##' @param observed_cts Observed time of death in real data. Lie between cts_lb and cts_ub.
##' @param cts_lb Lower bound of time of death.
##' @param cts_ub Lower bound of time of death.
##' @param n Number of time of death to generate.
##' @return A time of deach vector of length n
##' @author Wei Zong, Zhiguang Huo
##' @export
##' @examples
##' data(brain_data)
##' cts_kernel_generator(observed_cts=BA11_BA47.TOD, cts_lb = 0, cts_ub = 24, n = 100)

cts_kernel_generator = function(observed_cts, cts_lb = 0, cts_ub = 24, n){
  if(min(observed_cts)<cts_lb | max(observed_cts)>cts_ub){
    stop("observed_cts should lie between cts_lb and cts_ub")
  }
  period = abs(cts_ub - cts_lb)
  dens = density(observed_cts, n=1000, from = cts_lb, to = cts_ub)
  sample_cts = sample(dens$x, n, prob = dens$y)
  sample_cts2 = ifelse(sample_cts>cts_ub,sample_cts-period,sample_cts)
  sample_cts3 = ifelse(sample_cts2<cts_lb,sample_cts2+period,sample_cts2)
  return(sample_cts3)
}
