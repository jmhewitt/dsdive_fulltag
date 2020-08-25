lognormal_moment_fit = function(x) {
  
  # gamma mle parameters
  mle = fitdistr(x, 'gamma')$estimate
  
  # associated mean and sd
  mu = mle['shape'] / mle['rate']
  sigma_sq = mu / mle['rate']
  
  # log-normal parameters that match MLE mean and sd
  lognormal.params(mean = mu, var = sigma_sq)
}
