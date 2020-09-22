bootstrap_exposure = function(depths, times, dive_labels, exposure_time, 
                              response_lag, nsamples, deep_depth, window_length, 
                              statistic) {
  # Parameters:
  #  depths - time series of depth observations
  #  times - time at which depths are observed
  #  dive_labels - ids to associate depth time series with dive id's
  #  exposure_time - time at which exposure occured
  #  response_lag - number of hours before post-exposure window
  #  nsamples - number of pre/post samples to generate
  #  deep_depth - threshold below which dives will be considered deep
  #  window_length - number of hours to sample post-exposure
  #  statistic - function to apply to each sampled window; the function should 
  #   work directly on a time-series vector of depths
  
  # return early if tag was unexposed
  if(!is.finite(exposure_time)) {
    return(NULL)
  }
  
  
  #
  # define pre/post exposure windows
  #
  
  # pre-exposure observations
  pre_exposure_selector = times <= exposure_time
  
  # ideal post-exposure times to analyze
  exposure_window = exposure_time + hours(response_lag + c(0, window_length))
    
  # actual post-exposure observations to analyze
  exposed_selector = (!pre_exposure_selector) & (exposure_window[1] <= times) & 
    (times <= exposure_window[2]) 
  
  # return early there are no exposed observations to analyze 
  # (i.e., if response_lag is too large, or window_length is too short)
  if(sum(exposed_selector) == 0) {
    return(NULL)
  }
  
  #
  # extract pre/post exposure indices and dive labels
  #
  
  # index at which exposure occurs
  exposure_ind = which(!pre_exposure_selector)[1]
  # index of first exposed observation to analyze
  first_exposed = which(exposed_selector)[1]
  # number of exposed observations to analyze
  n_exposed = sum(exposed_selector)
  # number of observations after exposure, but before exposure window to analyze
  n_lag = sum((!pre_exposure_selector) & (times < exposure_window[1]))
  
  depths.exposed = window(x = depths, start = first_exposed, 
                          end = first_exposed + n_exposed - 1)
  labels.exposed = window(x = dive_labels, start = first_exposed, 
                          end = first_exposed + n_exposed - 1)
  
  depths.pre_exposed = window(x = depths, end = exposure_ind - 1, 
                          start = exposure_ind - n_exposed)
  labels.pre_exposed = window(x = dive_labels, end = exposure_ind - 1, 
                              start = exposure_ind - n_exposed)
  

  #
  # bootstrap sampling distribution
  #
  
  samples = replicate(n = nsamples, expr = {
    
    # sample two lagged time windows
    window_start_support = 1:(exposure_ind - 2 * n_exposed - n_lag)
    window_start = sample(x = window_start_support, size = 1)
    
    depths_first = window(x = depths, start = window_start, 
                          end = window_start + n_exposed - 1)
    labels_first = window(x = dive_labels, start = window_start, 
                          end = window_start + n_exposed - 1)
    
    depths_second = window(x = depths, 
                           end = window_start + 2 * n_exposed - 1 + n_lag, 
                           start = window_start + n_exposed + n_lag)
    labels_second = window(x = dive_labels, 
                           end = window_start + 2 * n_exposed - 1 + n_lag, 
                           start = window_start + n_exposed + n_lag)
    
    # evalute statistic
    statistic(depths_first, depths_second, labels_first, labels_second)
  })
  
  # repackage output
  if(is.list(samples)) {
    samples = do.call(rbind, samples)
  } else if(is.numeric(samples) & !is.matrix(samples)) {
    samples = matrix(samples, ncol = 1)
  }
  
  # test statistic
  test_stat = statistic(depths.pre_exposed, depths.exposed, 
                        labels.pre_exposed, labels.exposed)
  list(null.samples = samples,
       test = test_stat, 
       p = colMeans(
         sapply(1:length(test_stat), function(i) samples[, i] >= test_stat[i])
       )
      )
}
