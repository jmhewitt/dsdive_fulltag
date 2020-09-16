bootstrap_exposure = function(depths, times, dive_labels, exposure_time, 
                              nsamples, deep_depth, window_length, statistic) {
  # Parameters:
  #  depths - time series of depth observations
  #  times - time at which depths are observed
  #  dive_labels - ids to associate depth time series with dive id's
  #  exposure_time - time at which exposure occured
  #  nsamples - number of pre/post samples to generate
  #  deep_depth - threshold below which dives will be considered deep
  #  window_length - number of hours to sample post-exposure
  #  statistic - function to apply to each sampled window; the function should 
  #   work directly on a time-series vector of depths
  
  if(!is.finite(exposure_time)) {
    return(NULL)
  }
  
  # maximum range of post-exposure times to analyze
  exposure_window = exposure_time + hours(c(0,12))
  
  # index into pre-exposure observations
  pre_exposure_selector = times <= exposure_window[1]
    
  # index into post-exposure observations within exposure_window
  exposed_selector = (!pre_exposure_selector) & (times <= exposure_window[2])
  
  #
  # sequences of exposed and pre-exposed depths, and dive labels
  #
  
  first_exposed = which(exposed_selector)[1]
  n_exposed = sum(exposed_selector)
  
  depths.exposed = window(x = depths, start = first_exposed, 
                          end = first_exposed + n_exposed - 1)
  labels.exposed = window(x = dive_labels, start = first_exposed, 
                          end = first_exposed + n_exposed - 1)
  
  depths.pre_exposed = window(x = depths, end = first_exposed - 1, 
                          start = first_exposed - n_exposed)
  labels.pre_exposed = window(x = dive_labels, end = first_exposed - 1, 
                              start = first_exposed - n_exposed)
  

  # bootstrap sampling distribution
  samples = t(replicate(n = nsamples, expr = {
    
    # sample two consecutive time windows
    window_start_support = 1:(first_exposed - 2 * n_exposed)
    window_start = sample(x = window_start_support, size = 1)
    
    depths_first = window(x = depths, start = window_start, 
                          end = window_start + n_exposed - 1)
    labels_first = window(x = dive_labels, start = window_start, 
                          end = window_start + n_exposed - 1)
    
    depths_second = window(x = depths, end = window_start + 2 * n_exposed - 1, 
                           start = window_start + n_exposed)
    labels_second = window(x = dive_labels, 
                           end = window_start + 2 * n_exposed - 1, 
                           start = window_start + n_exposed)
    
    # evalute statistic
    statistic(depths_first, depths_second, labels_first, labels_second)
  }))
  
  # correct for univariate statistic functions
  if(nrow(samples) == 1) {
    samples = t(samples)
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
