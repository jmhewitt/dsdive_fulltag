sample_pairs = function(depths, depth_bins, times, exposure_time, response_lag, 
                        window_length, nsamples, diveIds) {
  # Parameters:
  #  depths - time series of depth observations
  #  depth_bins - time series of depth bin sequence
  #  times - time at which depths are observed
  #  diveIds - time series that segments the depth series into a dive sequence
  #  exposure_time - time at which exposure occured
  #  response_lag - number of hours before post-exposure window
  #  window_length - number of hours to sample post-exposure
  #
  # Return: 
  #  List of resampled pre/post sample pairs, and the target pair
  
  # fail analysis if tag was unexposed
  if(!is.finite(exposure_time)) {
    return(NULL)
  }
  
  # format data for processing
  dat = data.frame(depths = depths, depth_bins = depth_bins, times = times,
                   diveIds = diveIds)
  
  # data windows
  exposed_window = exposure_time + hours(response_lag + c(0, window_length))
  pre_exposed_window = exposure_time - hours(c(window_length, 0))
    
  # extract exposed observations to analyze
  dat.exposed = dat %>% dplyr::filter(exposed_window[1] <= times, 
                                      times < exposed_window[2])
  
  # extract unexposed observations to analyze
  dat.pre_exposed = dat %>% dplyr::filter(pre_exposed_window[1] <= times, 
                                          times < pre_exposed_window[2])
  
  # initialize output
  res = list(
    # initialize container for Monte Carlo pre/post pairs
    resampled = list(),
    # observed pre/post pair
    observed = list(
      pre = dat.pre_exposed,
      post = dat.exposed
    ),
    # time windows associated with pre/post pair
    windows = list(
      pre = pre_exposed_window,
      post = exposed_window
    ),
    # length of pre-exposure record available for "baselining" response
    pre_exposure_days = difftime(time2 = times[1], time1 = exposure_time,
                                 units = 'days'),
    # ideal amount of total duration of pre-post-lag analysis period
    analysis_duration_days = (response_lag + 2 * window_length) / 24
  )
  
  # failure conditions
  if(any(
    # there are no exposed observations to analyze
    nrow(dat.exposed) == 0,
    # unequal number of pre/post-exposure observations
    nrow(dat.exposed) != nrow(dat.pre_exposed),
    # exposed window is missing observations
    length(unique(diff(dat.exposed$times))) > 1,
    # unexposed window is missing observations
    length(unique(diff(dat.pre_exposed$times))) > 1
  )) {
    return(res)
  }
  
  # number of observations in each window
  n_per_window = nrow(dat.exposed)

  # total duration of pre-post-lag analysis period
  res$analysis_duration_days = difftime(time1 = tail(dat.exposed$times, 1),
                                        time2 = dat.pre_exposed$times[1],
                                        units = 'days')

  #
  # Monte Carlo sample baseline (i.e., unexposed) sampling distribution
  #

  # time range in which baseline, lagged pre/post window pairs may begin
  support = c(times[1], exposure_time - res$analysis_duration_days)

  # constants for sampling
  window_duration = hours(window_length)
    
  # draw bootstrap samples of test statistic
  res$resampled = replicate(n = nsamples, simplify = FALSE, expr = {

    # sample start time for lagged window pair
    start_by = runif(1, min = support[1],  max = support[2])
    start_time = max(times[times <= start_by])
    
    pre_end = start_time + window_duration
    post_start = pre_end + hours(response_lag)
    
    pre_sample = dat %>% dplyr::filter(start_time <= times, times < pre_end)
    post_sample = dat %>% dplyr::filter(post_start <= times,
                                        times < post_start + window_duration)
    
    # early-return if sample is invalid
    if(any(
      # there are no post observations to analyze
      nrow(post_sample) == 0,
      # there are no pre observations to analyze
      nrow(pre_sample) == 0,
      # unequal number of pre/post-exposure observations
      nrow(pre_sample) != nrow(post_sample),
      # number of observations doesn't match target observations
      nrow(pre_sample) !=  n_per_window,
      # exposed window is missing observations
      length(unique(diff(pre_sample$times))) > 1,
      # unexposed window is missing observations
      length(unique(diff(post_sample$times))) > 1
    )) {
      res = NULL
    } else {
      res = list(pre = pre_sample, post = post_sample)
    }
    
    res
  })
  
  res
}
