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
  
  # fail analysis if tag was unexposed
  if(!is.finite(exposure_time)) {
    return(NULL)
  }
  
  
  #
  # extract data for analysis
  #
  
  # format data for processing
  dat = data.frame(depths = depths, times = times, labels = dive_labels)
  
  # data windows
  exposed_window = exposure_time + hours(response_lag + c(0, window_length))
  pre_exposed_window = exposure_time - hours(c(window_length, 0))
    
  # extract exposed observations
  dat.exposed = dat %>% dplyr::filter(exposed_window[1] <= times, 
                                      times < exposed_window[2])
  
  # extract observations immediately before exposure
  dat.pre_exposed = dat %>% dplyr::filter(pre_exposed_window[1] <= times, 
                                          times < pre_exposed_window[2])
  
  
  # initialize output
  res = list(
    pre_exposed_window = pre_exposed_window,
    exposed_window = exposed_window,
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
  analysis_duration_days = difftime(time1 = tail(dat.exposed$times, 1),
                                    time2 = dat.pre_exposed$times[1],
                                    units = 'days')


  #
  # bootstrap baseline (i.e., unexposed) sampling distribution
  #

  # time range in which baseline, lagged pre/post window pairs may begin
  support = c(times[1], exposure_time - analysis_duration_days)

  # draw bootstrap samples of test statistic
  samples = replicate(n = nsamples, expr = {

    # sample start time for lagged window pair
    start_by = runif(1, min = support[1],  max = support[2])
    start_time = max(times[times <= start_by])
    
    pre_sample = dat %>% dplyr::filter(
      start_time <= times,
      times < start_time + hours(window_length)
    )
    
    post_sample = dat %>% dplyr::filter(
      start_time + hours(response_lag) + hours(window_length) <= times,
      times < start_time + hours(response_lag) + hours(2 * window_length)
    )
    
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
      return(NULL)
    }
    
    # evalute statistic
    statistic(x = pre_sample$depths, y = post_sample$depths,  
              x_labels = pre_sample$labels, y_labels = post_sample$labels)
  })
  
  # repackage output
  if(is.list(samples)) {
    samples = do.call(rbind, samples)
  } else if(is.numeric(samples) & !is.matrix(samples)) {
    samples = matrix(samples, ncol = 1)
  }
  
  # test statistic
  test_stat = statistic(x = dat.pre_exposed$depths,
                        x_labels = dat.pre_exposed$labels,
                        y = dat.exposed$depth,
                        y_labels = dat.exposed$labels)
  
  # add to output
  res$null.samples = samples
  res$test = test_stat
  res$p = colMeans(
    sapply(1:length(test_stat), function(i) samples[, i] >= test_stat[i])
  )
  res$analysis_duration_days = analysis_duration_days

  res
}
