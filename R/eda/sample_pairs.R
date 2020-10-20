sample_pairs = function(data, times, exposure_time, response_lag, 
                        window_length, nsamples, sparse = FALSE) {
  # Parameters:
  #  data - data frame in which each row contains all data for a single obs,
  #   or a function that returns such a data frame.  This is only needed if 
  #   sparse = TRUE.
  #  times - time at which depths are observed
  #  exposure_time - time at which exposure occured
  #  response_lag - number of hours before post-exposure window
  #  window_length - number of hours to sample post-exposure
  #  sparse - TRUE to output the exact samples, FALSE to output the indices
  #    in the data argument associated with each sample
  #
  # Return: 
  #  List of resampled pre/post sample pairs, and the target pair
  
  tscale = 3600
  
  # fail analysis if tag was unexposed
  if(!is.finite(exposure_time)) {
    return(NULL)
  }
  
  if(is.function(data)) {
    data_loader = data
    data = data() %>% dplyr::select(-times)
  } else {
    data_loader = NULL
  }
  
  # format data for processing
  dat = cbind(data, times = times)
  
  # data windows
  exposed_window = exposure_time + seconds(tscale * response_lag + 
                                             c(0, tscale * window_length))
  pre_exposed_window = exposure_time - seconds(c(tscale * window_length, 0))
  
  # indices associated with windows
  exposed_inds = (exposed_window[1] <= dat$times) & 
    (dat$times < exposed_window[2])
  pre_exposed_inds = (pre_exposed_window[1] <= dat$times) &
    (dat$times < pre_exposed_window[2])
    
  # extract and format observations to analyze
  if(sparse) {
   dat.exposed = c(start_ind = which(exposed_inds)[1],
                   end_ind = tail(which(exposed_inds), 1))
   dat.pre_exposed = c(start_ind = which(pre_exposed_inds)[1],
                       end_ind = tail(which(pre_exposed_inds), 1))
  } else {
    dat.exposed = dat[exposed_inds, , drop = FALSE]
    dat.pre_exposed = dat[pre_exposed_inds, , drop = FALSE]
  }
  
  # initialize output
  res = list(
    # initialize container for Monte Carlo pre/post pairs
    resampled = list(),
    # observed pre/post pair
    observed = list(
      pre = dat.pre_exposed,
      post = dat.exposed,
      sparse = sparse,
      data_loader = data_loader
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
    sum(exposed_inds) == 0,
    # unequal number of pre/post-exposure observations
    sum(exposed_inds) != sum(pre_exposed_inds)
    # # exposed window is missing observations
    # length(unique(diff(dat.exposed$times))) > 1,
    # # unexposed window is missing observations
    # length(unique(diff(dat.pre_exposed$times))) > 1
  )) {
    return(res)
  }
  
  # number of observations in each window
  n_per_window = sum(exposed_inds)

  # total duration of pre-post-lag analysis period
  res$analysis_duration_days = difftime(
    time1 = dat$times[tail(which(exposed_inds), 1)],
    time2 = dat$times[which(pre_exposed_inds)[1]],
    units = 'days'
  )

  #
  # Monte Carlo sample baseline (i.e., unexposed) sampling distribution
  #
  
  # failure conditions
  if(res$pre_exposure_days < res$analysis_duration_days) {
    return(res)
  }

  # time range in which baseline, lagged pre/post window pairs may begin
  support = c(times[1], exposure_time - res$analysis_duration_days)

  # constants for sampling
  window_duration = seconds(tscale * window_length)
    
  # draw bootstrap samples of test statistic
  res$resampled = replicate(n = nsamples, simplify = FALSE, expr = {

    # sample start time for lagged window pair
    start_by = runif(1, min = support[1],  max = support[2])
    start_time = max(times[times <= start_by])
    
    pre_end = start_time + window_duration
    post_start = pre_end + seconds(tscale * response_lag)
    
    # determine indices for pre/post samples
    pre_inds = (start_time <= dat$times) & (times < pre_end)
    post_inds = (post_start <= dat$times) & 
      (times < post_start + window_duration)
    
    # early-return if sample is invalid
    if(any(
      # there are no post observations to analyze
      sum(post_inds) == 0,
      # there are no pre observations to analyze
      sum(pre_inds) == 0,
      # unequal number of pre/post-exposure observations
      sum(pre_inds) != sum(post_inds),
      # number of observations doesn't match target observations
      sum(pre_inds) !=  n_per_window
      # # exposed window is missing observations
      # length(unique(diff(pre_sample$times))) > 1,
      # # unexposed window is missing observations
      # length(unique(diff(post_sample$times))) > 1
    )) {
      res = NULL
    } else { 
      # extract and format sample
      if(sparse) {
        pre_sample = c(start_ind = which(pre_inds)[1],
                       end_ind = tail(which(pre_inds), 1))
        post_sample = c(start_ind = which(post_inds)[1],
                        end_ind = tail(which(post_inds), 1))
      } else {
        pre_sample = dat[pre_inds, , drop = FALSE]
        post_sample = dat[post_inds, , drop = FALSE]
      }
      res = list(pre = pre_sample, post = post_sample, sparse = sparse,
                 data_loader = data_loader)
    }
    
    res
  })
  
  res
}
