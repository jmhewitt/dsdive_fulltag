sattag_eda = function(depths, depth_bins, dive_ids, dive_types, times, 
                      exposure_time, response_lag, window_length, nsamples, 
                      conditional_class = NULL, stat) {
  # 
  # response_lags and window_lengths will be crossed
  #
  # Parameters:
  #  depths - observed depth sequence
  #  times - sequence of observation times
  #  exposure_time - time at which exposure occured
  #  response_lag - number of hours before post-exposure window
  #  window_length - number of hours to sample post-exposure
  #  stat - function to summarize a segment of a window

  tscale = 3600
  
  # fail analysis if tag was unexposed
  if(!is.finite(exposure_time)) {
    return(NULL)
  }

  # format data for processing
  dat = data.frame(depths = depths, times = times, dive_id = dive_ids, 
                   depth_bins = depth_bins) %>% 
    dplyr::mutate(
      ddepths = c(0, diff(depths)),
      ddepths.sign = sign(ddepths),
      descending = ddepths.sign > 0,
      ascending = ddepths.sign < 0,
      no_change = ddepths.sign == 0
    ) %>% 
    dplyr::left_join(
      dive_types, by = c('dive_id'='diveId')
    )
  
  # data windows
  exposed_window = exposure_time + 
    seconds(tscale * response_lag + c(0, tscale * window_length))
  pre_exposed_window = exposure_time - seconds(c(tscale * window_length, 0))
  
  # indices associated with windows
  exposed_inds = (exposed_window[1] <= dat$times) & 
    (dat$times < exposed_window[2])
  pre_exposed_inds = (pre_exposed_window[1] <= dat$times) &
    (dat$times < pre_exposed_window[2])
  
  # determine dive class at time of exposure
  if(!is.null(conditional_class)) { 
    exposed_class = conditional_class(
      dat[which(exposed_inds)[1], , drop = FALSE]
    )
  }
  
  
  #
  # quick data validation
  #
  
  # number of observations in each window
  n_per_window = sum(exposed_inds)
  
  # length of pre-exposure record available for "baselining" response
  pre_exposure_days = difftime(time2 = times[1], time1 = exposure_time,
                               units = 'days')
  
  # ideal amount of total duration of pre-post-lag analysis period
  analysis_duration_days = (response_lag + 2 * window_length) / 24
  
  # failure conditions
  if(any(
    # there are no exposed observations to analyze
    n_per_window == 0,
    # unequal number of pre/post-exposure observations
    n_per_window != sum(pre_exposed_inds),
    # not enough pre-exposure data to build baseline distribution
    pre_exposure_days < analysis_duration_days
  )) { 
    return(NULL)
  }
  
  
  #
  # Monte Carlo sample baseline (i.e., unexposed) sampling distribution
  #
  
  # time range and associated inds where lagged pre/post window pairs may begin
  support = c(times[1], exposure_time - analysis_duration_days)
  start_inds = which((support[1] <= times) & (times < support[2]))
  
  # constants for sampling
  window_duration = seconds(tscale * window_length)
  
  # restrict indices at which pre/post window pairs may begin
  if(!is.null(conditional_class)) {
    # determine which window pairs match the conditional class
    correct_class = sapply(start_inds, function(sind) {
      # start time for lagged window pair
      start_time = times[sind]
      # indices for lagged window pair's post sample
      pre_end = start_time + window_duration
      post_start = pre_end + seconds(tscale * response_lag)
      post_inds = (post_start <= times) & (times < post_start + window_duration)
      # check to see if the post window class matches the conditional class
      rowind = which(post_inds)[1]
      ifelse(is.na(rowind), FALSE, conditional_class(
        dat[rowind, , drop = FALSE]
      )  == unlist(exposed_class))
    })
    
    # only keep pre/post windows that match the conditional class
    start_inds = start_inds[correct_class]
  }
  
  # draw bootstrap samples of test statistics
  resampled = do.call(rbind, replicate(n = nsamples, expr = {
    
    # sample start time for lagged window pair
    start_time = times[sample(x = start_inds, size = 1)]
    
    pre_end = start_time + window_duration
    post_start = pre_end + seconds(tscale * response_lag)
    
    # determine indices for pre/post samples
    pre_inds = (start_time <= times) & (times < pre_end)
    post_inds = (post_start <= times) & (times < post_start + window_duration)
    
    if(any(
      # there are no post observations to analyze
      sum(post_inds) == 0,
      # there are no pre observations to analyze
      sum(pre_inds) == 0,
      # unequal number of pre/post-exposure observations
      sum(pre_inds) != sum(post_inds),
      # number of observations doesn't match target observations
      sum(pre_inds) !=  n_per_window
    )) {
      NULL
    } else {
      # return difference in pre/post summary statistics
      stat(dat[post_inds, , drop = FALSE]) - stat(dat[pre_inds, , drop = FALSE])
    }
  }))
  
  # observed value of test statistic
  observed = stat(dat[exposed_inds, , drop = FALSE]) - 
    stat(dat[pre_exposed_inds, , drop = FALSE])
  
  # right-tail probability
  p.right = rowMeans(apply(X = resampled, MARGIN = 1, FUN = function(r) {
    r >= observed
  }))
  
  # left-tail probability
  p.left = rowMeans(apply(X = resampled, MARGIN = 1, FUN = function(r) {
    r <= observed
  }))
  
  # SSE upper-tail probability
  p.sse = rowMeans(apply(X = resampled, MARGIN = 1, FUN = function(r) {
    r^2 >= observed^2
  }))
  
  # package results
  list(
    # Monte Carlo samples from null distribution 
    resampled = resampled,
    # observed value of test statistic
    observed = observed,
    # p-values
    p.right = p.right,
    p.left = p.left,
    p.sse = p.sse,
    # copy-forward key arguments
    response_lag = response_lag,
    window_length = window_length,
    conditional = !is.null(conditional_class)
  )
   
}
