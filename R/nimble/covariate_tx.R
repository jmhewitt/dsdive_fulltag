covariate_tx = function(covariates, control = list()) {
  # Parameters:
  #  covariates - an ncovariates x ntimepoints matrix containing covariates to 
  #    be transformed
  #  control - a list containing parameters that help define the transformations
  
  # control$deep_depth = 800  # what is deep?
  # control$window_len = 3600 # how many seconds should prop_recent depth do?
  # control$obs_freq = 300    # how many seconds between observations?
  # control$poly_degree = 3   # polynomial order
  # control$lon               # proximal longitude of data collection
  # control$lat               # proximal latitude of data collection
  
  # validate minimum required contents for control list
  required = c('deep_depth', 'window_len', 'obs_freq', 'poly_degree')
  if(!all(required %in% names(control))) {
    stop(paste(
      'control list missing required element(s):',
      paste(setdiff(required, names(control)), collapse = ', ')
    ))
  }
  
  # derived covariate: proportion of recent observations at depth
  prop_recent_deep = sapply(1:ncol(covariates), function(i) {
    # window defining recent observations
    window = covariates['time', i] - c(control$window_len, 0)
    # data indices of recent observations
    window_inds = (window[1] <= covariates['time',]) & 
                  (covariates['time',] < window[2])
    # stop processing if we don't have a complete record of recent observations
    if(sum(window_inds) != control$window_len / control$obs_freq) {
      NA
    } else {
      mean(covariates['depth',window_inds] >= control$deep_depth)
    }
  })
  
  # polynomial-transform derived covariate (no intercept)
  prop_poly = t(outer(prop_recent_deep, 1:control$poly_degree, FUN = '^'))
  
  rownames(prop_poly) = paste(
    'prop_deep_poly_', 1:nrow(prop_poly), sep = ''
  )
  
  # derived covariate: total vertical movement over last 1h
  total_recent_vertical = sapply(1:ncol(covariates), function(i) {
    # window defining recent observations
    window = covariates['time', i] - c(control$window_len, 0)
    # data indices of recent observations
    window_inds = (window[1] <= covariates['time',]) & 
      (covariates['time',] < window[2])
    # stop processing if we don't have a complete record of recent observations
    if(sum(window_inds) != control$window_len / control$obs_freq) {
      NA
    } else {
      sum(abs(diff(covariates['depth',window_inds])))
    }
  })
  
  # polynomial-transform derived covariate (no intercept)
  vertical_poly = t(
    outer(total_recent_vertical, 1:control$poly_degree, FUN = '^')
  )
  
  rownames(vertical_poly) = paste(
    'total_vertical_poly_', 1:nrow(vertical_poly), sep = ''
  )
  
  # (re-)compute celestial covariates is lon/lat information is provided
  if(all(!is.null(control$lon), !is.null(control$lat))) {
    dates = as.POSIXct(
      covariates['time',], origin = '1970-01-01 00:00.00 UTC', tz = 'UTC'
    )
    covariates['daytime',] = daytime(
      date = dates, lat = control$lat, lon = control$lon
    )
    covariates['moonlit',] = moonlit(
      date = dates, lat = control$lat, lon = control$lon
    )
  }
  
  # transform celestial covariates to scientifically meaningful factors
  celestial = rbind(
    daytime = covariates['daytime',],
    night_dark = 
      (covariates['daytime',] == FALSE) & (covariates['moonlit',] == FALSE),
    night_moonlit = 
      (covariates['daytime',] == FALSE) & (covariates['moonlit',] == TRUE)
  )
  
  # assemble final covariate matrix
  rbind(
    intercept = 1,
    celestial,
    prop_poly,
    vertical_poly
  )
}
