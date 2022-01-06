covariate_tx = function(covariates, control = list()) {
  # Parameters:
  #  covariates - an ncovariates x ntimepoints matrix containing covariates to 
  #    be transformed
  #  control - a list containing parameters that help define the transformations
  
  # control$deep_depth = 800  # what is deep?
  # control$window_len = 3600 # how many seconds should prop_recent depth do?
  # control$obs_freq = 300    # how many seconds between observations?
  # control$spline_degree = 3 # spline order
  
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
  
  # spline-transform derived covariate
  prop_poly = t(bernsteinPoly(
    x = prop_recent_deep, 
    degree = control$spline_degree, 
    intercept = TRUE,
    Boundary.knots = c(0,1)
  ))
  
  rownames(prop_poly) = paste(
    'prop_deep_poly_', 1:nrow(prop_poly), sep = ''
  )
  
  # assemble final covariate matrix
  rbind(
    covariates[c('daytime','moonlit'),],
    prop_poly
  )
}
