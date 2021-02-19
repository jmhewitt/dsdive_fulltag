time_in_state = nimble::nimbleFunction(
  run = function(x = double(1), n_timepoints = double(0)) {
    # compute a time-in-state covariate
    #
    # Parameters:
    #  x - sequence of states
    #  n_timepoints - number of observed depth bins
    # 
    # Return:
    #  vector "out" that indicates, at each index out[i], the number of 
    #  indices for which x[i] has not changed by time i.
    
    returnType(double(1))
    
    # initialize return
    time_in_state <- numeric(n_timepoints, init = FALSE)
    
    # initial value
    x_cur <- x[1]
    time_in_state[1] <- 0
    
    # enumerate runs
    for(i in 2:n_timepoints) {
      if(x[i] == x_cur) {
        # state has not changed; increment counters
        time_in_state[i] <- time_in_state[i-1] + 1
      } else {
        # state has changed; reset counters
        time_in_state[i] <- 0
        x_cur <- x[i]
      }
    }
    
    return(time_in_state)
  }
)
