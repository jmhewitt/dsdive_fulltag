ffbs = nimble::nimbleFunction(
  run = function(L = double(2), B = double(3), a0 = double(1), 
                 n_states = integer(0), n_timepoints = integer(0)) {
    # Forward-filtering backwards-sampling, as described in Rao and Teh (2013),
    # Algorithm 4.
    # 
    # Parameters:
    #   L - matrix of likelihoods, where each column is L^t.
    #   B - array of transition matrices, where each array entry B[,,i] is a
    #     row-stochastic matrix that defines the "from i -> to j" probabilities
    #     at each timepoint.
    #   a0 - initial distribution over latent states
    #   n_states - number of latent states
    #   n_timepoints - number of observations, i.e., ncol(L).
    # 
    # Return:
    #   a conditional sample of latent states in 1:
    
    returnType(integer(1))
    
    # initialize forward-filtering distributions
    a <- matrix(nrow = n_states, ncol = n_timepoints, init = FALSE)
    a[,1] <- a0
    
    # forward-filter
    for(t in 2:n_timepoints) {
      a[,t] <- t(B[,,t-1]) %*% (a[,t-1] * L[,t-1])
      a[,t] <- a[,t] / sum(a[,t])
    }
    
    # initialize sample
    s <- integer(n_timepoints)
    
    # backwards-sample
    b <- a[,n_timepoints] * L[,n_timepoints]
    s[n_timepoints] <- rcat(n = 1, prob = b)
    for(i in 1:(n_timepoints-1)) {  
      # nimble does not compile backwards-loops, but this is what we need
      t <- n_timepoints - i 
      # backward sample
      b <- a[,t] * L[,t] * B[,s[t+1],t]
      s[t] <- rcat(n = 1, prob = b)
    }
      
    return(s)
  }
)
