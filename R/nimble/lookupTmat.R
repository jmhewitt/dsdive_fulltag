lookupTmat = nimble::nimbleFunction(
  run = function(pi_ind = integer(0), lambda_ind = integer(0), 
                 n_pi = integer(0), n_lambda = integer(0), n_bins = integer(0), 
                 tmats = double(1)) {
    # Parameters:
    #  pi_ind - the index of the pre-specified descent preference value
    #  lambda_ind - the index of the pre-specified speed value
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    #  n_bins - the number of depth bins
    #  tmats - the transition matrices, stored in flat format
    
    returnType(double(2))
    
    # elements w/in a single transition matrix
    Ksq <- n_bins^2
    
    # index w/in target matrix (stored in col. major format)
    index <- 1:Ksq
    
    # return after offsetting wrt. grouping vars
    tm <- matrix(
      tmats[
        index + Ksq * (n_pi * (lambda_ind - 1) + pi_ind - 1)
      ], 
      nrow = n_bins
    )
    
    return(tm)
  }
)
