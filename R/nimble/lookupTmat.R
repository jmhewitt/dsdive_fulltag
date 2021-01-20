lookupTmat = nimble::nimbleFunction(
  run = function(movement_type = integer(0), pi_ind = integer(0), 
                 lambda_ind = integer(0), 
                 n_pi = integer(1), n_lambda = integer(1), n_bins = integer(0), 
                 tmats = double(1)) {
    # Parameters:
    #  movement_type - the type of movement (i.e., ascent, forage, descent) 
    #    the probability lookup is for
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
    
    # offset wrt. movement type
    if(movement_type > 1) {
      index <- index + Ksq * n_pi[1] * n_lambda[1]
    }
    if(movement_type > 2) {
      index <- index + Ksq * n_pi[2] * n_lambda[2]
    }
    
    # return after offsetting wrt. grouping vars
    tm <- matrix(
      tmats[
        index + Ksq * (n_pi[movement_type] * (lambda_ind - 1) + pi_ind - 1)
      ], 
      nrow = n_bins
    )
    
    return(tm)
  }
)
