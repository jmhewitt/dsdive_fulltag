lookupProb = nimble::nimbleFunction(
  run = function(pi_ind = double(0), 
                 lambda_ind = double(0), i = double(0), j = double(0), 
                 n_pi = double(0), n_lambda = double(0), n_bins = double(0), 
                 tmats = double(1)) {
    # Parameters:
    #  pi_ind - the index of the pre-specified descent preference value
    #  lambda_ind - the index of the pre-specified speed value
    #  i - the depth bin being transitioned from
    #  j - the depth bin being transitioned to
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    #  n_bins - the number of depth bins
    #  tmats - the transition matrices, stored in flat format
    
    returnType(double(0))
    
    # elements w/in a single transition matrix
    Ksq <- n_bins^2
    
    # index w/in target matrix (stored in col. major format)
    index <- n_bins * (j-1) + i
    
    # return after offsetting wrt. grouping vars
    p <- tmats[
      index + Ksq * (n_pi * (lambda_ind - 1) + pi_ind - 1)
      ]
    
    return(p)
  }
)
