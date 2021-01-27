#
# support for CTMC model across depth bins
#

dtagsegment = nimble::nimbleFunction(
  run = function(x = integer(1), stages = integer(1), n_segments = integer(0),
                 segment_info = integer(2), tmats = double(1), 
                 n_bins = integer(0),
                 stage_map = integer(1), alpha = double(0), beta = double(0),
                 pi_discretization = double(2), 
                 n_pi = integer(1), n_lambda = integer(1),
                 lambda_discretization = double(2),
                 log = logical(0, default = 0)) {
    # Parameters:
    #  x - sequence of observed depth bins
    #  stages - sequence of dive stages
    #  n_segments - number of sattag segments in the data
    #  segment_info - matrix, each row of which describes a segment's support
    #  tmats - family of pre-computed transition matrices, in flat format
    #  n_bins - the number of depth bins
    #  stage_map - maps the stage value to an associated movement type in tmats
    #  alpha - transformed coefficient for descent probability
    #  beta - transformed coefficient for vertical speed
    #  pi_discretization - matrix describing descent probability grids
    #  lambda_discretization - matrix describing grid for descent probabilities
    #  log - if TRUE, then the log-likelihood is returned
    
    returnType(double(0))
    
    # initialize log-likelihood
    ll <- 0
    
    # aggregate likelihood over tag segments
    for(seg_num in 1:n_segments) {
      # only process segments with at least one transition
      if(segment_info[seg_num,2] > 1) {
        # process each transition in segment; we go one fewer than the 
        # segment length b/c this accounts for the last transition
        seg_start = segment_info[seg_num,1]
        last_tx = segment_info[seg_num,1] + segment_info[seg_num,2] - 2
        for(ind in seg_start:last_tx) {
          # extract movement type associated with dive stage
          mtype <- stage_map[stages[ind]]
          # discretize pi and lambdas
          pi_ind  <- closestIndex(
            value = ilogit(alpha),
            minval = pi_discretization[mtype, 1], 
            stepsize = pi_discretization[mtype, 2], 
            nvals = pi_discretization[mtype, 3]
          )
          lambda_ind  <- closestIndex(
            value = exp(beta),
            minval = lambda_discretization[mtype, 1], 
            stepsize = lambda_discretization[mtype, 2], 
            nvals = lambda_discretization[mtype, 3]
          )
          # aggregate likelihood for transition
          ll <- ll + log(lookupProb(
            movement_type = mtype, pi_ind = pi_ind, lambda_ind = lambda_ind, 
            i = x[ind], j = x[ind + 1], n_pi = n_pi, n_lambda = n_lambda, 
            n_bins = n_bins, tmats = tmats
          ))
          # end transition processing
        }
        # end conditional segment processing
      } 
      # end segment processing
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)