#
# support for CTMC model across depth bins
#

stageTxMats = nimble::nimbleFunction(
  run = function(betas = double(2), covariates = double(2), 
                 n_timepoints = integer(0)) {
    # discrete-time stage transition matrices, conditional on covariates
    # 
    # Parameters: 
    #   betas - matrix of coefficients for multinomial logit distributions.
    #    each column holds the coefficient for the logit of one pair of 
    #    allowed state transitions.  the columns should define coeffs. for
    #    logit(stage_from->stage_to) = log(pi_stage_to/pi_stage_from) in the 
    #    following order:
    #    logit(deep_descent->deep_forage), 
    #    logit(deep_forage->deep_ascent), 
    #    logit(deep_ascent->deep_descent),
    #    logit(deep_ascent->shallow_descent), 
    #    logit(deep_ascent->free_surface),
    #    logit(shallow_descent->shallow_ascent), 
    #    logit(shallow_ascent->deep_descent), 
    #    logit(shallow_ascent->shallow_descent), 
    #    logit(shallow_ascent->free_surface),
    #    logit(free_surface->deep_descent), 
    #    logit(free_surface->shallow_descent)
    #   covariates - matrix of covariates used to generate stage transition 
    #     matrices.  each column specifies covariates for a timepoint.
    #   n_timepoints - number of timepoints, each of which will get a stage 
    #     transition matrix.
    # 
    # Return:
    #   array of stage transition matrices, one for each timepoint. the row/col
    #   order both follow:
    #     1) deep_descent
    #     2) deep_forage
    #     3) deep_ascent
    #     4) shallow_descent
    #     5) shallow_ascent
    #     6) free_surface
    
    returnType(double(3))
    
    m <- array(dim = c(6, 6, n_timepoints))
    
    for(i in 1:n_timepoints) {
      p <- multinomialLogitProbs(betas = betas[, 1, drop = FALSE], 
                                 x = covariates[,i])
      m[1,2,i] <- p[1]  # deep_descent -> deep_forage
      m[1,1,i] <- p[2]  # deep_descent -> deep_descent
      p <- multinomialLogitProbs(betas = betas[, 2, drop = FALSE], 
                                 x = covariates[,i])
      m[2,3,i] <- p[1]  # deep_forage -> deep_ascent
      m[2,2,i] <- p[2]  # deep_forage -> deep_forage
      p <- multinomialLogitProbs(betas = betas[, 3:5, drop = FALSE], 
                                 x = covariates[,i])
      m[3,1,i] <- p[1]  # deep_ascent -> deep_descent
      m[3,4,i] <- p[2]  # deep_ascent -> shallow_descent
      m[3,6,i] <- p[3]  # deep_ascent -> free_surface
      m[3,3,i] <- p[4]  # deep_ascent -> deep_ascent
      p <- multinomialLogitProbs(betas = betas[, 6, drop = FALSE], 
                                 x = covariates[,i])
      m[4,5,i] <- p[1]  # shallow_descent -> shallow_ascent
      m[4,4,i] <- p[2]  # shallow_descent -> shallow_descent
      p <- multinomialLogitProbs(betas = betas[, 7:9, drop = FALSE], 
                                 x = covariates[,i])
      m[5,1,i] <- p[1]  # shallow_ascent -> deep_descent
      m[5,4,i] <- p[2]  # shallow_ascent -> shallow_descent
      m[5,6,i] <- p[3]  # shallow_ascent -> free_surface
      m[5,5,i] <- p[4]  # shallow_ascent -> shallow_ascent
      p <- multinomialLogitProbs(betas = betas[, 10:11, drop = FALSE], 
                                 x = covariates[,i])
      m[6,1,i] <- p[1]  # free_surface -> deep_descent
      m[6,4,i] <- p[2]  # free_surface -> shallow_descent
      m[6,6,i] <- p[3]  # free_surface -> free_surface
    }
    
    return(m)
  }
)

dstageLik = nimble::nimbleFunction(
  run = function(x = integer(1), segment_num = integer(0),
                 segment_info = integer(2), tmats = double(1),
                 n_bins = integer(0), n_stages = integer(0),
                 stage_map = integer(1), alpha = double(2), beta = double(2),
                 covariates = double(2), pi_discretization = double(2),
                 n_pi = integer(1), n_lambda = integer(1),
                 lambda_discretization = double(2)) {
    # distribution for depth bin transitions as a function of latent stages
    # given movement parameters. used as the unnormalized likelihood for 
    # latent stages, conditional on observed movement (i.e., transitions).  the 
    # likelihood for the final stage in the target segment has a uniform 
    # distribution since no "final" transition is observed, which would provide
    # information for learning the final stage.
    #
    # Parameters:
    #  x - sequence of observed depth bins
    #  segment_num - the sattag segment to process
    #  segment_info - matrix, each row of which describes a segment's support
    #  tmats - family of pre-computed transition matrices, in flat format
    #  n_bins - the number of depth bins
    #  n_stages - the number of possible stages
    #  stage_map - maps the stage value to an associated movement type in tmats
    #  alpha - matrix of coefficients for descent probability, each column 
    #   specifies coefficients for a different stage
    #  beta - matrix of coefficients for vertical speed, each column 
    #   specifies coefficients for a different stage
    #  covariates - matrix of covariates associated with depth bins. each column
    #   contains all covariates for each timepoint
    #  pi_discretization - matrix describing descent probability grids
    #  lambda_discretization - matrix describing grid for descent probabilities
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    
    #
    # Return:
    #  a matrix, where each column outlines the unnormalized distribution,
    #  across stages, of an observed depth bin transition.  so each entry
    #  describes the probability of observing a specific transition under a
    #  different value for the latent stage at the time of the transition.

    returnType(double(2))
    
    res <- matrix(nrow = n_stages, ncol = segment_info[segment_num, 2], 
                  init = FALSE)
    
    # only process segments with at least one transition
    if(segment_info[segment_num,2] > 1) {
      # process each transition in segment; we go one fewer than the 
      # segment length b/c this accounts for the last transition
      seg_start = segment_info[segment_num,1]
      last_tx = segment_info[segment_num,1] + segment_info[segment_num,2] - 2
      ind_range = seg_start:last_tx
      for(it in 1:length(ind_range)) {
        ind <- ind_range[it]
        # consider movement from all stages
        for(s in 1:n_stages) {
          # extract movement type associated with dive stage
          mtype <- stage_map[s]
          # compute pi and lambdas
          pi_val <- ilogit(inprod(alpha[,s], covariates[,ind]))
          lambda_val <- exp(inprod(beta[,s], covariates[,ind]))
          # discretize pi and lambdas
          pi_ind  <- closestIndex(
            value = pi_val,
            minval = pi_discretization[mtype, 1], 
            stepsize = pi_discretization[mtype, 2], 
            nvals = pi_discretization[mtype, 3]
          )
          lambda_ind  <- closestIndex(
            value = lambda_val,
            minval = lambda_discretization[mtype, 1], 
            stepsize = lambda_discretization[mtype, 2], 
            nvals = lambda_discretization[mtype, 3]
          )
          # aggregate likelihood for transition
          res[s,it] <- lookupProb(
            movement_type = mtype, pi_ind = pi_ind, lambda_ind = lambda_ind, 
            i = x[ind], j = x[ind + 1], n_pi = n_pi, n_lambda = n_lambda, 
            n_bins = n_bins, tmats = tmats
          )
          # end transition processing for stage s
        }
        # end processing observation
      }
      # add a flat likelihood for the final stage transition
      for(i in 1:n_stages) {
        res[i,segment_info[segment_num, 2]] <- 1
      }
      # end conditional segment processing
    } 
    
    return(res)
  }
)

dbins = nimble::nimbleFunction(
  run = function(x = integer(1), stages = integer(1), n_segments = integer(0),
                 segment_info = integer(2), tmats = double(1), 
                 n_bins = integer(0),
                 stage_map = integer(1), alpha = double(2), beta = double(2),
                 covariates = double(2), pi_discretization = double(2), 
                 n_pi = integer(1), n_lambda = integer(1),
                 lambda_discretization = double(2),
                 log = logical(0, default = 0)) {
    # likelihood for depth bin transitions, conditional on latent stages and 
    # movement parameters (alpha, beta)
    # 
    # Parameters:
    #  x - sequence of observed depth bins
    #  stages - sequence of dive stages
    #  n_segments - number of sattag segments in the data
    #  segment_info - matrix, each row of which describes a segment's support
    #  tmats - family of pre-computed transition matrices, in flat format
    #  n_bins - the number of depth bins
    #  stage_map - maps the stage value to an associated movement type in tmats
    #  alpha - matrix of coefficients for descent probability, each column 
    #   specifies coefficients for a different stage
    #  beta - matrix of coefficients for vertical speed, each column 
    #   specifies coefficients for a different stage
    #  covariates - matrix of covariates associated with depth bins. each column
    #   contains all covariates for each timepoint
    #  pi_discretization - matrix describing descent probability grids
    #  lambda_discretization - matrix describing grid for descent probabilities
    #  log - if TRUE, then the log-likelihood is returned
    #
    # Return:
    #  the (log-)likelihood value for the observed sequence of depth bin tx's.
    
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
          s <- stages[ind]
          mtype <- stage_map[s]
          # compute pi and lambdas
          pi_val <- ilogit(inprod(alpha[,s], covariates[,ind]))
          lambda_val <- exp(inprod(beta[,s], covariates[,ind]))
          # discretize pi and lambdas
          pi_ind  <- closestIndex(
            value = pi_val,
            minval = pi_discretization[mtype, 1], 
            stepsize = pi_discretization[mtype, 2], 
            nvals = pi_discretization[mtype, 3]
          )
          lambda_ind  <- closestIndex(
            value = lambda_val,
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
