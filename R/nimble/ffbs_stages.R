ffbs_stages = nimble::nimbleFunction(
  run = function(n_segments = integer(0), segment_info = integer(2), 
                 depths = integer(1), 
                 transition_matrices = double(1),
                 n_bins = integer(0), n_stages = integer(0),
                 stage_map = integer(1), alpha = double(2), beta = double(2),
                 covariates = double(2), pi_discretization = double(2),
                 n_pi = integer(1), n_lambda = integer(1),
                 lambda_discretization = double(2), betas_tx = double(2),
                 stage_supports = double(2)) {
    # use forward-filtering backwards-sampling to impute stages
    #
    # Parameters:
    #  n_segments - number of segments for which to impute stages
    #  segments - matrix, each row of which defines the starting index, length, 
    #   and subject id for each segment in the data
    #  depths - sequence of observed depth bins
    #  segment_info - matrix, each row of which describes a segment's support
    #  transition_matrices - family of pre-computed tx. matrices, in flat format
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
    #  betas_tx - matrix of coefficients for multinomial logit distributions 
    #   that govern transitions between latent stages.  see documentation for 
    #   stageTxMats for complex description of the required format.
    #  stage_supports - matrix, each column of which specifies support for 
    #   latent stage distribution at each timepoint
    
    returnType(integer(1))
    
    # initialize return
    sampled_stages <- integer(sum(segment_info[1:n_segments, 2]), init = FALSE)
    
    # initial distribution over latent stages
    a0 <- numeric(n_stages)
    a0[1:n_stages] <- 1/n_stages
    
    # impute stages for each segment
    for(segment_num in 1:n_segments) {
    
      # extract basic information about segment
      start_ind <- segment_info[segment_num, 1]
      n_timepoints <-  segment_info[segment_num, 2]
      end_ind <- start_ind + n_timepoints - 1
      
      # compute likelihood for each latent stage for each depth bin tx. observed
      lmat <- dstageLik(
        x = depths, segment_num = segment_num, segment_info = segment_info, 
        tmats = transition_matrices, n_bins = n_bins, n_stages = n_stages, 
        stage_map = stage_map, alpha = alpha, beta = beta, 
        covariates = covariates, pi_discretization = pi_discretization, 
        n_pi = n_pi, n_lambda = n_lambda, 
        lambda_discretization = lambda_discretization
      )
      
      # transition matrices for latent stages
      stx <- stageTxMats(
        betas = betas_tx, covariates = covariates, n_timepoints = n_timepoints
      )
      
      # sample stages
      sampled_stages[start_ind:end_ind] <- ffbs(
        L = lmat * stage_supports[, start_ind:end_ind], B = stx, 
        a0 = a0, n_states = n_stages, n_timepoints = n_timepoints
      )
      
    }
    
    return(sampled_stages)
  }
)
