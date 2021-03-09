ffbs_stages = nimble::nimbleFunction(
  run = function(depths = double(1), n_timepoints = double(0),
                 transition_matrices = double(1),
                 n_bins = double(0), n_stages = double(0),
                 stage_map = double(1), alpha = double(2), beta = double(2),
                 covariates = double(2), pi_discretization = double(2),
                 n_pi = double(1), n_lambda = double(1),
                 lambda_discretization = double(2), betas_tx = double(2),
                 stage_supports = double(2), surface_bin = double(1)) {
    # use forward-filtering backwards-sampling to impute stages
    #
    # Parameters:
    #  depths - sequence of observed depth bins
    #  n_timepoints - number of observed depth bins
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
    #   surface_bin - vector of indicators for whether animal is in surface bin
    
    returnType(double(1))
    
    # initialize return
    sampled_stages <- numeric(n_timepoints, init = FALSE)
    
    # initial distribution over latent stages
    a0 <- numeric(n_stages)
    a0[1:n_stages] <- 1/n_stages
    
    # extract basic information about segment
    start_ind <- 1
    end_ind <- n_timepoints
    
    # compute likelihood for each latent stage for each depth bin tx. observed
    lmat <- dstageLik(
      x = depths, n_timepoints = n_timepoints,
      tmats = transition_matrices, n_bins = n_bins, n_stages = n_stages, 
      stage_map = stage_map, alpha = alpha, beta = beta, 
      covariates = covariates, pi_discretization = pi_discretization, 
      n_pi = n_pi, n_lambda = n_lambda, 
      lambda_discretization = lambda_discretization
    )
    
    # transition matrices for latent stages
    stx <- stageTxMats(
      betas = betas_tx, covariates = covariates, surface_bin = surface_bin,
      n_timepoints = n_timepoints
    )
    
    # sample stages
    sampled_stages[1:n_timepoints] <- ffbs(
      L = lmat * stage_supports, B = stx, a0 = a0, n_states = n_stages, 
      n_timepoints = n_timepoints
    )
      
    return(sampled_stages)
  }
)
