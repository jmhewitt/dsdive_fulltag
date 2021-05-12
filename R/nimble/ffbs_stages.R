ffbs_stages = nimble::nimbleFunction(
  run = function(depths = double(1), n_timepoints = double(0), 
                 n_stages = double(0), stage_defs = double(2), 
                 lambda_inds = double(1), n_bins = double(0), 
                 n_pi = double(0), n_lambda = double(0),
                 transition_matrices = double(1),
                 betas_tx = double(3), covariates = double(2)) {
    # use forward-filtering backwards-sampling to impute stages
    #
    # Parameters:
    #  depths - sequence of observed depth bins
    #  n_timepoints - number of observed depth bins
    #  n_stages - the number of possible stages
    #  stage_defs - map between stage number and pi/lambda values used
    #  lambda_inds - indices of discretized lambda values in model
    #  n_bins - the number of depth bins
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    #  transition_matrices - family of pre-computed tx. matrices, in flat format
    
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
      x = depths, n_timepoints = n_timepoints, n_stages = n_stages, 
      stage_defs = stage_defs, lambda_inds = lambda_inds, n_bins = n_bins,
      n_pi = n_pi, n_lambda = n_lambda, tmats = transition_matrices
    )
    
    # transition matrices for latent stages
    stx <- stageTxMats(
      betas = betas_tx, covariates = covariates, n_timepoints = n_timepoints
    )
    
    # sample stages
    sampled_stages[1:n_timepoints] <- ffbs(
      L = lmat, B = stx, a0 = a0, n_states = n_stages, 
      n_timepoints = n_timepoints
    )
      
    return(sampled_stages)
  }
)

if(!exists('ffbs_stages')) {
  ffbs_stages_c = nimble::compileNimble(ffbs_stages)
}