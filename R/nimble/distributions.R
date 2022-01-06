#
# support for CTMC model across depth bins
#

dcovariates = nimble::nimbleFunction(
  run = function(x = double(2), time_in_stage_row = double(0), 
                 sizes = double(1), probs = double(1), 
                 nrow = double(0), ncol = double(0),
                 log = logical(0, default = 0)) {
    # dummy, flat likelihood for covariates.  mostly used to specify 
    # model relationships used to dynamically update covariates through sampling
    #
    # Parameters:
    #  x - sequence of covariate vectors, one timepoint per column
    #  stages - sequence of stages
    #  time_in_stage_row - row in which time-in-stage covariate is stored
    #  size - stage-specific sizes of negative binomial distributions that 
    #    control the semi-markov stage transitions
    #  prob - stage-specific success probs of negative binomial distributions
    #    that control the semi-markov stage transitions
    #  log - if TRUE, then the log-likelihood is returned
    returnType(double(0))
    if(log) { return(0) } else { return(1) }
  }
)

rcovariates = nimble::nimbleFunction(
  run = function(n = integer(0), time_in_stage_row = double(0), 
                 sizes = double(1), probs = double(1), 
                 nrow = double(0), ncol = double(0)) {
    returnType(double(2))
    m <- matrix(nrow = nrow, ncol = ncol)
    return(m)
  }
)

dstages = nimble::nimbleFunction(
  run = function(x = double(1), beta_tx = double(3), covariates = double(2),
                 n_timepoints = double(0), log = logical(0, default = 0)) {
    # likelihood for sequence of latent stages, conditional on 
    # transition parameter coefficients
    #
    # Parameters:
    #  x - sequence of stages
    #  beta_tx - coefs. to define stage transition matrix via multinomial logits
    #  covariates - covariates used to define stage transition matrix
    #  n_timepoints - number of timepoints, each of which will get a stage 
    #    transition matrix.
    #  log - if TRUE, then the log-likelihood is returned
    #
    # Return:
    #  the (log-)likelihood value for the sequence of stage transitions
    
    returnType(double(0))
    
    # initialize log-likelihood
    ll <- 0
    
    # only process segments with at least one transition
    if(n_timepoints > 1) {
      # process each transition in segment; we go one fewer than the
      # segment length b/c this accounts for the last transition
      for(ind in 2:n_timepoints) {
        # compute stage transition distribution for timepoint:
        #   the stage s_{ij}(t_k) is drawn using covariates at time t_k, and 
        #   the stage is Markov wrt. the previous stage
        stx <- stageTxVec(stageFrom = x[ind-1], betas = beta_tx, 
                          covariates = covariates[, ind], 
                          log = TRUE)
        # aggregate probability of transition
        ll <- ll + stx[x[ind]]
        
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)

stageTxVec = nimble::nimbleFunction(
  run = function(stageFrom = double(0), betas = double(3), 
                 covariates = double(1), log = logical(0)) {
    # discrete-time stage transition vector, conditional on covariates and 
    # the stage being transitioned from
    # 
    # Parameters: 
    #   stageFrom - the row of the complete stage transition matrix to compute
    #   betas - array of coefficients for multinomial logit distributions.
    #    each column holds the coefficient for the logit of one pair of 
    #    allowed state transitions.  the columns should define coefs. for
    #    logit(stage_from->stage_to) = log(pi_stage_to/pi_stage_from), ordered 
    #    such that the last stage is the reference stage, and the other columns 
    #    define the other stages
    #   covariates - vector of covariates used to generate stage transition 
    #     probabilities.
    #   log - if TRUE, returns an approximation to the log probabilities
    # 
    # Return:
    #   vector of stage transition probabilities
    
    returnType(double(1))
    
    n_stages <- length(betas[1,,1])
    
    m <- multinomialLogitProbs(
      betas = betas[, stageFrom, 1:(n_stages-1)],
      x = covariates, log = log
    )
    
    return(m)
  }
)

stageTxMats = nimble::nimbleFunction(
  run = function(betas = double(3), covariates = double(2),
                 n_timepoints = double(0)) {
    # discrete-time stage transition matrices, conditional on covariates
    # 
    # Parameters: 
    #   betas - array of coefficients for multinomial logit distributions.
    #    see stageTxVec documentation for details about the specification of the 
    #    coefficient matrix.
    #   covariates - matrix of covariates used to generate stage transition 
    #     matrices.  each column specifies covariates for a timepoint.
    #   n_timepoints - number of timepoints, each of which will get a stage 
    #     transition matrix.
    # 
    # Return:
    #   array of stage transition matrices, one for each timepoint. 
    
    returnType(double(3))
    
    n_stages <- length(betas[1,,1])
    
    m <- array(dim = c(n_stages, n_stages, n_timepoints), init = FALSE)
    
    for(i in 1:n_timepoints) {
      for(j in 1:n_stages) {
        m[j,,i] <- stageTxVec(stageFrom = j, betas = betas, 
                              covariates = covariates[,i], 
                              log = FALSE)
      }
    }
    
    return(m)
  }
)


dstageLik = nimble::nimbleFunction(
  run = function(x = double(1), n_timepoints = double(0), n_stages = double(0),
                 stage_defs = double(2), lambda_inds = double(1), 
                 n_bins = double(0), n_pi = double(0), n_lambda = double(0),
                 tmats = double(1)) {
    # distribution for depth bin transitions as a function of latent stages
    # given movement parameters. used as the unnormalized likelihood for 
    # latent stages, conditional on observed movement (i.e., transitions).  the 
    # likelihood for the final stage in the target segment has a uniform 
    # distribution since no "final" transition is observed, which would provide
    # information for learning the final stage.
    #
    # the output quantifies the evidence for the value of stage s_{ij}(t_k) 
    # given an observation of the depth bin it influences \ell_{ij}(t_k+1), 
    # and is useful for forward-filtering backward-sampling of the stage
    # process.  so, here, the covariates are associated with the start of a
    # movement interval.
    #
    # Parameters:
    #  x - sequence of observed depth bins
    #  n_timepoints - number of observations
    #  n_stages - the number of possible stages
    #  stage_defs - map between stage number and pi/lambda values used
    #  lambda_inds - indices of discretized lambda values in model
    #  n_bins - the number of depth bins
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    #  tmats - family of pre-computed transition matrices, in flat format
    #
    # Return:
    #  a matrix, where each column outlines the unnormalized distribution,
    #  across stages, of an observed depth bin transition.  so each entry
    #  describes the probability of observing a specific transition under a
    #  different value for the latent stage at the time of the transition.

    returnType(double(2))
    
    res <- matrix(nrow = n_stages, ncol = n_timepoints, init = TRUE)
    
    # only process segments with at least one transition
    if(n_timepoints > 1) {
      # process each transition in segment; we go one fewer than the 
      # segment length b/c this accounts for the last transition
      for(ind in 1:(n_timepoints - 1)) {
        # consider movement from all stages
        for(s in 1:n_stages) {
          # extract pi/lambda values associated with dive stage
          pi_ind <- stage_defs[s, 1]
          lambda_ind <- lambda_inds[stage_defs[s, 2]]
          # aggregate likelihood for transition
          res[s,ind] <- lookupProb(
            pi_ind = pi_ind, lambda_ind = lambda_ind, i = x[ind], 
            j = x[ind + 1], n_pi = n_pi, n_lambda = n_lambda, n_bins = n_bins, 
            tmats = tmats
          )
          # end transition processing for stage s
        }
        # end processing observation
      }
      # add a flat likelihood for the final stage distribution
      for(i in 1:n_stages) {
        res[i,n_timepoints] <- 1
      }
      # end conditional segment processing
    } 
    
    return(res)
  }
)

dbins = nimble::nimbleFunction(
  run = function(x = double(1), n_timepoints = double(0), n_pi = double(0), 
                 n_lambda = double(0), n_bins = double(0), tmats = double(1),
                 stages = double(1), stage_defs = double(2), 
                 lambda_inds = double(1), log = logical(0, default = 0)) {
    # likelihood for depth bin transitions, conditional on latent stages and 
    # movement parameters (alpha, beta)
    # 
    # Parameters:
    #  x - sequence of observed depth bins
    #  n_timepoints - length of x
    #  stages - sequence of dive stages
    #  stage_defs - map between stage number and pi/lambda values used
    #  n_pi - number of discretized pi values used in tmats
    #  n_lambda - number of discretized lambda values used in tmats
    #  n_bins - the number of depth bins
    #  tmats - family of pre-computed transition matrices, in flat format
    #  lambda_inds - indices of discretized lambda values in model
    #  log - if TRUE, then the log-likelihood is returned
    #
    # Return:
    #  the (log-)likelihood value for the observed sequence of depth bin tx's.
    
    returnType(double(0))
    
    # initialize log-likelihood
    ll <- 0
    
    # only process segments with at least one transition
    if(n_timepoints > 1) {
      # process each transition in segment; we go one fewer than the
      # segment length b/c this accounts for the last transition
      for(ind in 1:(n_timepoints - 1)) {
        # extract pi/lambda values
        pi_ind <- stage_defs[stages[ind], 1]
        lambda_ind <- lambda_inds[stage_defs[stages[ind], 2]]
        # aggregate likelihood
        lp <- log(lookupProb(
          pi_ind = pi_ind, lambda_ind = lambda_ind, i = x[ind], j = x[ind + 1], 
          n_pi = n_pi, n_lambda = n_lambda, n_bins = n_bins, tmats = tmats
        ))
        ll <- ll + lp
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)
  
dist_str_bins = paste(
  'dbins(n_timepoints, n_pi, n_lambda, n_bins, tmats, stages, stage_defs, ',
        'lambda_inds)',
  sep = ''
)

dist_str_stages = 'dstages(beta_tx, covariates, n_timepoints)'
 
 
nimble::registerDistributions(list(
  
  dstages = list(
    BUGSdist = dist_str_stages,
    Rdist = dist_str_stages,
    types = c('value = double(1)', 'beta_tx = double(3)', 
              'covariates = double(2)', 'n_timepoints = double(0)'),
    discrete = TRUE,
    pqAvail = FALSE
  ),
  
  dbins = list(
    BUGSdist = dist_str_bins,
    Rdist = dist_str_bins,
    types = c('value = double(1)', 
              'n_timepoints = double(0)', 'n_pi = double(0)', 
              'n_lambda = double(0)', 'n_bins = double(0)', 'tmats = double(1)',
              'stages = double(1)', 'stage_defs = double(2)', 
              'lambda_inds = double(1)'),
    discrete = TRUE,
    pqAvail = FALSE
  )
  
))

dflatmat = nimble::nimbleFunction(
  run = function(x = double(2), nrow = double(0), ncol = double(0),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    if(log) { return(0) } else { return(1) }
  }
)

rflatmat = nimble::nimbleFunction(
  run = function(n = integer(0), nrow = double(0), ncol = double(0)) {
    returnType(double(2))
    m <- matrix(nrow = nrow, ncol = ncol)
    return(m)
  }
)

dflatvec = nimble::nimbleFunction(
  run = function(x = double(1), length = double(0),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    if(log) { return(0) } else { return(1) }
  }
)

rflatvec = nimble::nimbleFunction(
  run = function(n = integer(0), length = double(0)) {
    returnType(double(1))
    v <- numeric(length = length)
    return(v)
  }
)
