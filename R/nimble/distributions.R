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
  run = function(x = double(1), betas = double(2), covariates = double(2), 
                 stage_supports = double(2),
                 n_timepoints = double(0), log = logical(0, default = 0)) {
    # likelihood for sequence of latent stages, conditional on covariates and 
    # transition parameter coefficients
    #
    # Parameters:
    #  x - sequence of stages
    #  betas - matrix of coefficients for multinomial logit distributions. see 
    #    documentation for stageTxMats function for more details.
    #  covariates - matrix of covariates used to generate stage transition 
    #    matrices.  each column specifies covariates for a timepoint.
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
        # # aggregate support for stage
        # ll <- ll + log(stage_supports[x[ind], ind])
        # compute stage transition distribution for timepoint:
        #   the stage s_{ij}(t_k) is drawn using covariates at time t_k, and 
        #   the stage is Markov wrt. the previous stage
        stx <- stageTxVec(stageFrom = x[ind-1], betas = betas, 
                          covariates = covariates[, ind], log = TRUE)
        # aggregate probability of transition
        ll <- ll + stx[x[ind]]
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)

stageTxVec = nimble::nimbleFunction(
  run = function(stageFrom = double(0), betas = double(2), 
                 covariates = double(1), log = logical(0)) {
    # discrete-time stage transition vector, conditional on covariates and 
    # the stage being transitioned from
    # 
    # Parameters: 
    #   stageFrom - the row of the complete stage transition matrix to compute
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
    #   covariates - vector of covariates used to generate stage transition 
    #     probabilities.
    #   log - if TRUE, returns an approximation to the log probabilities
    # 
    # Return:
    #   vector of stage transition probabilities order of transition-to 
    #   probabilities follows:
    #     1) deep_descent
    #     2) deep_forage
    #     3) deep_ascent
    #     4) shallow_descent
    #     5) shallow_ascent
    #     6) free_surface
    
    returnType(double(1))
    
    m <- numeric(6)
    
    if(stageFrom == 1) {
      p <- multinomialLogitProbs(betas = betas[, 1, drop = FALSE], 
                                 x = covariates, log = log)
      m[2] <- p[1]  # deep_descent -> deep_forage
      m[1] <- p[2]  # deep_descent -> deep_descent
    } else if(stageFrom == 2) {
      p <- multinomialLogitProbs(betas = betas[, 2, drop = FALSE], 
                                 x = covariates, log = log)
      m[3] <- p[1]  # deep_forage -> deep_ascent
      m[2] <- p[2]  # deep_forage -> deep_forage
    } else if(stageFrom == 3) {
      p <- multinomialLogitProbs(betas = betas[, 3:5, drop = FALSE], 
                                 x = covariates, log = log)
      m[1] <- p[1]  # deep_ascent -> deep_descent
      m[4] <- p[2]  # deep_ascent -> shallow_descent
      m[6] <- p[3]  # deep_ascent -> free_surface
      m[3] <- p[4]  # deep_ascent -> deep_ascent
    } else if(stageFrom == 4) {
      p <- multinomialLogitProbs(betas = betas[, 6, drop = FALSE], 
                                 x = covariates, log = log)
      m[5] <- p[1]  # shallow_descent -> shallow_ascent
      m[4] <- p[2]  # shallow_descent -> shallow_descent
    } else if(stageFrom == 5) {
      p <- multinomialLogitProbs(betas = betas[, 7:9, drop = FALSE], 
                                 x = covariates, log = log)
      m[1] <- p[1]  # shallow_ascent -> deep_descent
      m[4] <- p[2]  # shallow_ascent -> shallow_descent
      m[6] <- p[3]  # shallow_ascent -> free_surface
      m[5] <- p[4]  # shallow_ascent -> shallow_ascent
    } else if(stageFrom == 6) {
      p <- multinomialLogitProbs(betas = betas[, 10:11, drop = FALSE], 
                                 x = covariates, log = log)
      m[1] <- p[1]  # free_surface -> deep_descent
      m[4] <- p[2]  # free_surface -> shallow_descent
      m[6] <- p[3]  # free_surface -> free_surface
    }
    
    return(m)
  }
)

stageTxMats = nimble::nimbleFunction(
  run = function(betas = double(2), covariates = double(2), 
                 n_timepoints = double(0)) {
    # discrete-time stage transition matrices, conditional on covariates
    # 
    # Parameters: 
    #   betas - matrix of coefficients for multinomial logit distributions.
    #    see stageTxVec documentation for details about the specification of the 
    #    coefficient matrix.
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
    
    m <- array(dim = c(6, 6, n_timepoints), init = FALSE)
    
    for(i in 1:n_timepoints) {
      for(j in 1:6) {
        m[j,,i] <- stageTxVec(stageFrom = j, betas = betas, 
                              covariates = covariates[,i], log = FALSE)
      }
    }
    
    return(m)
  }
)

dstageLik = nimble::nimbleFunction(
  run = function(x = double(1), n_timepoints = double(0), tmats = double(1),
                 n_bins = double(0), n_stages = double(0),
                 stage_map = double(1), alpha = double(2), beta = double(2),
                 covariates = double(2), pi_discretization = double(2),
                 n_pi = double(1), n_lambda = double(1),
                 lambda_discretization = double(2)) {
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
    
    res <- matrix(nrow = n_stages, ncol = n_timepoints, init = TRUE)
    
    # only process segments with at least one transition
    if(n_timepoints > 1) {
      # process each transition in segment; we go one fewer than the 
      # segment length b/c this accounts for the last transition
      for(ind in 1:(n_timepoints - 1)) {
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
          res[s,ind] <- lookupProb(
            movement_type = mtype, pi_ind = pi_ind, lambda_ind = lambda_ind, 
            i = x[ind], j = x[ind + 1], n_pi = n_pi, n_lambda = n_lambda, 
            n_bins = n_bins, tmats = tmats
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
  run = function(x = double(1), stages = double(1), n_timepoints = double(0),
                 stage_map = double(1), n_bins = double(0), tmats = double(1),
                 alpha = double(2), beta = double(2), covariates = double(2), 
                 pi_discretization = double(2), n_pi = double(1), 
                 n_lambda = double(1), lambda_discretization = double(2),
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
    
    # only process segments with at least one transition
    if(n_timepoints > 1) {
      # process each transition in segment; we go one fewer than the
      # segment length b/c this accounts for the last transition
      for(ind in 1:(n_timepoints - 1)) {
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
        lp <- log(lookupProb(
          movement_type = mtype, pi_ind = pi_ind, lambda_ind = lambda_ind,
          i = x[ind], j = x[ind + 1], n_pi = n_pi, n_lambda = n_lambda,
          n_bins = n_bins, tmats = tmats
        ))
        ll <- ll + lp
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)
  
dist_str_bins = paste(
  'dbins(stages, n_timepoints, stage_map, n_bins, alpha, beta, ',
        'covariates, tmats, pi_discretization, n_pi, n_lambda, ',
        'lambda_discretization)',
  sep = ''
)

dist_str_stages = 'dstages(betas, covariates, stage_supports, n_timepoints)'


nimble::registerDistributions(list(
  
  dstages = list(
    BUGSdist = dist_str_stages,
    Rdist = dist_str_stages,
    types = c('value = double(1)', 'betas = double(2)', 
              'covariates = double(2)', 'stage_supports = double(2)',
              'n_timepoints = double(0)'),
    discrete = TRUE,
    pqAvail = FALSE
  ),
  
  dbins = list(
    BUGSdist = dist_str_bins,
    Rdist = dist_str_bins,
    types = c('value = double(1)', 'stages = double(1)', 
              'n_timepoints = double(0)', 
              'stage_map = double(1)', 'n_bins = double(0)', 
              'alpha = double(2)', 'beta = double(2)', 
              'covariates = double(2)', 'tmats = double(1)', 
              'pi_discretization = double(2)', 'n_pi = double(1)',
              'lambda_discretization = double(2)', 'n_lambda = double(1)'),
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
