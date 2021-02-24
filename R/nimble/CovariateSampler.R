sampler_Covariate = nimble::nimbleFunction(
  
  contains = nimble::sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    # for nimble samplers: nodes impacted by update
    calcNodes <- model$getDependencies(target)
    
    #
    # detect dynamic covariate dependencies induced via covariates distribution
    #
    
    # where to find the dynamic covariate controlling semi-markov behavior
    time_in_stage_row <- model$getParamExpr(target, 'time_in_stage_row')
    
    # nodes holding parameters for semi-markov behavior
    size_nodes <- deparse(model$getParamExpr(target, 'sizes'))
    prob_nodes <- deparse(model$getParamExpr(target, 'probs'))
    
    # nodes holding stages, which control semi-markov behavior
    stage_nodes <- calcNodes[
      which(model$getDistribution(calcNodes) == 'dstages')
    ]

    #
    # derived dimensions
    #
    
    # number of timepoints
    n_timepoints <- ncol(model[[target]])
    
    time_in_stage <- time_in_state(x = model[[stage_nodes]],
                                   n_timepoints = n_timepoints)

  },
  
  run = function() {
    
    # update time in stage information
    time_in_stage <- time_in_state(x = model[[stage_nodes]],
                                   n_timepoints = n_timepoints)
    
    # update dynamic covariates individually, since the hazard function depends
    # on state
    for(i in 1:n_timepoints) {
      model[[target]][time_in_stage_row,i] <<- hbinom(
        x = time_in_stage[i],
        size = model[[size_nodes]][model[[stage_nodes]][i]],
        prob = model[[prob_nodes]][model[[stage_nodes]][i]],
        # prob = .02,
        log = TRUE
      )
    }
    
    # update log probability
    calculate(model, calcNodes)
    
    # keep the model and mvsaved objects consistent
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list( reset = function() {} )
  
)
