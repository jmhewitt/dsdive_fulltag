sampler_Stage = nimble::nimbleFunction(
  
  contains = nimble::sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    # for nimble samplers: nodes impacted by update
    calcNodes <- model$getDependencies(target)
    
    #
    # detect FFBS dependencies induced via depths distribution
    #

    # depth nodes, which constrain filtering distribution for latent stages
    depth_nodes <- calcNodes[which(model$getDistribution(calcNodes) == 'dbins')]
    
    # number of depth bins model is defined over
    n_bins <- model$getParamExpr(depth_nodes, 'n_bins')
    
    # node containing depth bin and stage transition matrices
    transition_matrices <- deparse(model$getParamExpr(depth_nodes, 'tmats'))
    txmat_stages_nodes <- deparse(model$getParamExpr(target, 'txmat'))
    
    # node containing modeled lambda values
    lambda_inds_nodes <- deparse(model$getParamExpr(depth_nodes, 'lambda_inds'))
    
    # nodes containing parameter discretization and definition
    stage_defs_nodes <- deparse(model$getParamExpr(depth_nodes, 'stage_defs'))
    n_pi <- model$getParamExpr(depth_nodes, 'n_pi')
    n_lambda <- model$getParamExpr(depth_nodes, 'n_lambda')

    
    #
    # derived dimensions
    #
    
    # number of latent stage types in model
    n_stages <- nrow(model[[stage_defs_nodes]])
    
    # number of latent stages to sample
    n_timepoints <- length(model[[depth_nodes]])
  
  },
  
  run = function() {
  
    # sample stages
    model[[target]] <<- ffbs_stages(
      depths = model[[depth_nodes]],
      n_timepoints = n_timepoints,
      n_stages = n_stages,
      stage_defs = model[[stage_defs_nodes]],
      lambda_inds = model[[lambda_inds_nodes]],
      n_bins = n_bins, 
      n_pi = n_pi,
      n_lambda = n_lambda,
      transition_matrices = model[[transition_matrices]],
      txmat_stages = model[[txmat_stages_nodes]]
    )
    
    # update log probability
    calculate(model, calcNodes)
    
    # keep the model and mvsaved objects consistent
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list( reset = function() {} )
  
)
