sampler_Stage = nimble::nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    # for nimble samplers: nodes impacted by update
    calcNodes <- model$getDependencies(target)
    
    #
    # define "external" dependencies
    #
    
    # node holding stage transition probability coefficients 
    betas_tx_node <- extractControlElement(
      control, 'betas_tx_node', error = 'Must specify betas_tx_node'
    )
    
    # support for latent stage at each timepoint
    stage_supports <- extractControlElement(
      control, 'stage_supports', error = 'Must specify stage_supports'
    )
    
    #
    # detect FFBS dependencies induced via depths distribution
    #
    
    # depth nodes, which constrain filtering distribution for latent stages
    depth_nodes <- calcNodes[which(model$getDistribution(calcNodes) == 'dbins')]
    
    # number of depth bins model is defined over
    n_bins <- model$getParamExpr(depth_nodes, 'n_bins')
    
    # node containing depth bin transition matrices
    transition_matrices <- deparse(model$getParamExpr(depth_nodes, 'tmats'))
    alpha_nodes_as_scalar <- deparse(model$getParamExpr(depth_nodes, 'alpha'))
    
    # nodes containing movement parameters and covariates
    alpha_node <- deparse(model$getParamExpr(depth_nodes, 'alpha'))
    beta_node <- deparse(model$getParamExpr(depth_nodes, 'beta'))
    covariate_node <- deparse(model$getParamExpr(depth_nodes, 'covariates'))

    # nodes containing parameter discretization and definition
    stage_map_nodes <- deparse(model$getParamExpr(depth_nodes, 'stage_map'))
    n_pi_node <- deparse(model$getParamExpr(depth_nodes, 'n_pi'))
    n_lambda_node <- deparse(model$getParamExpr(depth_nodes, 'n_lambda'))
    pi_discretization_node <- deparse(
      model$getParamExpr(depth_nodes, 'pi_discretization')
    )
    lambda_discretization_node <- deparse(
      model$getParamExpr(depth_nodes, 'lambda_discretization')
    )
    
    #
    # derived dimensions
    #
    
    # number of latent stage types in model
    n_stages <- length(model[[stage_map_nodes]])
    
    # number of latent stages to sample
    n_timepoints <- length(model[[depth_nodes]])
    
    # number of covariates in model
    n_covariates <- nrow(model[[covariate_node]])
  
  },
  
  run = function() {
  
    # sample stages
    model[[target]] <<- ffbs_stages(
      depths = model[[depth_nodes]],
      n_timepoints = n_timepoints,
      transition_matrices = model[[transition_matrices]],
      n_bins = n_bins, 
      n_stages = n_stages,
      stage_map = model[[stage_map_nodes]],
      alpha = matrix(values(model, alpha_node), nrow = n_covariates), 
      beta = matrix(values(model ,beta_node), nrow = n_covariates),
      covariates = model[[covariate_node]],
      pi_discretization = model[[pi_discretization_node]],
      n_pi = model[[n_pi_node]],
      n_lambda = model[[n_lambda_node]],
      lambda_discretization = model[[lambda_discretization_node]],
      betas_tx = matrix(values(model, betas_tx_node), nrow = n_covariates),
      stage_supports = stage_supports
    )
    
    # update log probability
    calculate(model, calcNodes)
    
    # keep the model and mvsaved objects consistent
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list( reset = function() {} )
  
)
