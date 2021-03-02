nimble_targets = list(
  
  # covariates that should not be estimated as these are used to 
  # encode definitions of various aspects of diving behavior
  tar_target(structural_covariates, c('non_surface_bin', 'surface_bin', 
                                      'all_shallow_depths_since_surface')),
  
  # initial population-level regression coefficients for stage transitions
  tar_target(
    name = init_stage_tx_coefficients,
    command = {
      # stage transition matrix, used to derive multinomial logit coefs.
      o = c('deep_descent', 'deep_forage', 'deep_ascent', 'shallow_descent', 
            'shallow_ascent', 'free_surface')
      ideal_stage_txmat = matrix(0, nrow = 6, ncol = 6, 
                                 dimnames = list(rownames = o, colnames = o))
      ideal_stage_txmat['deep_descent', 'deep_forage'] = 1/5
      ideal_stage_txmat['deep_forage', 'deep_ascent'] = 1/20
      ideal_stage_txmat['deep_ascent', 'deep_descent'] = 1/15 * 1/5 * 95/100
      ideal_stage_txmat['deep_ascent', 'shallow_descent'] = 1/15 * 4/5
      ideal_stage_txmat['deep_ascent', 'free_surface'] = 1/15 * 1/5 * 5/100
      ideal_stage_txmat['shallow_descent', 'shallow_ascent'] = 1/5
      ideal_stage_txmat['shallow_ascent', 'deep_descent'] = 1/5 * 1/5 * 95/100
      ideal_stage_txmat['shallow_ascent', 'shallow_descent'] = 1/5 * 4/5
      ideal_stage_txmat['shallow_ascent', 'free_surface'] = 1/5 * 1/5 * 5/100
      ideal_stage_txmat['free_surface', 'deep_descent'] = 4/5 * 9/10 
      ideal_stage_txmat['free_surface', 'shallow_descent'] = 1/5 * 9/10
      diag(ideal_stage_txmat) = 1 - rowSums(ideal_stage_txmat)
      
      # return multinomial logit coefficients for ideal stage transition matrix
      rbind(
        intercept = c(
          deep_descent__deep_forage = 
            log(ideal_stage_txmat['deep_descent', 'deep_forage']) - 
            log(ideal_stage_txmat['deep_descent', 'deep_descent']),
          deep_forage__deep_ascent =
            log(ideal_stage_txmat['deep_forage', 'deep_ascent']) - 
            log(ideal_stage_txmat['deep_forage', 'deep_forage']),
          deep_ascent__deep_descent =
            log(ideal_stage_txmat['deep_ascent', 'deep_descent']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          deep_ascent__shallow_descent =
            log(ideal_stage_txmat['deep_ascent', 'shallow_descent']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          deep_ascent__free_surface =
            log(ideal_stage_txmat['deep_ascent', 'free_surface']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          shallow_descent__shallow_ascent =
            log(ideal_stage_txmat['shallow_descent', 'shallow_ascent']) - 
            log(ideal_stage_txmat['shallow_descent', 'shallow_descent']),
          shallow_ascent__deep_descent =
            log(ideal_stage_txmat['shallow_ascent', 'deep_descent']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          shallow_ascent__shallow_descent =
            log(ideal_stage_txmat['shallow_ascent', 'shallow_descent']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          shallow_ascent__free_surface =
            log(ideal_stage_txmat['shallow_ascent', 'free_surface']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          free_surface__deep_descent =
            log(ideal_stage_txmat['free_surface', 'deep_descent']) - 
            log(ideal_stage_txmat['free_surface', 'free_surface']),
          free_surface__shallow_descent =
            log(ideal_stage_txmat['free_surface', 'shallow_descent']) - 
            log(ideal_stage_txmat['free_surface', 'free_surface'])
        ),
        depth = rep(0, 11),
        deep_depth =  rep(0, 11),
        shallow_depth = rep(0, 11),
        non_surface_bin = c( 
          deep_descent__deep_forage = 0,
          deep_forage__deep_ascent = 0,
          # don't restart a dive at depth
          deep_ascent__deep_descent = -1e2, 
          # don't switch dive types at depth
          deep_ascent__shallow_descent = -1e2,
          # don't start free_surface periods at depth
          deep_ascent__free_surface = -1e2,
          shallow_descent__shallow_ascent = 0,
          # don't switch dive types at depth
          shallow_ascent__deep_descent = -1e2,
          # don't restart a dive at depth
          shallow_ascent__shallow_descent = -1e2,
          # don't start free_surface periods at depth
          shallow_ascent__free_surface = -1e2,
          # free_surface periods to end at depths
          free_surface__deep_descent = 1e2,
          free_surface__shallow_descent = 1e2
        ),
        surface_bin = c( 
          deep_descent__deep_forage = 0,
          deep_forage__deep_ascent = 0,
          # end ascent phase at surface
          deep_ascent__deep_descent = 1e2, 
          deep_ascent__shallow_descent = 1e2,
          deep_ascent__free_surface = 1e2,
          shallow_descent__shallow_ascent = 0,
          # end ascent phase at surface
          shallow_ascent__deep_descent = 1e2,
          shallow_ascent__shallow_descent = 1e2,
          shallow_ascent__free_surface = 1e2,
          free_surface__deep_descent = 0,
          free_surface__shallow_descent = 0
        ),
        time_since_surface = rep(0, 11),
        all_shallow_depths_since_surface = c( 
          deep_descent__deep_forage = 0,
          # don't start ending a deep dive before a deep depth is reached
          deep_forage__deep_ascent = -1e2,
          deep_ascent__deep_descent = 0, 
          deep_ascent__shallow_descent = 0,
          deep_ascent__free_surface = 0,
          shallow_descent__shallow_ascent = 0,
          shallow_ascent__deep_descent = 0,
          shallow_ascent__shallow_descent = 0,
          shallow_ascent__free_surface = 0,
          free_surface__deep_descent = 0,
          free_surface__shallow_descent = 0
        ),
        daytime = rep(0, 11),
        moonlit = rep(0, 11)
      )
    }
  ),
  
  # associate the dive stage from which the stage tx. coefficients leave from
  tar_target(
    name = stages_tx_from,
    command = {
      effect_names = colnames(init_stage_tx_coefficients)
      res = sapply(
        X = strsplit(x = effect_names, split = '__'), 
        FUN = function(stages_from_to) {
          # get stage number associated stage string
          stages[ stages_from_to[[1]] ]
        }
      )
      names(res) = effect_names
      res
    }
  ),
  
  # initial population-level regression coefficients for depth bin movement
  tar_target(
    name = init_movement_coefficients,
    command = list(
      alpha = rbind(
        intercept = qlogis(c(
          deep_descent = .99,
          deep_forage = .5,
          deep_ascent = .01,
          shallow_descent = .97,
          shallow_ascent = .03,
          free_surface = .5
        )),
        depth = rep(0, 6),
        deep_depth = rep(0, 6),
        shallow_depth = rep(0, 6),
        non_surface_bin = rep(0, 6),
        surface_bin = rep(0, 6),
        time_since_surface = rep(0, 6),
        all_shallow_depths_since_surface = rep(0, 6),
        daytime = rep(0, 6),
        moonlit = rep(0, 6)
      ),
      beta = rbind(
        intercept = log(c(
          deep_descent = 1.6,
          deep_forage = .03,
          deep_ascent = 1.4,
          shallow_descent = .6,
          shallow_ascent = .6,
          free_surface = .2
        )),
        depth = rep(0, 6),
        deep_depth = rep(0, 6),
        shallow_depth = rep(0, 6),
        non_surface_bin = rep(0, 6),
        surface_bin = rep(0, 6),
        time_since_surface = rep(0, 6),
        all_shallow_depths_since_surface = rep(0, 6),
        daytime = rep(0, 6),
        moonlit = rep(0, 6)
      )
    )
  ),
  
  tar_target(
    name = nim_pkg, 
    command = flatten_tags(
      tag_list = imputed_dive, 
      transition_matrices = transition_matrices, 
      template_bins = template_bins, 
      movement_types = movement_types, 
      pi_discretization = parameter_discretization$pi, 
      lambda_discretization = parameter_discretization$lambda,
      init_movement_coefficients = init_movement_coefficients,
      init_stage_tx_coefficients = init_stage_tx_coefficients,
      stages_tx_from = stages_tx_from,
      stages = stages
    )
  ),
    )
  ),
  
  tar_target(
    name = nim_fit,
    command = fit(
      nim_pkg = nim_pkg, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 30 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir,
      structural_covariates = structural_covariates
    )
  )
  
)
