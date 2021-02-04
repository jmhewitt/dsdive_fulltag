nimble_targets = list(
  
  # initial population-level regression coefficients for stage transitions
  tar_target(
    name = init_stage_tx_coefficients,
    command = {
      # stage transition matrix, used to derive multinomial logit coefs.
      o = c('deep_descent', 'deep_forage', 'deep_ascent', 'shallow_descent', 
            'shallow_ascent', 'free_surface')
      ideal_stage_txmat = matrix(0, nrow = 6, ncol = 6, 
                                 dimnames = list(rownames = o, colnames = o))
      ideal_stage_txmat['deep_descent', 'deep_forage'] = 1/15
      ideal_stage_txmat['deep_forage', 'deep_ascent'] = 1/20
      ideal_stage_txmat['deep_ascent', 'deep_descent'] = 1/15 * 1/5 * 95/100
      ideal_stage_txmat['deep_ascent', 'shallow_descent'] = 1/15 * 4/5
      ideal_stage_txmat['deep_ascent', 'free_surface'] = 1/15 * 1/5 * 5/100
      ideal_stage_txmat['shallow_descent', 'shallow_ascent'] = 1/7
      ideal_stage_txmat['shallow_ascent', 'deep_descent'] = 1/7 * 1/5 * 95/100
      ideal_stage_txmat['shallow_ascent', 'shallow_descent'] = 1/7 * 4/5
      ideal_stage_txmat['shallow_ascent', 'free_surface'] = 1/7 * 1/5 * 5/100
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
        non_surface_bin = c( # discourage certain stage tx's when not on surface
          deep_descent__deep_forage = 0,
          deep_forage__deep_ascent = 0,
          deep_ascent__deep_descent = -1e2, 
          deep_ascent__shallow_descent = -1e2,
          deep_ascent__free_surface = 0,
          shallow_descent__shallow_ascent = 0,
          shallow_ascent__deep_descent = -1e2,
          shallow_ascent__shallow_descent = -1e2,
          shallow_ascent__free_surface = 0,
          free_surface__deep_descent = 0,
          free_surface__shallow_descent = 0
        )
      )
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
        non_surface_bin = rep(0, 11)
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
        non_surface_bin = rep(0, 11)
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
      init_stage_tx_coefficients = init_stage_tx_coefficients
    )
  ),
  
  tar_target(
    name = nim_fit,
    command = fit(nim_pkg)
  )
  
)
