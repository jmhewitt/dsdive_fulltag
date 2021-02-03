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
          log(ideal_stage_txmat['deep_descent', 'deep_forage']) - 
            log(ideal_stage_txmat['deep_descent', 'deep_descent']),
          log(ideal_stage_txmat['deep_forage', 'deep_ascent']) - 
            log(ideal_stage_txmat['deep_forage', 'deep_forage']),
          log(ideal_stage_txmat['deep_ascent', 'deep_descent']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          log(ideal_stage_txmat['deep_ascent', 'shallow_descent']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          log(ideal_stage_txmat['deep_ascent', 'free_surface']) - 
            log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
          log(ideal_stage_txmat['shallow_descent', 'shallow_ascent']) - 
            log(ideal_stage_txmat['shallow_descent', 'shallow_descent']),
          log(ideal_stage_txmat['shallow_ascent', 'deep_descent']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          log(ideal_stage_txmat['shallow_ascent', 'shallow_descent']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          log(ideal_stage_txmat['shallow_ascent', 'free_surface']) - 
            log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
          log(ideal_stage_txmat['free_surface', 'deep_descent']) - 
            log(ideal_stage_txmat['free_surface', 'free_surface']),
          log(ideal_stage_txmat['free_surface', 'shallow_descent']) - 
            log(ideal_stage_txmat['free_surface', 'free_surface'])
        ),
        depth = rep(0, 11),
        deep_depth =  rep(0, 11),
        shallow_depth = rep(0, 11)
      )
    }
  ),
  
  # initial population-level regression coefficients for depth bin movement
  tar_target(
    name = init_movement_coefficients,
    command = list(
      alpha = rbind(
        intercept = qlogis(c(
          .99, # deep_descent
          .5,  # deep_forage
          .01, # deep_ascent
          .97, # shallow_descent
          .03, # shallow_ascent
          .5   # free_surface
        )),
        depth = rep(0, 6), # slopes to ignore depth
        deep_depth = rep(0, 6), # slopes to ignore deep depth indicator
        shallow_depth = rep(0, 6) # slopes to ignore shallow depth indicator
      ),
      beta = rbind(
        intercept = log(c(
          1.6,  # deep_descent
          .03,  # deep_forage
          1.4,  # deep_ascent
          .6,   # shallow_descent
          .6,   # shallow_ascent
          .2    # free_surface
        )),
        depth = rep(0, 6), # slopes to ignore depth
        deep_depth = rep(0, 6), # slopes to ignore deep depth indicator
        shallow_depth = rep(0, 6) # slopes to ignore shallow depth indicator
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
